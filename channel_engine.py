"""Channel calculation engine. SI units internally; frequencies are strictly > 0.

Port order throughout: [left conductors..., right conductors...]. A differential
pair is two adjacent conductors P,N. S matrices use real, equal reference Z0.
See CALCULATION_GUIDE.md for approximations and data requirements.
"""
from __future__ import annotations

import copy
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6
C0 = 1 / math.sqrt(EPS0 * MU0)

# Illustrative engineering starting points, NOT manufacturer-qualified materials.
MATERIALS = {
    "FR4 example": {"dk": 4.2, "df": 0.018, "reference_ghz": 1.0},
    "Low-loss example": {"dk": 3.5, "df": 0.004, "reference_ghz": 1.0},
    "PTFE-like example": {"dk": 2.2, "df": 0.001, "reference_ghz": 1.0},
    "Custom": {"dk": 4.0, "df": 0.01, "reference_ghz": 1.0},
}


def default_line():
    return {
        "name": "PCB line", "kind": "line", "length_mm": 50.0,
        "material": "FR4 example", **MATERIALS["FR4 example"],
        "structure": "microstrip", "width_mm": 0.18, "gap_mm": 0.18,
        "lane_gap_mm": 0.8, "height_mm": 0.18, "top_height_mm": 0.18,
        "thickness_mm": 0.035, "conductivity": 5.8e7,
        "roughness_um": 0.0, "dispersion": True,
        "dielectric_low_hz": 1e3, "dielectric_high_hz": 1e14,
        "extra_lengths_mm": [],
    }


def default_via():
    return {
        "name": "Via (approximate)", "kind": "via", "model": "through",
        "length_mm": 1.0, "stub_mm": 0.5, "barrel_mm": 0.25,
        "pad_mm": 0.5, "antipad_mm": 0.8, "dk": 4.2,
        "stub_z_ohm": 60.0, "series_r_ohm": 0.02, "mutual_fraction": 0.1,
    }


def default_file(kind="connector"):
    return {
        "name": "Connector" if kind == "connector" else "Measured line",
        "kind": kind, "path": "", "port_order": [],
        "reference_length_mm": 50.0, "length_mm": 50.0,
        "uniform_line_confirmed": False,
        "includes": "Specify connector / mating pair / launches / fixture",
        "source": "User-supplied model", "replicate_pair": False,
    }


def default_project():
    p = {
        "version": 1, "signal": "differential", "lanes": 1,
        "interface": "Generic", "target_ohm": 100.0, "zref_ohm": 50.0,
        "fstart_ghz": 0.01, "fstop_ghz": 20.0, "points": 401,
        "evaluation_ghz": 1.0, "blocks": [default_line()],
    }
    p["blocks"][0]["width_mm"] = synthesize_width(p["blocks"][0], p)
    return p


@dataclass
class Network:
    f: np.ndarray
    s: np.ndarray
    z0: float = 50.0
    notes: list[str] = field(default_factory=list)

    def __post_init__(self):
        self.f = np.asarray(self.f, dtype=float)
        self.s = np.asarray(self.s, dtype=complex)
        if self.s.shape != (len(self.f), self.s.shape[-1], self.s.shape[-1]):
            raise ValueError("S array must be frequency x port x port")
        if not np.all(np.isfinite(self.s)) or not np.all(np.isfinite(self.f)):
            raise ValueError("Network contains NaN/Inf")
        if np.any(np.diff(self.f) <= 0) or np.any(self.f < 0):
            raise ValueError("Frequencies must be nonnegative and strictly increasing")
        if not np.isfinite(self.z0) or self.z0 <= 0:
            raise ValueError("Reference impedance must be positive")


def through(f, conductors, z0):
    s = np.zeros((len(f), 2 * conductors, 2 * conductors), complex)
    s[:, :conductors, conductors:] = np.eye(conductors)
    s[:, conductors:, :conductors] = np.eye(conductors)
    return Network(f, s, z0)


def cascade(a: Network, b: Network):
    """Redheffer connection: connect a.right to b.left, retain all reflections."""
    if a.s.shape != b.s.shape or not np.array_equal(a.f, b.f) or a.z0 != b.z0:
        raise ValueError("Cascade requires matching frequencies, port counts and Z0")
    n = a.s.shape[1] // 2
    a11, a12, a21, a22 = a.s[:, :n, :n], a.s[:, :n, n:], a.s[:, n:, :n], a.s[:, n:, n:]
    b11, b12, b21, b22 = b.s[:, :n, :n], b.s[:, :n, n:], b.s[:, n:, :n], b.s[:, n:, n:]
    x = np.linalg.solve(np.eye(n) - a22 @ b11, a21)
    y = np.linalg.solve(np.eye(n) - a22 @ b11, a22 @ b12)
    out = np.empty_like(a.s)
    out[:, :n, :n] = a11 + a12 @ b11 @ x
    out[:, :n, n:] = a12 @ (b12 + b11 @ y)
    out[:, n:, :n] = b21 @ x
    out[:, n:, n:] = b22 + b21 @ y
    return Network(a.f, out, a.z0, list(dict.fromkeys(a.notes + b.notes)))


def renormalize(net, znew):
    """Equal real reference renormalization without a singular Z conversion."""
    if net.z0 == znew:
        return net
    r = (znew - net.z0) / (znew + net.z0)
    ident = np.eye(net.s.shape[1])
    den, num = ident - r * net.s, net.s - r * ident
    sout = np.linalg.solve(np.swapaxes(den, -1, -2), np.swapaxes(num, -1, -2))
    return Network(net.f, np.swapaxes(sout, -1, -2), znew, net.notes.copy())


def interpolate(net, f):
    """Linear Cartesian interpolation inside available bandwidth only."""
    if f[0] < net.f[0] or f[-1] > net.f[-1]:
        raise ValueError(f"Requested band exceeds file band {net.f[0]/1e9:g}..{net.f[-1]/1e9:g} GHz")
    out = np.empty((len(f),) + net.s.shape[1:], complex)
    for i in range(out.shape[1]):
        for j in range(out.shape[2]):
            v = net.s[:, i, j]
            out[:, i, j] = np.interp(f, net.f, v.real) + 1j * np.interp(f, net.f, v.imag)
    return Network(f, out, net.z0, net.notes.copy())


def read_touchstone(path):
    """Read Touchstone 1.x full-matrix S/RI,MA,DB. Reject unsupported 2.0.

    Touchstone 2-port order is S11,S21,S12,S22; >=3 ports use row-major.
    Noise sections, unequal references and mixed-mode 2.0 are not guessed.
    """
    path = Path(path)
    match = re.search(r"\.s(\d+)p$", path.name, re.I)
    if not match:
        raise ValueError("Expected a Touchstone .sNp file")
    ports = int(match.group(1))
    if ports < 2 or ports % 2:
        raise ValueError("Channel files must have an even port count >= 2")
    unit, parameter, fmt, z0 = "GHZ", "S", "MA", 50.0
    numbers, comments = [], []
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        if "!" in line:
            comment = line.split("!", 1)[1]
            comments.append(comment)
        line = line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            raise ValueError("Touchstone 2.0 is not supported. Export full-matrix single-ended Touchstone 1.x with equal real Z0.")
        if line.startswith("#"):
            words = line[1:].upper().split()
            if len(words) < 3:
                raise ValueError("Incomplete Touchstone option line")
            unit, parameter, fmt = words[:3]
            if "R" in words:
                z0 = float(words[words.index("R") + 1])
            continue
        numbers.extend(float(v.replace("D", "E").replace("d", "e")) for v in line.split())
    if parameter != "S" or unit not in {"HZ", "KHZ", "MHZ", "GHZ"} or fmt not in {"RI", "MA", "DB"}:
        raise ValueError("Only S parameter RI/MA/DB formats are supported")
    stride = 1 + 2 * ports * ports
    if not numbers or len(numbers) % stride:
        raise ValueError("Malformed full matrix or unsupported noise section")
    a = np.array(numbers).reshape(-1, stride)
    f = a[:, 0] * {"HZ": 1, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}[unit]
    u, v = a[:, 1::2], a[:, 2::2]
    vals = u + 1j * v if fmt == "RI" else (10 ** (u / 20) if fmt == "DB" else u) * np.exp(1j * np.deg2rad(v))
    s = vals.reshape(-1, ports, ports)
    if ports == 2:
        s = s.transpose(0, 2, 1)
    notes = [f"Imported {path.name}; confirm SINGLE-ENDED port mapping and included fixtures."]
    if any("mixed" in x.lower() for x in comments):
        raise ValueError("File comments mention mixed mode. Convert to single-ended before importing.")
    return Network(f, s, z0, notes)


def port_labels(project):
    n = int(project["lanes"])
    labels = [f"L{i+1}_{p}" for i in range(n) for p in ("P", "N")] if project["signal"] == "differential" else [f"L{i+1}" for i in range(n)]
    return [f"{side}_{x}" for side in ("IN", "OUT") for x in labels]


def write_touchstone(net, path, labels=None):
    path = Path(path)
    ports = net.s.shape[1]
    path = path.with_suffix(f".s{ports}p")
    with path.open("w", encoding="utf-8") as out:
        out.write("! Channel Model Generator - single-ended, equal real references\n")
        for i, label in enumerate(labels or []):
            out.write(f"! Port {i+1}: {label}\n")
        for note in net.notes:
            out.write("! " + note.replace("\n", " ") + "\n")
        out.write(f"# Hz S RI R {net.z0:.12g}\n")
        for freq, s in zip(net.f, net.s):
            values = s.T.ravel() if ports == 2 else s.ravel()
            pairs = [f"{z.real:.12e} {z.imag:.12e}" for z in values]
            # Wrap at one matrix row. Continuation lines have no frequency.
            for row in range(ports):
                prefix = f"{freq:.12e} " if row == 0 else "  "
                out.write(prefix + " ".join(pairs[row*ports:(row+1)*ports]) + "\n")
    return path


def effective_permittivity(block, f):
    """Finite-band logarithmic Debye continuum; e^(jwt), eps=eps'-j eps''.

    Anchor real Dk and loss tangent at reference frequency; low/high poles
    bound the dispersion. Constant Dk/Df option is frequency-domain only.
    """
    dk, df = float(block["dk"]), float(block["df"])
    if dk <= 1 or df < 0:
        raise ValueError("Dk must exceed 1 and Df must be nonnegative")
    if block.get("dispersion", True) and df > 0:
        lo = float(block.get("dielectric_low_hz", 1e3))
        hi = float(block.get("dielectric_high_hz", 1e14))
        ref = float(block.get("reference_ghz", 1.0)) * 1e9
        if not 0 < lo < ref < hi:
            raise ValueError("Dielectric poles must bracket reference frequency")
        qref = np.log((hi + 1j * ref) / (lo + 1j * ref))
        strength = dk * df / (-qref.imag)
        eps_inf = dk - strength * qref.real
        if eps_inf < 1:
            raise ValueError("Dispersion inputs produce epsilon_infinity < 1")
        eps = eps_inf + strength * np.log((hi + 1j * f) / (lo + 1j * f))
    else:
        eps = np.full(len(f), dk * (1 - 1j * df), complex)
    if block["structure"] == "microstrip":
        u = block["width_mm"] / block["height_mm"]
        fill = (1 + (1 + 12 / u) ** -0.5 + (0.04 * (1-u)**2 if u < 1 else 0)) / 2
        eps = 1 + fill * (eps - 1)
    return eps


def capacitance_air(block, project, panels=6):
    """Approximate 2-D boundary-element capacitance matrix per unit length.

    Rectangular conductors are divided into straight charge panels. Point
    collocation uses analytic ground-plane image Green functions. Self panel
    log singularity is integrated through equivalent distance ds/(2e).
    This is an engineering approximation, not a qualified EM field solver.
    """
    lanes = int(project["lanes"])
    diff = project["signal"] == "differential"
    w, t, h, gap, lane_gap = [float(block[k]) * 1e-3 for k in ("width_mm", "thickness_mm", "height_mm", "gap_mm", "lane_gap_mm")]
    if min(w, t, h, gap, lane_gap) <= 0 or lanes < 1:
        raise ValueError("Geometry dimensions and lane count must be positive")
    strip = block["structure"] == "stripline"
    if block["structure"] not in {"microstrip", "stripline"}:
        raise ValueError("Unknown line structure")
    ht = float(block.get("top_height_mm", block["height_mm"])) * 1e-3
    if strip and ht <= 0:
        raise ValueError("Top dielectric height must be positive")
    group_width = 2*w + gap if diff else w
    centers = [lane*(group_width + lane_gap) + offset for lane in range(lanes) for offset in ([w/2, 1.5*w + gap] if diff else [w/2])]
    x, y, ds, owners = [], [], [], []
    for owner, center in enumerate(centers):
        # h is dielectric clearance from bottom ground to bottom copper.
        for yy in (h, h+t):
            for k in range(panels):
                x.append(center-w/2 + (k+0.5)*w/panels); y.append(yy); ds.append(w/panels); owners.append(owner)
        for xx in (center-w/2, center+w/2):
            for k in range(2):
                x.append(xx); y.append(h + (k+0.5)*t/2); ds.append(t/2); owners.append(owner)
    x, y, ds = np.array(x), np.array(y), np.array(ds)
    dx, dy = x[:, None]-x, y[:, None]-y
    rr = dx*dx + dy*dy
    np.fill_diagonal(rr, (ds/(2*math.e))**2)
    if strip:
        H = h+t+ht
        # Stable identity: cosh(a)-cos(b)=2[sinh(a/2)^2+sin(b/2)^2].
        a = np.pi * dx / H
        if np.max(np.abs(a)) > 600:
            raise ValueError("Cross section aspect ratio too large; reduce lane count/gaps")
        denom = np.sinh(a/2)**2 + np.sin(np.pi*dy/(2*H))**2
        np.fill_diagonal(denom, np.sinh(np.pi*ds/(4*math.e*H))**2)
        numer = np.sinh(a/2)**2 + np.sin(np.pi*(y[:, None]+y)/(2*H))**2
        potential = np.log(numer/denom)/(4*np.pi*EPS0)
    else:
        image_r2 = dx*dx + (y[:, None]+y)**2
        potential = np.log(image_r2/rr)/(4*np.pi*EPS0)
    incidence = np.eye(len(centers))[owners]
    cap = incidence.T @ np.linalg.solve(potential, incidence)
    cap = (cap + cap.T)/2
    if np.min(np.linalg.eigvalsh(cap)) <= 0:
        raise ValueError("Electrostatic approximation produced invalid capacitance")
    return cap


def line_modes(block, project, f):
    ca = capacitance_air(block, project)
    c, u = np.linalg.eigh(ca)
    eps = effective_permittivity(block, f)
    w = 2*np.pi*f[:, None]
    width, thick = block["width_mm"]*1e-3, block["thickness_mm"]*1e-3
    sigma = float(block.get("conductivity", 5.8e7))
    rough = float(block.get("roughness_um", 0))*1e-6
    if sigma <= 0 or rough < 0:
        raise ValueError("Invalid conductivity or roughness")
    rdc = 1/(sigma*width*thick)
    skin = np.sqrt(np.pi*f*MU0/sigma)/(2*(width+thick))
    delta = np.sqrt(1/(np.pi*f*MU0*sigma))
    rough_factor = 1 + 2/np.pi*np.arctan(1.4*(rough/delta)**2)
    # Smooth DC/skin crossover; return-plane/proximity losses are not solved.
    rac = skin*rough_factor
    resistance = np.sqrt(rdc**2 + rac**2)
    z = resistance[:, None] + 1j*w*(MU0*EPS0/c[None, :])
    y = 1j*w*c[None, :]*eps[:, None]
    gamma = np.sqrt(z*y)
    zc = np.sqrt(z/y)
    return gamma, zc, u


def characteristic_impedances(block, project, frequency=None):
    f = np.array([frequency or float(project.get("evaluation_ghz", 1))*1e9])
    _, zc, u = line_modes(block, project, f)
    zm = (u * zc[0]) @ u.T
    if project["signal"] == "differential":
        return np.array([zm[2*i, 2*i]+zm[2*i+1, 2*i+1]-2*zm[2*i, 2*i+1] for i in range(project["lanes"])])
    return np.diag(zm)


def synthesize_width(block, project):
    """Bisection: match mean real pair impedance at evaluation frequency.
    Multiple lanes can have unequal edge/interior impedances; report all values.
    """
    b = copy.deepcopy(block)
    target = float(project["target_ohm"])
    low, high = 0.005, max(2.0, 20*float(b["height_mm"]))
    def error(width):
        b["width_mm"] = width
        return float(np.mean(characteristic_impedances(b, project).real)) - target
    if error(low)*error(high) > 0:
        raise ValueError("Target impedance cannot be bracketed by 0.005 mm..maximum width")
    for _ in range(35):
        mid = (low+high)/2
        if error(mid) > 0:
            low = mid
        else:
            high = mid
    return (low+high)/2


def physical_line(block, project, f):
    length = float(block["length_mm"])*1e-3
    if length < 0:
        raise ValueError("Line length cannot be negative")
    gamma, zc, u = line_modes(block, project, f)
    z0 = float(project["zref_ohm"])
    # Stable hyperbolic equivalent using exp(-gamma*l), avoids overflow.
    p = np.exp(-gamma*length)
    rho = (zc-z0)/(zc+z0)
    den = 1-rho*rho*p*p
    reflection = rho*(1-p*p)/den
    transmission = p*(1-rho*rho)/den
    r = np.einsum("im,fm,jm->fij", u, reflection, u)
    t = np.einsum("im,fm,jm->fij", u, transmission, u)
    s = np.concatenate([np.concatenate([r, t], axis=2), np.concatenate([t, r], axis=2)], axis=1)
    notes = ["Physical line: approximate 2D charge-panel capacitance, scalar effective dielectric, quasi-TEM modal RLGC.",
             "Conductor return-plane/proximity loss, weave anisotropy, radiation and solder mask not solved; crosstalk is approximate."]
    if not block.get("dispersion", True):
        notes.append("Constant Dk/Df selected: broadband time-domain causality is not established.")
    result = Network(f, s, z0, notes)
    extras = block.get("extra_lengths_mm", [])
    if extras:
        if len(extras) != project["lanes"] or min(extras) < 0:
            raise ValueError("extra_lengths_mm must contain one nonnegative length per pair/lane")
        # Extra lengths are uncoupled tail sections, NOT a changed coupled region.
        tail = through(f, u.shape[0], z0)
        size = 2 if project["signal"] == "differential" else 1
        for lane, extra in enumerate(extras):
            bb = copy.deepcopy(block); bb["length_mm"] = extra; bb["extra_lengths_mm"] = []
            pp = {**project, "lanes": 1}
            sub = physical_line(bb, pp, f)
            ix = list(range(lane*size, (lane+1)*size))
            ix += [k+u.shape[0] for k in ix]
            tail.s[:, np.array(ix)[:, None], np.array(ix)] = sub.s
        result = cascade(result, tail)
        result.notes.append("Lane extra lengths modeled as independent tails following the common coupled section.")
    return result


def lumped_series(f, impedance, z0):
    n = impedance.shape[-1]
    ident = np.broadcast_to(np.eye(n), impedance.shape)
    t = np.linalg.solve(2*z0*ident+impedance, 2*z0*ident)
    r = ident-t
    return Network(f, np.block([[r, t], [t, r]]), z0)


def lumped_shunt(f, admittance, z0):
    n = admittance.shape[-1]
    ident = np.broadcast_to(np.eye(n), admittance.shape)
    t = np.linalg.solve(2*ident+z0*admittance, 2*ident)
    r = t-ident
    return Network(f, np.block([[r, t], [t, r]]), z0)


def via_network(block, project, f):
    """Illustrative series RL + pad shunt C + open-stub model.
    Not a 3D extraction. Ground-via/antipad fields and inter-lane XT absent.
    """
    h, d, pad, anti = [float(block[k]) for k in ("length_mm", "barrel_mm", "pad_mm", "antipad_mm")]
    stub = float(block["stub_mm"])
    dk = float(block["dk"])
    if min(h, d, pad, anti) <= 0 or anti <= pad or stub < 0 or dk <= 1:
        raise ValueError("Via needs positive dimensions, antipad > pad, nonnegative stub, Dk > 1")
    # Common approximate via inductance and capacitance formulae, inch inputs.
    L = 5.08*(h/25.4)*(math.log(4*h/d)+1)*1e-9
    C = 1.41*dk*(h/25.4)*(pad/25.4)/((anti-pad)/25.4)*1e-12
    if L <= 0:
        raise ValueError("Via aspect ratio outside lumped formula range")
    n = project["lanes"]*(2 if project["signal"] == "differential" else 1)
    lm = np.eye(n)*L
    k = float(block.get("mutual_fraction", 0.1))
    if not 0 <= k < 1:
        raise ValueError("Mutual fraction must be in [0,1)")
    if project["signal"] == "differential":
        for pair in range(project["lanes"]):
            lm[2*pair, 2*pair+1] = lm[2*pair+1, 2*pair] = k*L
    omega = 2*np.pi*f
    rr = float(block.get("series_r_ohm", 0.02))
    zs = float(block.get("stub_z_ohm", 60))
    if rr < 0 or zs <= 0:
        raise ValueError("Invalid via R/stub impedance")
    z = rr*np.eye(n)[None] + 1j*omega[:, None, None]*lm
    # Tiny attenuation regularizes an open-stub resonance, explicitly approximate.
    gamma = (0.001 + 1j*omega*np.sqrt(dk)/C0)*stub*1e-3
    ystub = np.tanh(gamma)/zs
    y = (1j*omega*C+ystub)[:, None, None]*np.eye(n)
    z0 = float(project["zref_ohm"])
    result = cascade(cascade(lumped_shunt(f, y/2, z0), lumped_series(f, z, z0)), lumped_shunt(f, y/2, z0))
    result.notes.append("Via is illustrative lumped RL/C + open stub; pair mutual L is user assumption; inter-lane via crosstalk NOT modeled. Use EM/measured data for high-speed decisions.")
    return result


def wave_chain(s):
    n = s.shape[0]//2
    a, b, c, d = s[:n, :n], s[:n, n:], s[n:, :n], s[n:, n:]
    if np.linalg.cond(c) > 1e10:
        raise ValueError("Transmission matrix too ill-conditioned for length extraction")
    ci = np.linalg.solve(c, np.eye(n))
    return np.block([[ci, -ci@d], [a@ci, b-a@ci@d]])


def chain_to_s(t):
    n = t.shape[0]//2
    a, b, c, d = t[:n, :n], t[:n, n:], t[n:, :n], t[n:, n:]
    ai = np.linalg.solve(a, np.eye(n))
    return np.block([[c@ai, d-c@ai@b], [ai, -ai@b]])


def resize_uniform(net, reference_mm, target_mm):
    """T(l_new)=exp(log(T(l_ref))*l_new/l_ref).

    Track eigenmodes by vector overlap plus eigenvalue distance, unwrap phases.
    Unique first-frequency branch is assumed. Multi-lane/degenerate modes are
    experimental; refuse poorly-conditioned eigenvectors or coarse phase steps.
    """
    if reference_mm <= 0 or target_mm < 0:
        raise ValueError("Reference length must be >0 and requested length >=0")
    if len(net.f) < 3 or net.f[0] <= 0:
        raise ValueError("Measured extraction requires >=3 positive frequency samples; omit DC")
    if target_mm == reference_mm:
        return Network(net.f, net.s.copy(), net.z0, net.notes + ["Measured reference length unchanged."])
    if target_mm == 0:
        return through(net.f, net.s.shape[1]//2, net.z0)
    vals, vecs = [], []
    previous = None
    for s in net.s:
        ev, vv = np.linalg.eig(wave_chain(s))
        if np.min(np.abs(ev)) < 1e-14 or np.linalg.cond(vv) > 1e8:
            raise ValueError("Unstable eigenmode extraction; use a shorter/de-embedded reference")
        if previous is not None:
            pv, pq = previous
            overlap = abs(pq.conj().T @ vv)
            cost = 1-overlap + 0.15*abs(pv[:, None]-ev[None, :])/np.maximum(abs(pv[:, None]), 1e-12)
            # Global greedy assignment; ambiguous eigenmodes remain a stated limit.
            order = np.full(len(ev), -1, int)
            for _ in ev:
                i, j = np.unravel_index(np.argmin(cost), cost.shape)
                order[i] = j; cost[i, :] = np.inf; cost[:, j] = np.inf
            ev, vv = ev[order], vv[:, order]
        vals.append(ev); vecs.append(vv); previous = (ev, vv)
    vals, vecs = np.array(vals), np.array(vecs)
    phases = np.unwrap(np.angle(vals), axis=0)
    if np.max(abs(phases[0])) > 0.65:
        raise ValueError("Reference begins too high for safe phase anchoring; supply lower-frequency data")
    if np.max(abs(np.diff(phases, axis=0))) > np.pi*0.8:
        raise ValueError("Frequency sampling too coarse or eigenmode tracking unstable")
    ratio = target_mm/reference_mm
    out = np.empty_like(net.s)
    rebuilt = np.empty_like(net.s)
    for k in range(len(net.f)):
        logev = np.log(abs(vals[k])) + 1j*phases[k]
        vi = np.linalg.inv(vecs[k])
        out[k] = chain_to_s((vecs[k]*np.exp(logev*ratio)) @ vi)
        rebuilt[k] = chain_to_s((vecs[k]*np.exp(logev)) @ vi)
    error = float(np.max(abs(rebuilt-net.s)))
    if error > 1e-5:
        raise ValueError(f"Reference reconstruction failed: max error {error:g}")
    if np.max(np.linalg.svd(out, compute_uv=False)) > 1.005:
        raise ValueError("Length-resized result fails passivity tolerance; extraction is unsuitable")
    notes = net.notes + [f"Uniform measured-line wave-chain scaling {reference_mm:g} -> {target_mm:g} mm; reference reconstruction error {error:.3g}.",
        "Length extraction assumes electrically short starting frequency, dense sampling and uniform de-embedded line; reconstruction alone does NOT verify these assumptions."]
    if net.s.shape[1] > 4:
        notes.append("Multi-lane measured length scaling is experimental; modal degeneracy/branch ambiguity requires independent validation.")
    return Network(net.f, out, net.z0, notes)


def imported_block(block, project, f):
    net = read_touchstone(block["path"])
    nports = net.s.shape[1]
    order = block.get("port_order", [])
    if not order:
        raise ValueError(f"{block['name']}: specify port_order, e.g. [1,2,3,4] = IN P,N; OUT P,N")
    if sorted(order) != list(range(1, nports+1)):
        raise ValueError("port_order must be a one-based permutation of every file port")
    ix = np.array(order)-1
    net = Network(net.f, net.s[:, ix][:, :, ix], net.z0, net.notes)
    net = renormalize(net, float(project["zref_ohm"]))
    if block["kind"] == "measured_line":
        if not block.get("uniform_line_confirmed", False):
            raise ValueError("Confirm uniform_line_confirmed only for uniform, de-embedded reference traces")
        # Extract on original dense samples before interpolating to output grid.
        net = resize_uniform(net, float(block["reference_length_mm"]), float(block["length_mm"]))
    net = interpolate(net, f)
    expected = len(port_labels(project))
    if nports != expected:
        unit = 4 if project["signal"] == "differential" else 2
        if nports != unit or not block.get("replicate_pair", False):
            raise ValueError(f"Expected {expected} ports, got {nports}. Explicit independent replication is optional for one lane.")
        nc = expected//2; size = unit//2
        ss = np.zeros((len(f), expected, expected), complex)
        for lane in range(project["lanes"]):
            ii = list(range(lane*size, (lane+1)*size)); ii += [x+nc for x in ii]
            ss[:, np.array(ii)[:, None], np.array(ii)] = net.s
        net = Network(f, ss, net.z0, net.notes + ["Single-lane file replicated independently: inter-lane crosstalk NOT PROVIDED."])
    net.notes.append(f"{block['name']} source: {block.get('source', 'user')}; includes: {block.get('includes', 'unspecified')}")
    return net


def mixed_mode(net):
    """Orthonormal power waves. Per side order D1,C1,D2,C2,...
    Derived physical references: differential 2*Z0, common Z0/2.
    """
    ports = net.s.shape[1]
    if ports % 4:
        raise ValueError("Mixed mode requires paired conductors on both sides")
    q = np.zeros((ports, ports))
    for i in range(0, ports, 2):
        q[i:i+2, i:i+2] = np.array([[1, -1], [1, 1]])/np.sqrt(2)
    return Network(net.f, q @ net.s @ q.T, net.z0, net.notes.copy())


def validate_project(project):
    if project.get("signal") not in {"single", "differential"}:
        raise ValueError("signal must be single or differential")
    n = project["lanes"]
    if not isinstance(n, int) or not 1 <= n <= 16:
        raise ValueError("Pair/lane count must be an integer from 1 to 16")
    if not 0 < project["fstart_ghz"] < project["fstop_ghz"]:
        raise ValueError("Need 0 < start frequency < stop frequency")
    if not isinstance(project["points"], int) or not 3 <= project["points"] <= 20001:
        raise ValueError("Frequency points must be 3..20001")
    if project["zref_ohm"] <= 0 or project["target_ohm"] <= 0 or project["evaluation_ghz"] <= 0:
        raise ValueError("Impedances and evaluation frequency must be positive")


def calculate(project, progress=None):
    validate_project(project)
    f = np.linspace(project["fstart_ghz"]*1e9, project["fstop_ghz"]*1e9, project["points"])
    nc = len(port_labels(project))//2
    total = through(f, nc, float(project["zref_ohm"]))
    details = []
    for i, b in enumerate(project["blocks"]):
        if progress:
            progress(f"Calculating {i+1}/{len(project['blocks'])}: {b['name']}")
        try:
            if b["kind"] == "line":
                net = physical_line(b, project, f)
            elif b["kind"] == "via":
                net = via_network(b, project, f)
            elif b["kind"] in {"connector", "sparameter", "measured_line"}:
                net = imported_block(b, project, f)
            else:
                raise ValueError("Unknown block kind")
            block_peak = float(np.max(np.linalg.svd(net.s, compute_uv=False)))
            if block_peak > 1.005:
                raise ValueError(f"Passive block check failed: maximum singular value {block_peak:g}")
            total = cascade(total, net)
        except Exception as exc:
            raise ValueError(f"Block {i+1} [{b.get('name')}]: {exc}") from exc
        details.append((b["name"], net))
    peak = float(np.max(np.linalg.svd(total.s, compute_uv=False)))
    total.notes.append(f"Maximum singular value = {peak:.8g}; sampled passivity check only, causality not certified.")
    if peak > 1.005:
        raise ValueError(f"Channel violates sampled passivity: maximum singular value {peak:g}")
    total.notes.append("Cross-block spatial coupling is not included; closely interacting discontinuities need a joint EM/S-parameter block.")
    return total, details


def load_project(path):
    path = Path(path)
    project = json.loads(path.read_text(encoding="utf-8-sig"))
    for block in project["blocks"]:
        if block.get("path") and not Path(block["path"]).is_absolute():
            block["path"] = str((path.parent/block["path"]).resolve())
    validate_project(project)
    return project


def save_project(project, path):
    Path(path).write_text(json.dumps(project, indent=2, ensure_ascii=False), encoding="utf-8")
