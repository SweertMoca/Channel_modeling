"""Run GUI: python Channel_model_generator.py
Headless: python Channel_model_generator.py --project project.json --output channel.s4p
Self-test: python Channel_model_generator.py --self-test

Only NumPy is required beyond the Python standard library. Tkinter is required
for the GUI. Calculation details live in channel_engine.py and CALCULATION_GUIDE.md.
"""
from __future__ import annotations

import argparse
import copy
import json
import queue
import threading
from pathlib import Path

import numpy as np

from channel_engine import (
    MATERIALS, calculate, characteristic_impedances, default_file, default_line,
    default_project, default_via, load_project, mixed_mode, port_labels,
    read_touchstone, save_project, synthesize_width, write_touchstone,
)

ROOT = Path(__file__).resolve().parent

LABELS = {
    "length_mm": "선로 길이 / via 사용 길이 [mm]", "reference_length_mm": "측정 기준 선로 길이 [mm]",
    "width_mm": "선폭 [mm]", "gap_mm": "Pair 내부 edge 간격 [mm]",
    "lane_gap_mm": "인접 lane 외곽 edge 간격 [mm]", "height_mm": "하부 ground~동박 하단 [mm]",
    "top_height_mm": "동박 상단~상부 ground [mm]", "thickness_mm": "동박 두께 [mm]",
    "dk": "Dk", "df": "Df", "reference_ghz": "Dk/Df 기준 주파수 [GHz]",
    "conductivity": "동박 전도도 [S/m]", "roughness_um": "동박 거칠기 [µm]",
    "dielectric_low_hz": "유전체 분산 하한 pole [Hz]", "dielectric_high_hz": "유전체 분산 상한 pole [Hz]",
    "dispersion": "인과적 유전체 분산 사용", "extra_lengths_mm": "Lane별 독립 추가 길이 [mm], JSON 배열",
    "stub_mm": "잔여 stub [mm]", "barrel_mm": "Via barrel 직경 [mm]", "pad_mm": "Pad 직경 [mm]",
    "antipad_mm": "Antipad 직경 [mm]", "stub_z_ohm": "Stub 가정 임피던스 [Ω]",
    "series_r_ohm": "Via 직렬 저항 [Ω]", "mutual_fraction": "Pair mutual L / self L [0..1)",
    "path": "Single-ended Touchstone 파일", "port_order": "입력→출력 포트 순서 (1부터), JSON 배열",
    "uniform_line_confirmed": "균일 선로·fixture 제거·저주파 시작 조건 확인",
    "replicate_pair": "1-lane 파일을 독립 복제 (lane 간 XT 없음)",
    "includes": "파일에 포함된 구조 / 측정 기준면", "source": "모델 출처",
    "name": "블록 이름", "structure": "단면 구조", "material": "재질 프리셋 (예시값)",
    "model": "Via 형상 프리셋",
}


class App:
    def __init__(self, tk, ttk, root):
        self.tk, self.ttk, self.root = tk, ttk, root
        self.project = default_project()
        self.selection = 0
        self.fields = {}
        self.events = queue.Queue()
        self.busy = False
        self.last = None
        self.details = []
        self.result_project = None
        self.pending = None
        root.title("Channel Model Generator | 선로 · Connector · Via")
        root.geometry("1380x920")
        root.minsize(1050, 720)
        self.status = tk.StringVar(value="초기 선폭은 목표 임피던스에 맞춰 생성했습니다. 형상을 바꾸면 임피던스만 갱신됩니다.")
        self.build()
        self.refresh_globals()
        self.refresh_blocks()
        root.after(100, self.poll)

    def build(self):
        tk, ttk = self.tk, self.ttk
        bar = ttk.Frame(self.root, padding=8); bar.pack(fill="x")
        for text, cmd in [("새 프로젝트", self.new), ("열기", self.open), ("저장", self.save),
                          ("계산", self.run), ("S-parameter 내보내기", self.export), ("계산 설명", self.help)]:
            ttk.Button(bar, text=text, command=cmd).pack(side="left", padx=3)
        gf = ttk.LabelFrame(self.root, text="채널 공통 설정 — lane = 한 방향 pair (PCIe TX/RX는 별도 pair)", padding=8)
        gf.pack(fill="x", padx=8)
        self.globals = {}
        config = [
            ("interface", "인터페이스", ["Generic", "eDP", "DP", "PCIe Gen5+"]),
            ("signal", "신호", ["differential", "single"]), ("lanes", "Pair / lane 수", None),
            ("target_ohm", "목표 Ω", None), ("zref_ohm", "포트 기준 Ω", None),
            ("fstart_ghz", "시작 GHz", None), ("fstop_ghz", "끝 GHz", None),
            ("points", "주파수 점수", None), ("evaluation_ghz", "임피던스 평가 GHz", None),
        ]
        for i, (key, label, choices) in enumerate(config):
            row, col = divmod(i, 5)
            cell = ttk.Frame(gf); cell.grid(row=row, column=col, sticky="ew", padx=6, pady=3)
            ttk.Label(cell, text=label).pack(anchor="w")
            var = tk.StringVar(); self.globals[key] = var
            widget = ttk.Combobox(cell, textvariable=var, values=choices, state="readonly", width=18) if choices else ttk.Entry(cell, textvariable=var, width=15)
            widget.pack(fill="x")
            var.trace_add("write", lambda *_: self.schedule_impedance())
            if key == "signal":
                widget.bind("<<ComboboxSelected>>", self.signal_changed)
        ttk.Label(gf, text="인터페이스는 프로젝트 태그입니다. 목표 Ω/대역은 설계 조건에 맞게 지정하세요. 규격 적합성 판정 도구가 아닙니다.").grid(row=2, column=0, columnspan=5, sticky="w")
        pane = ttk.Panedwindow(self.root, orient="horizontal"); pane.pack(fill="both", expand=True, padx=8, pady=8)
        left = ttk.Frame(pane, width=250); pane.add(left, weight=1)
        ttk.Label(left, text="채널 블록 (위 → 아래 연결)").pack(anchor="w")
        self.blocklist = tk.Listbox(left, exportselection=False, height=12)
        self.blocklist.pack(fill="both", expand=True)
        self.blocklist.bind("<<ListboxSelect>>", self.select_block)
        for text, kind in [("+ 형상 기반 선로", "line"), ("+ 측정 기반 선로", "measured_line"),
                           ("+ Connector / 파일", "connector"), ("+ Via 근사 모델", "via"),
                           ("+ Via / 기타 S-parameter", "sparameter")]:
            ttk.Button(left, text=text, command=lambda k=kind: self.add(k)).pack(fill="x", pady=2)
        self.ttk.Button(left, text="+ 저장된 부품 프리셋", command=self.choose_preset).pack(fill="x", pady=2)
        self.ttk.Button(left, text="현재 파일 블록을 프리셋 저장", command=self.register_preset).pack(fill="x", pady=2)
        for text, cmd in [("위로", lambda: self.move(-1)), ("아래로", lambda: self.move(1)),
                          ("복제", self.duplicate), ("삭제", self.delete)]:
            ttk.Button(left, text=text, command=cmd).pack(fill="x", pady=1)
        middle = ttk.Frame(pane); pane.add(middle, weight=3)
        self.impedance = tk.StringVar(value="")
        ttk.Label(middle, textvariable=self.impedance, wraplength=490, foreground="#135783").pack(fill="x", pady=5)
        controls = ttk.Frame(middle); controls.pack(fill="x")
        ttk.Button(controls, text="목표 임피던스로 재설계", command=self.match).pack(side="left")
        ttk.Button(controls, text="블록 설정 적용", command=self.apply).pack(side="left", padx=4)
        canvas = tk.Canvas(middle, highlightthickness=0)
        scroll = ttk.Scrollbar(middle, orient="vertical", command=canvas.yview)
        scroll.pack(side="right", fill="y"); canvas.pack(fill="both", expand=True)
        canvas.configure(yscrollcommand=scroll.set)
        self.form = ttk.Frame(canvas, padding=6)
        win = canvas.create_window((0, 0), window=self.form, anchor="nw")
        self.form.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(win, width=e.width))
        right = ttk.Frame(pane); pane.add(right, weight=4)
        self.tabs = ttk.Notebook(right); self.tabs.pack(fill="both", expand=True)
        graph = ttk.Frame(self.tabs); self.tabs.add(graph, text="결과 그래프")
        select = ttk.Frame(graph); select.pack(fill="x")
        self.metric = tk.StringVar(value="Insertion loss")
        metric = ttk.Combobox(select, textvariable=self.metric, values=["Insertion loss", "Return loss", "NEXT", "FEXT", "D to C", "Group delay"], state="readonly", width=18)
        metric.pack(side="left"); metric.bind("<<ComboboxSelected>>", lambda e: self.plot())
        self.result_source = tk.StringVar(value="Entire channel")
        self.source_box = ttk.Combobox(select, textvariable=self.result_source, values=["Entire channel"], state="readonly", width=22)
        self.source_box.pack(side="left", padx=5); self.source_box.bind("<<ComboboxSelected>>", lambda e: self.plot())
        self.plotcanvas = tk.Canvas(graph, background="#ffffff", highlightthickness=0)
        self.plotcanvas.pack(fill="both", expand=True); self.plotcanvas.bind("<Configure>", lambda e: self.plot())
        logtab = ttk.Frame(self.tabs); self.tabs.add(logtab, text="모델 가정 / 계산 로그")
        self.log = tk.Text(logtab, wrap="word", font=("Consolas", 10)); self.log.pack(fill="both", expand=True)
        ttk.Label(self.root, textvariable=self.status, padding=6, wraplength=1300).pack(fill="x")

    def error(self, exc):
        from tkinter import messagebox
        messagebox.showerror("입력 / 계산 오류", str(exc), parent=self.root)
        self.status.set(str(exc))

    def refresh_globals(self):
        for key, var in self.globals.items():
            var.set(str(self.project[key]))

    def collect_globals(self):
        p = copy.deepcopy(self.project)
        for key, var in self.globals.items():
            value = var.get()
            p[key] = int(value) if key in {"lanes", "points"} else value if key in {"signal", "interface"} else float(value)
        return p

    def collect_block(self):
        if not self.project["blocks"]:
            return None
        b = copy.deepcopy(self.project["blocks"][self.selection])
        for key, (var, old) in self.fields.items():
            value = var.get()
            if isinstance(old, bool):
                b[key] = bool(value)
            elif isinstance(old, list):
                b[key] = json.loads(value)
                if not isinstance(b[key], list):
                    raise ValueError(key + " must be a JSON array")
            elif isinstance(old, (int, float)):
                b[key] = float(value)
            else:
                b[key] = value
        return b

    def commit(self):
        p = self.collect_globals()
        b = self.collect_block()
        if b is not None:
            p["blocks"][self.selection] = b
        self.project = p
        self.last = None
        self.plotcanvas.delete("all")

    def refresh_blocks(self):
        self.blocklist.delete(0, "end")
        for i, block in enumerate(self.project["blocks"]):
            self.blocklist.insert("end", f"{i+1}. {block['name']} [{block['kind']}]")
        if self.project["blocks"]:
            self.selection = max(0, min(self.selection, len(self.project["blocks"])-1))
            self.blocklist.selection_set(self.selection)
        self.render_form()

    def render_form(self):
        for child in self.form.winfo_children():
            child.destroy()
        self.fields = {}
        if not self.project["blocks"]:
            return
        b = self.project["blocks"][self.selection]
        for row, (key, value) in enumerate((k, v) for k, v in b.items() if k != "kind"):
            self.ttk.Label(self.form, text=LABELS.get(key, key)).grid(row=2*row, column=0, sticky="w", pady=(6, 0))
            if isinstance(value, bool):
                var = self.tk.BooleanVar(value=value)
                widget = self.ttk.Checkbutton(self.form, variable=var)
            else:
                var = self.tk.StringVar(value=json.dumps(value) if isinstance(value, list) else str(value))
                choices = list(MATERIALS) if key == "material" else ["microstrip", "stripline"] if key == "structure" else ["through", "backdrilled", "blind", "microvia"] if key == "model" else None
                widget = self.ttk.Combobox(self.form, textvariable=var, values=choices, state="readonly") if choices else self.ttk.Entry(self.form, textvariable=var)
                if key == "material":
                    widget.bind("<<ComboboxSelected>>", self.material_changed)
                if key == "model":
                    widget.bind("<<ComboboxSelected>>", self.via_changed)
            widget.grid(row=2*row+1, column=0, sticky="ew")
            self.fields[key] = (var, value)
            var.trace_add("write", lambda *_: self.schedule_impedance())
            if key == "path":
                self.ttk.Button(self.form, text="파일 선택 / 포트 확인", command=self.browse_model).grid(row=2*row+1, column=1, padx=3)
        self.form.columnconfigure(0, weight=1)
        self.schedule_impedance()

    def schedule_impedance(self):
        if self.pending:
            self.root.after_cancel(self.pending)
        self.pending = self.root.after(450, self.update_impedance)
        if hasattr(self, "last") and self.last is not None:
            self.status.set("입력값이 변경되었습니다. 내보내기 전에 다시 계산하세요.")

    def update_impedance(self):
        self.pending = None
        try:
            b, p = self.collect_block(), self.collect_globals()
            if not b:
                self.impedance.set(""); return
            if b["kind"] != "line":
                self.impedance.set("파일 모델은 형상/임피던스를 자동 변경하지 않습니다." if b["kind"] != "via" else "Via: 초기 설계용 근사 모델. Gen5+ 검증에는 측정/EM 파일을 사용하세요.")
                return
            if not 1 <= p["lanes"] <= 16:
                raise ValueError("Lane count: 1..16")
            z = characteristic_impedances(b, p)
            values = ", ".join(f"L{i+1}: {v.real:.2f}{v.imag:+.2f}j Ω" for i, v in enumerate(z))
            self.impedance.set(f"목표 {p['target_ohm']:g} Ω | {p['evaluation_ghz']:g} GHz\n{values}\n형상 유지 · 2D 근사 계산값")
        except Exception as exc:
            self.impedance.set("임피던스 입력 확인: " + str(exc))

    def material_changed(self, _=None):
        material = self.fields["material"][0].get()
        for key, value in MATERIALS[material].items():
            self.fields[key][0].set(str(value))

    def signal_changed(self, _=None):
        self.globals["target_ohm"].set("100" if self.globals["signal"].get() == "differential" else "50")
        if self.project["blocks"] and self.project["blocks"][self.selection]["kind"] == "line":
            self.match()

    def register_preset(self):
        from tkinter import simpledialog
        try:
            self.commit()
            block = self.project["blocks"][self.selection]
            if block["kind"] not in {"connector", "sparameter"}:
                raise ValueError("Connector / Via S-parameter 파일 블록을 선택하세요.")
            net = read_touchstone(block["path"])
            if sorted(block.get("port_order", [])) != list(range(1, net.s.shape[1]+1)):
                raise ValueError("포트 매핑을 먼저 지정하세요.")
            name = simpledialog.askstring("로컬 프리셋 저장", "부품명 / 품번 / 포함 구조:", initialvalue=block["name"])
            if not name:
                return
            library = ROOT/"model_library"; library.mkdir(exist_ok=True)
            catalog_path = library/"user_presets.json"
            catalog = json.loads(catalog_path.read_text(encoding="utf-8")) if catalog_path.exists() else []
            entry = copy.deepcopy(block); entry["name"] = name
            entry["path"] = str(Path(entry["path"]).resolve())
            catalog.append(entry)
            catalog_path.write_text(json.dumps(catalog, indent=2, ensure_ascii=False), encoding="utf-8")
            self.status.set("로컬 프리셋 등록 완료. 원본 S-parameter 파일 경로를 참조합니다.")
        except Exception as exc:
            self.error(exc)

    def choose_preset(self):
        from tkinter import messagebox
        try:
            path = ROOT/"model_library"/"user_presets.json"
            catalog = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
            if not catalog:
                messagebox.showinfo("프리셋 라이브러리", "검증/배포 허가된 제조사 Touchstone 파일은 아직 포함되어 있지 않습니다.\n\nConnector 또는 Via 파일 블록을 추가한 후,\n'현재 파일 블록을 프리셋 저장'으로 등록하세요.\n\nmodel_library/README.md에 제조사 자료 경로를 안내했습니다.")
                return
            dialog = self.tk.Toplevel(self.root); dialog.title("저장된 부품 프리셋")
            box = self.tk.Listbox(dialog, width=85, height=15); box.pack(fill="both", expand=True)
            for entry in catalog:
                box.insert("end", entry["name"] + (" [파일 없음]" if not Path(entry["path"]).exists() else ""))
            def select():
                try:
                    picked = box.curselection()
                    if not picked:
                        return
                    entry = copy.deepcopy(catalog[picked[0]])
                    read_touchstone(entry["path"])
                    self.commit()
                    index = self.selection+1 if self.project["blocks"] else 0
                    self.project["blocks"].insert(index, entry); self.selection = index
                    self.refresh_blocks(); dialog.destroy()
                except Exception as exc:
                    self.error(exc)
            self.ttk.Button(dialog, text="채널에 추가", command=select).pack(pady=5)
        except Exception as exc:
            self.error(exc)

    def via_changed(self, _=None):
        model = self.fields["model"][0].get()
        examples = {"through": (1.0, 0.5, 0.25, 0.5, 0.8), "backdrilled": (1.0, 0.1, 0.25, 0.5, 0.8),
                    "blind": (0.3, 0.0, 0.15, 0.35, 0.6), "microvia": (0.1, 0.0, 0.08, 0.2, 0.4)}
        for key, value in zip(("length_mm", "stub_mm", "barrel_mm", "pad_mm", "antipad_mm"), examples[model]):
            self.fields[key][0].set(str(value))

    def select_block(self, _=None):
        indices = self.blocklist.curselection()
        if not indices or indices[0] == self.selection:
            return
        try:
            self.commit(); self.selection = indices[0]; self.render_form()
        except Exception as exc:
            self.blocklist.selection_clear(0, "end"); self.blocklist.selection_set(self.selection); self.error(exc)

    def apply(self):
        try:
            self.commit(); self.refresh_blocks(); self.status.set("설정을 적용했습니다. 계산 버튼으로 S-parameter를 갱신하세요.")
        except Exception as exc:
            self.error(exc)

    def add(self, kind):
        try:
            self.commit()
            b = default_line() if kind == "line" else default_via() if kind == "via" else default_file(kind)
            if kind == "sparameter":
                b["name"] = "Via / other measured block"
            if kind == "line":
                b["width_mm"] = synthesize_width(b, self.project)
            index = self.selection+1 if self.project["blocks"] else 0
            self.project["blocks"].insert(index, b); self.selection = index; self.refresh_blocks()
        except Exception as exc:
            self.error(exc)

    def move(self, step):
        try:
            self.commit(); other = self.selection+step
            if 0 <= other < len(self.project["blocks"]):
                bs = self.project["blocks"]; bs[self.selection], bs[other] = bs[other], bs[self.selection]
                self.selection = other; self.refresh_blocks()
        except Exception as exc:
            self.error(exc)

    def duplicate(self):
        try:
            self.commit()
            if self.project["blocks"]:
                b = copy.deepcopy(self.project["blocks"][self.selection]); b["name"] += " copy"
                self.project["blocks"].insert(self.selection+1, b); self.selection += 1; self.refresh_blocks()
        except Exception as exc:
            self.error(exc)

    def delete(self):
        try:
            self.commit()
            if self.project["blocks"]:
                del self.project["blocks"][self.selection]; self.refresh_blocks()
        except Exception as exc:
            self.error(exc)

    def match(self):
        try:
            self.commit()
            b = self.project["blocks"][self.selection]
            if b["kind"] != "line":
                raise ValueError("형상 기반 선로에서만 재설계할 수 있습니다.")
            b["width_mm"] = synthesize_width(b, self.project)
            self.render_form(); self.status.set("평가 주파수에서 lane 평균 Re(Z)를 목표에 맞췄습니다. Lane별 편차는 상단에 표시됩니다.")
        except Exception as exc:
            self.error(exc)

    def browse_model(self):
        from tkinter import filedialog, messagebox
        path = filedialog.askopenfilename(filetypes=[("Touchstone", "*.s*p"), ("All", "*.*")])
        if not path:
            return
        try:
            net = read_touchstone(path)
            self.fields["path"][0].set(path)
            self.fields["port_order"][0].set(json.dumps(list(range(1, net.s.shape[1]+1))))
            messagebox.showinfo("포트 매핑 확인 필요", f"{net.s.shape[1]} ports, {net.z0:g} Ω reference\n{net.f[0]/1e9:g}..{net.f[-1]/1e9:g} GHz\n\n제안한 포트 순서를 반드시 확인하세요.\n[입력 L1_P,L1_N,L2_P,L2_N,..., 출력 ...]\n파일의 물리 포트 번호를 이 순서로 배열합니다.\n측정 기준 길이는 파일에서 자동 추론하지 않습니다.")
        except Exception as exc:
            self.error(exc)

    def new(self):
        self.project = default_project(); self.selection = 0; self.last = None
        self.refresh_globals(); self.refresh_blocks(); self.plot()

    def open(self):
        from tkinter import filedialog
        path = filedialog.askopenfilename(initialdir=ROOT, filetypes=[("Project JSON", "*.json")])
        if path:
            try:
                self.project = load_project(path); self.selection = 0; self.last = None
                self.refresh_globals(); self.refresh_blocks(); self.plot()
            except Exception as exc:
                self.error(exc)

    def save(self):
        from tkinter import filedialog
        try:
            self.commit()
            path = filedialog.asksaveasfilename(initialdir=ROOT, defaultextension=".json", filetypes=[("Project JSON", "*.json")])
            if path:
                save_project(self.project, path); self.status.set("저장: " + path)
        except Exception as exc:
            self.error(exc)

    def run(self):
        if self.busy:
            return
        try:
            self.commit(); snapshot = copy.deepcopy(self.project)
            self.busy = True; self.status.set("계산 중...")
            def worker():
                try:
                    net, details = calculate(snapshot, lambda s: self.events.put(("progress", s)))
                    self.events.put(("done", (net, details, snapshot)))
                except Exception as exc:
                    self.events.put(("error", str(exc)))
            threading.Thread(target=worker, daemon=True).start()
        except Exception as exc:
            self.error(exc)

    def poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "progress":
                    self.status.set(value)
                elif kind == "error":
                    self.busy = False; self.error(value)
                else:
                    self.busy = False
                    self.last, self.details, self.result_project = value
                    self.source_box.configure(values=["Entire channel"]+[f"{i+1}: {name}" for i, (name, _) in enumerate(self.details)])
                    self.result_source.set("Entire channel")
                    self.log.delete("1.0", "end")
                    self.log.insert("end", "\n\n".join(self.last.notes))
                    self.status.set(f"계산 완료: {self.last.s.shape[1]} ports, {len(self.last.f)} frequencies. 모델 가정 탭을 확인하세요.")
                    self.plot()
        except queue.Empty:
            pass
        self.root.after(100, self.poll)

    def plot(self):
        c = self.plotcanvas; c.delete("all")
        if self.last is None:
            c.create_text(20, 25, text="계산 후 S-parameter 응답을 표시합니다.", anchor="nw", fill="#536273"); return
        net = self.last
        if self.result_source.get() != "Entire channel":
            try:
                net = self.details[int(self.result_source.get().split(":")[0])-1][1]
            except (ValueError, IndexError):
                pass
        p = self.result_project
        diff = p["signal"] == "differential"
        s = mixed_mode(net).s if diff else net.s
        n = s.shape[1]//2; stride = 2 if diff else 1
        metric = self.metric.get(); traces = []
        for lane in range(p["lanes"]):
            k = lane*stride
            if metric == "Insertion loss":
                traces.append((f"L{lane+1}", s[:, n+k, k]))
            elif metric == "Return loss":
                traces.append((f"L{lane+1}", s[:, k, k]))
            elif metric == "D to C" and diff:
                traces.append((f"L{lane+1} D→C", s[:, n+k+1, k]))
            elif metric in {"NEXT", "FEXT"} and lane > 0:
                traces.append((f"L1→L{lane+1}", s[:, (n if metric == "FEXT" else 0)+k, 0]))
            elif metric == "Group delay":
                traces.append((f"L{lane+1}", s[:, n+k, k]))
        if not traces:
            c.create_text(20, 25, text="이 지표는 현재 신호/lane 구성에서 사용할 수 없습니다.", anchor="nw"); return
        x = net.f/1e9
        ys = [(name, -np.gradient(np.unwrap(np.angle(v)), 2*np.pi*net.f)*1e12 if metric == "Group delay" else 20*np.log10(np.maximum(abs(v), 1e-8))) for name, v in traces]
        ymin = min(float(np.min(y)) for _, y in ys); ymax = max(float(np.max(y)) for _, y in ys)
        if metric != "Group delay":
            ymin = max(-160, math_floor(ymin/10)*10); ymax = max(0, ymax)
        span = max(ymax-ymin, 1); ymin -= span*.05; ymax += span*.05
        W, H = max(c.winfo_width(), 250), max(c.winfo_height(), 220)
        left, right, top, bottom = 65, W-20, 45, H-90
        for k in range(6):
            xx = left+(right-left)*k/5; yy = top+(bottom-top)*k/5
            c.create_line(xx, top, xx, bottom, fill="#e8edf2")
            c.create_line(left, yy, right, yy, fill="#e8edf2")
            c.create_text(xx, bottom+15, text=f"{x[0]+(x[-1]-x[0])*k/5:.2g}")
            c.create_text(left-7, yy, text=f"{ymax-(ymax-ymin)*k/5:.3g}", anchor="e")
        c.create_text(left, 18, text=metric+ (" [ps]" if metric == "Group delay" else " — 20 log10 |S| [dB]"), anchor="w")
        c.create_text((left+right)/2, bottom+35, text="Frequency [GHz]")
        colors = ["#1267a5", "#d85b27", "#21835b", "#8c50ac", "#bc355b", "#98721c"]
        for i, (name, y) in enumerate(ys):
            xx = left+(x-x[0])/(x[-1]-x[0])*(right-left)
            yy = top+(ymax-y)/(ymax-ymin)*(bottom-top)
            coords = np.column_stack([xx, yy]).ravel().tolist()
            c.create_line(*coords, fill=colors[i%len(colors)], width=2)
            c.create_text(left+(i%4)*100, bottom+55+(i//4)*15, text=name, anchor="w", fill=colors[i%len(colors)])

    def export(self):
        from tkinter import filedialog
        try:
            if self.last is None:
                raise ValueError("먼저 계산하세요.")
            current = self.collect_globals()
            b = self.collect_block()
            if b is not None:
                current["blocks"][self.selection] = b
            if current != self.result_project:
                raise ValueError("계산 후 설정이 변경되었습니다. 다시 계산하세요.")
            suffix = f".s{self.last.s.shape[1]}p"
            path = filedialog.asksaveasfilename(initialdir=ROOT, defaultextension=suffix, filetypes=[("Touchstone", "*"+suffix)])
            if path:
                path = export_bundle(self.last, self.result_project, path)
                self.status.set("저장 완료: " + str(path) + " + project/report JSON")
        except Exception as exc:
            self.error(exc)

    def help(self):
        import webbrowser
        webbrowser.open((ROOT/"CALCULATION_GUIDE.md").as_uri())


def math_floor(value):
    return int(np.floor(value))


def export_bundle(net, project, path):
    path = write_touchstone(net, path, port_labels(project))
    save_project(project, path.with_suffix(".project.json"))
    report = {"ports": port_labels(project), "zref_single_ended_ohm": net.z0,
              "minimum_frequency_hz": float(net.f[0]), "maximum_frequency_hz": float(net.f[-1]),
              "max_singular_value": float(np.max(np.linalg.svd(net.s, compute_uv=False))),
              "notes": net.notes}
    path.with_suffix(".report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        import unittest
        suite = unittest.defaultTestLoader.discover(str(ROOT), pattern="test_channel_engine.py")
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        raise SystemExit(0 if result.wasSuccessful() else 1)
    if args.output:
        p = load_project(args.project) if args.project else default_project()
        net, _ = calculate(p, print)
        print(export_bundle(net, p, args.output))
        return
    import tkinter as tk
    from tkinter import ttk
    root = tk.Tk()
    app = App(tk, ttk, root)
    if args.project:
        app.project = load_project(args.project); app.refresh_globals(); app.refresh_blocks()
    root.mainloop()


if __name__ == "__main__":
    main()
