"""Numerical regression/physics checks; run with --self-test or unittest."""
import copy
import tempfile
import unittest
from pathlib import Path

import numpy as np

import channel_engine as e


class ChannelTests(unittest.TestCase):
    def setUp(self):
        self.p = e.default_project()
        self.p['points'] = 151
        self.f = np.linspace(1e6, 20e9, 301)
        self.b = self.p['blocks'][0]

    def test_initial_matching_and_width_response(self):
        z = e.characteristic_impedances(self.b, self.p).real[0]
        self.assertAlmostEqual(z, 100, places=5)
        b = {**self.b, 'width_mm': self.b['width_mm']*1.2}
        self.assertLess(e.characteristic_impedances(b, self.p).real[0], z)
        self.p['signal'] = 'single'; self.p['target_ohm'] = 50
        self.b['width_mm'] = e.synthesize_width(self.b, self.p)
        self.assertAlmostEqual(e.characteristic_impedances(self.b, self.p).real[0], 50, places=5)

    def test_zero_length(self):
        net = e.physical_line({**self.b, 'length_mm': 0}, self.p, self.f)
        np.testing.assert_allclose(net.s, e.through(self.f, 2, 50).s, atol=1e-12)

    def test_cascade_equal_total_length(self):
        a = e.physical_line({**self.b, 'length_mm': 27}, self.p, self.f)
        b = e.physical_line({**self.b, 'length_mm': 43}, self.p, self.f)
        full = e.physical_line({**self.b, 'length_mm': 70}, self.p, self.f)
        np.testing.assert_allclose(e.cascade(a, b).s, full.s, atol=1e-11)

    def test_measured_length_conversion(self):
        ref = e.physical_line({**self.b, 'length_mm': 40}, self.p, self.f)
        actual = e.physical_line({**self.b, 'length_mm': 67}, self.p, self.f)
        resized = e.resize_uniform(ref, 40, 67)
        np.testing.assert_allclose(resized.s, actual.s, atol=1e-9)

    def test_multilane_length_and_passivity(self):
        self.p['lanes'] = 2
        self.b['lane_gap_mm'] = 0.25
        net = e.physical_line(self.b, self.p, self.f)
        self.assertEqual(net.s.shape, (301, 8, 8))
        self.assertLessEqual(np.max(np.linalg.svd(net.s, compute_uv=False)), 1+1e-10)
        np.testing.assert_allclose(net.s, net.s.transpose(0, 2, 1), atol=1e-12)
        self.assertGreater(np.max(abs(e.mixed_mode(net).s[:, 6, 0])), 1e-6)
        expected = e.physical_line({**self.b, 'length_mm': 80}, self.p, self.f)
        np.testing.assert_allclose(e.resize_uniform(net, 50, 80).s, expected.s, atol=1e-8)

    def test_stripline_and_material_presets(self):
        for name, material in e.MATERIALS.items():
            b = {**self.b, **material, 'material': name, 'structure': 'stripline'}
            b['width_mm'] = e.synthesize_width(b, self.p)
            self.assertAlmostEqual(e.characteristic_impedances(b, self.p).real[0], 100, places=5)
            net = e.physical_line(b, self.p, self.f)
            self.assertLessEqual(np.max(np.linalg.svd(net.s, compute_uv=False)), 1+1e-10)

    def test_reference_renormalization_roundtrip(self):
        a = e.physical_line(self.b, self.p, self.f)
        b = e.renormalize(e.renormalize(a, 75), 50)
        np.testing.assert_allclose(a.s, b.s, atol=1e-12)
        np.testing.assert_allclose(e.renormalize(e.through(self.f, 2, 50), 75).s, e.through(self.f, 2, 75).s, atol=1e-12)

    def test_touchstone_roundtrip_and_external_order(self):
        rng = np.random.default_rng(9)
        with tempfile.TemporaryDirectory() as tmp:
            for ports in [2, 4, 8]:
                s = .05*(rng.normal(size=(3, ports, ports))+1j*rng.normal(size=(3, ports, ports)))
                original = e.Network([1e6, 2e6, 3e6], s, 50)
                path = e.write_touchstone(original, Path(tmp)/f'test.s{ports}p')
                np.testing.assert_allclose(e.read_touchstone(path).s, s, atol=1e-12)
            # Independent hand-authored asymmetric files detect matching writer/reader bugs.
            path = Path(tmp)/'hand.s2p'
            path.write_text('# MHz S RI R 50\n1 1 0 2 0 3 0 4 0\n')
            np.testing.assert_array_equal(e.read_touchstone(path).s[0], [[1,3],[2,4]])
            path = Path(tmp)/'hand.s4p'
            path.write_text('# MHz S RI R 50\n1 '+' '.join(f'{i} 0' for i in range(16)))
            np.testing.assert_array_equal(e.read_touchstone(path).s[0], np.arange(16).reshape(4,4))
            path = Path(tmp)/'db.s2p'
            path.write_text('# GHz S DB R 75\n1 -20 0 -6 90 -6 -90 -20 180\n')
            a = e.read_touchstone(path)
            self.assertEqual(a.z0, 75)
            self.assertAlmostEqual(a.s[0,1,0].imag, 10**(-6/20))

    def test_via_and_full_chain_passivity(self):
        self.p['lanes'] = 2
        self.p['blocks'] = [copy.deepcopy(self.b), e.default_via(), copy.deepcopy(self.b)]
        n, _ = e.calculate(self.p)
        self.assertLessEqual(np.max(np.linalg.svd(n.s, compute_uv=False)), 1+1e-10)
        np.testing.assert_allclose(n.s, n.s.transpose(0,2,1), atol=1e-10)

    def test_import_mapping_and_replication(self):
        with tempfile.TemporaryDirectory() as tmp:
            net = e.physical_line(self.b, self.p, self.f)
            perm = [2,0,3,1]
            shuffled = e.Network(self.f, net.s[:,perm][:,:,perm], 50)
            path = e.write_touchstone(shuffled, Path(tmp)/'shuffled.s4p')
            block = e.default_file(); block['path'] = str(path)
            block['port_order'] = (np.argsort(perm)+1).tolist()
            np.testing.assert_allclose(e.imported_block(block, self.p, self.f).s, net.s, atol=1e-11)
            self.p['lanes'] = 2
            with self.assertRaises(ValueError):
                e.imported_block(block, self.p, self.f)
            block['replicate_pair'] = True
            replicated = e.imported_block(block, self.p, self.f)
            np.testing.assert_array_equal(replicated.s[:,2,0], 0)

    def test_extra_lengths(self):
        self.p['lanes'] = 2
        net = e.physical_line({**self.b, 'extra_lengths_mm': [0, 10]}, self.p, self.f)
        self.assertLessEqual(np.max(np.linalg.svd(net.s, compute_uv=False)), 1+1e-10)

    def test_invalid_inputs_and_no_extrapolation(self):
        net = e.physical_line(self.b, self.p, self.f)
        with self.assertRaises(ValueError):
            e.interpolate(net, np.array([1, 21e9]))
        with self.assertRaises(ValueError):
            e.resize_uniform(net, 0, 10)
        with self.assertRaises(ValueError):
            e.via_network({**e.default_via(), 'antipad_mm': .1}, self.p, self.f)

    def test_independent_matched_delay_oracle(self):
        # Independent ideal 2-port with known delay; not generated by physical_line.
        delay = 0.4e-9
        ref = e.through(self.f, 1, 50)
        wave = np.exp(-1j*2*np.pi*self.f*delay)
        ref.s[:,0,1] = wave; ref.s[:,1,0] = wave
        resized = e.resize_uniform(ref, 50, 125)
        np.testing.assert_allclose(resized.s[:,1,0], np.exp(-1j*2*np.pi*self.f*delay*2.5), atol=1e-11)
        np.testing.assert_allclose(resized.s[:,0,0], 0, atol=1e-12)

    def test_active_input_rejected_even_with_lossy_trace(self):
        with tempfile.TemporaryDirectory() as tmp:
            active = e.through(self.f, 2, 50)
            active.s *= 1.1
            path = e.write_touchstone(active, Path(tmp)/'active.s4p')
            block = e.default_file()
            block.update(path=str(path), port_order=[1,2,3,4])
            self.p['blocks'] = [block, {**self.b, 'length_mm': 1000}]
            with self.assertRaisesRegex(ValueError, 'Passive block check failed'):
                e.calculate(self.p)


if __name__ == '__main__':
    unittest.main(verbosity=2)
