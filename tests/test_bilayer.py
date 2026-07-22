from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.io import read

from extension_bsf import (
    BilayerError,
    build_bilayer,
    generate_scan,
    periodic_shifts,
    read_monolayer,
    write_bilayer,
)
from extension_bsf.cli import main
import extension_bsf.bilayer as bilayer_module


FIXTURE = Path(__file__).parent / "fixtures" / "scaled_cartesian_selective.POSCAR"
BOUNDARY_FIXTURE = (
    Path(__file__).parent / "fixtures" / "negative_scale_direct_boundary.POSCAR"
)


class BilayerTests(unittest.TestCase):
    def test_ase_honors_scale_cartesian_and_selective_dynamics(self) -> None:
        atoms = read_monolayer(FIXTURE)
        np.testing.assert_allclose(atoms.cell.lengths(), [2.0, 2.0, 20.0])
        np.testing.assert_allclose(atoms.positions[:, 2], [4.0, 5.0])
        self.assertTrue(atoms.constraints)

    def test_build_has_requested_envelope_gap_and_vacuum(self) -> None:
        atoms = read_monolayer(FIXTURE)
        result = build_bilayer(
            atoms, interlayer_distance=3.2, vacuum=20.0, shift=(1.5, -1.0)
        )
        self.assertEqual(len(result.atoms), 4)
        self.assertEqual(result.atoms.constraints[0].get_indices().tolist(), [1, 3])
        self.assertGreaterEqual(result.minimum_interlayer_distance, 3.2 - 1.0e-8)
        self.assertAlmostEqual(result.atoms.cell.lengths()[2], 25.2)
        np.testing.assert_array_equal(result.atoms.arrays["layer_id"], [0, 0, 1, 1])

    def test_negative_scale_direct_layer_is_unwrapped_across_z(self) -> None:
        atoms = read_monolayer(BOUNDARY_FIXTURE)
        np.testing.assert_allclose(atoms.cell.lengths(), [2.0, 2.0, 20.0])
        result = build_bilayer(atoms, interlayer_distance=3.0, vacuum=20.0)
        self.assertAlmostEqual(result.atoms.cell.lengths()[2], 27.0)
        self.assertGreaterEqual(result.internal_interlayer_distance, 3.0 - 1.0e-10)

    def test_ambiguous_periodic_unwrap_fails_closed(self) -> None:
        atoms = Atoms(
            "CC",
            scaled_positions=[[0.0, 0.0, 0.25], [0.5, 0.5, 0.75]],
            cell=[2.0, 2.0, 20.0],
            pbc=True,
        )
        with self.assertRaisesRegex(BilayerError, "cannot uniquely unwrap"):
            build_bilayer(atoms, interlayer_distance=3.0, vacuum=20.0)

    def test_collision_limit_fails_closed(self) -> None:
        atoms = Atoms("C", positions=[[0.0, 0.0, 4.0]], cell=[2.0, 2.0, 12.0], pbc=True)
        with self.assertRaisesRegex(BilayerError, "full-periodic cross-layer"):
            build_bilayer(
                atoms,
                interlayer_distance=0.5,
                vacuum=10.0,
                minimum_distance=1.0,
            )

    def test_sheared_vacuum_fails_closed(self) -> None:
        atoms = Atoms(
            "C",
            positions=[[0.0, 0.0, 0.0]],
            cell=[[2, 0, 0], [0, 2, 0], [1, 0, 10]],
            pbc=True,
        )
        with self.assertRaisesRegex(BilayerError, "sheared-vacuum"):
            build_bilayer(atoms, interlayer_distance=3.0, vacuum=10.0)

    def test_periodic_grid_has_no_duplicate_endpoint(self) -> None:
        shifts = periodic_shifts(9, 5)
        self.assertEqual(len(shifts), 45)
        self.assertEqual(len(set(shifts)), 45)
        self.assertTrue(all(0.0 <= x < 1.0 and 0.0 <= y < 1.0 for x, y in shifts))
        with self.assertRaisesRegex(BilayerError, "positive integer"):
            periodic_shifts(2.5, 3)  # type: ignore[arg-type]

    def test_write_round_trip_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "POSCAR"
            manifest = write_bilayer(
                FIXTURE,
                output,
                interlayer_distance=3.2,
                vacuum=20.0,
                shift=(0.5, 0.25),
            )
            parsed = read(output, format="vasp")
            self.assertEqual(len(parsed), 4)
            self.assertEqual(parsed.constraints[0].get_indices().tolist(), [1, 3])
            self.assertEqual(manifest["validation"]["total_atom_count"], 4)
            self.assertEqual(manifest["parameters"]["max_cross_layer_pairs"], 4_000_000)
            self.assertEqual(
                manifest["input"]["sha256"],
                bilayer_module._digest_bytes(FIXTURE.read_bytes()),
            )
            self.assertEqual(len(manifest["validation"]["layer_id_by_output_index"]), 4)
            stored = json.loads((Path(temporary) / "POSCAR.manifest.json").read_text())
            self.assertEqual(stored["output"]["sha256"], manifest["output"]["sha256"])
            with self.assertRaisesRegex(BilayerError, "overwrite"):
                write_bilayer(FIXTURE, output, interlayer_distance=3.2, vacuum=20.0)

    def test_scan_manifest_and_cli_failure_status(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "scan"
            manifest = generate_scan(
                FIXTURE,
                root,
                grid=(2, 3),
                interlayer_distance=3.2,
                vacuum=20.0,
            )
            self.assertEqual(manifest["entry_count"], 6)
            self.assertEqual(len(list(root.glob("shift_*/POSCAR"))), 6)
            self.assertEqual(
                main(
                    [
                        "bilayer",
                        "missing",
                        str(root / "bad"),
                        "--interlayer-distance",
                        "3",
                        "--vacuum",
                        "10",
                    ]
                ),
                2,
            )

    def test_scan_force_replaces_directory_without_stale_points(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "scan"
            generate_scan(
                FIXTURE,
                root,
                grid=(2, 2),
                interlayer_distance=3.2,
                vacuum=20.0,
            )
            self.assertEqual(len(list(root.glob("shift_*/POSCAR"))), 4)
            generate_scan(
                FIXTURE,
                root,
                grid=(1, 1),
                interlayer_distance=3.2,
                vacuum=20.0,
                force=True,
            )
            self.assertEqual(len(list(root.glob("shift_*/POSCAR"))), 1)
            generate_scan(
                FIXTURE,
                root,
                grid=(2, 3),
                interlayer_distance=3.2,
                vacuum=20.0,
                force=True,
            )
            self.assertEqual(len(list(root.glob("shift_*/POSCAR"))), 6)

    def test_scan_generation_is_transactional_on_middle_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "scan"
            generate_scan(
                FIXTURE,
                root,
                grid=(1, 1),
                interlayer_distance=3.2,
                vacuum=20.0,
            )
            original_manifest = (root / "scan.manifest.json").read_bytes()
            real_write = bilayer_module.write_bilayer
            calls = 0

            def fail_second(*args, **kwargs):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise BilayerError("injected middle failure")
                return real_write(*args, **kwargs)

            with patch.object(bilayer_module, "write_bilayer", side_effect=fail_second):
                with self.assertRaisesRegex(BilayerError, "injected middle failure"):
                    generate_scan(
                        FIXTURE,
                        root,
                        grid=(2, 2),
                        interlayer_distance=3.2,
                        vacuum=20.0,
                        force=True,
                    )
            self.assertEqual(
                (root / "scan.manifest.json").read_bytes(), original_manifest
            )
            self.assertEqual(len(list(root.glob("shift_*/POSCAR"))), 1)

    def test_scan_resource_limits_fail_before_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "scan"
            with self.assertRaisesRegex(BilayerError, "max_scan_points"):
                generate_scan(
                    FIXTURE,
                    root,
                    grid=(2, 2),
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    max_scan_points=3,
                )
            self.assertFalse(root.exists())
            with self.assertRaisesRegex(BilayerError, "max_total_atoms"):
                generate_scan(
                    FIXTURE,
                    root,
                    grid=(2, 2),
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    max_total_atoms=15,
                )
            self.assertFalse(root.exists())

    def test_cross_layer_pair_ceiling_precedes_allocation_and_distance_is_blocked(
        self,
    ) -> None:
        atoms = read_monolayer(FIXTURE)
        with patch.object(
            bilayer_module, "find_mic", side_effect=AssertionError("allocated")
        ) as mic_mock:
            with self.assertRaisesRegex(BilayerError, "max_cross_layer_pairs"):
                build_bilayer(
                    atoms,
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    max_cross_layer_pairs=3,
                )
            mic_mock.assert_not_called()
        real_find_mic = bilayer_module.find_mic
        with (
            patch.object(bilayer_module, "_DISTANCE_BLOCK_PAIRS", 1),
            patch.object(bilayer_module, "find_mic", wraps=real_find_mic) as mic_mock,
        ):
            result = build_bilayer(
                atoms,
                interlayer_distance=3.2,
                vacuum=20.0,
                max_cross_layer_pairs=4,
            )
        self.assertGreaterEqual(result.minimum_interlayer_distance, 3.2 - 1.0e-10)
        self.assertGreaterEqual(mic_mock.call_count, 12)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "scan"
            with self.assertRaisesRegex(BilayerError, "max_cross_layer_pairs"):
                generate_scan(
                    FIXTURE,
                    root,
                    grid=(1, 1),
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    max_cross_layer_pairs=3,
                )
            self.assertFalse(root.exists())

    def test_output_aliases_and_scan_containment_fail_before_writes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.POSCAR"
            source.write_bytes(FIXTURE.read_bytes())
            original = source.read_bytes()
            with self.assertRaisesRegex(BilayerError, "aliases protected"):
                write_bilayer(
                    source,
                    source,
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    force=True,
                )

            symlink_output = root / "output-link.POSCAR"
            symlink_output.symlink_to(source)
            with self.assertRaisesRegex(BilayerError, "aliases protected"):
                write_bilayer(
                    source,
                    symlink_output,
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    force=True,
                )

            sidecar_source = root / "POSCAR.manifest.json"
            sidecar_source.write_bytes(FIXTURE.read_bytes())
            with self.assertRaisesRegex(BilayerError, "aliases protected"):
                write_bilayer(
                    sidecar_source,
                    root / "POSCAR",
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    force=True,
                )

            scan_root = root / "scan"
            scan_root.mkdir()
            nested_source = scan_root / "input.POSCAR"
            nested_source.write_bytes(FIXTURE.read_bytes())
            nested_original = nested_source.read_bytes()
            with self.assertRaisesRegex(BilayerError, "must not equal or contain"):
                generate_scan(
                    nested_source,
                    scan_root,
                    grid=(1, 1),
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    force=True,
                )
            self.assertEqual(nested_source.read_bytes(), nested_original)

            linked_scan_root = root / "linked-scan"
            linked_scan_root.mkdir()
            linked_source = linked_scan_root / "input.POSCAR"
            linked_source.symlink_to(source)
            with self.assertRaisesRegex(BilayerError, "must not equal or contain"):
                generate_scan(
                    linked_source,
                    linked_scan_root,
                    grid=(1, 1),
                    interlayer_distance=3.2,
                    vacuum=20.0,
                    force=True,
                )
            self.assertTrue(linked_source.is_symlink())
            self.assertEqual(source.read_bytes(), original)

    def test_poscar_manifest_pair_publish_rolls_back_on_failures(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "POSCAR"
            sidecar = root / "POSCAR.manifest.json"
            write_bilayer(
                FIXTURE,
                output,
                interlayer_distance=3.2,
                vacuum=20.0,
            )
            old_pair = (output.read_bytes(), sidecar.read_bytes())
            real_replace = bilayer_module.os.replace
            failed = False

            def fail_manifest_publish(source, destination):
                nonlocal failed
                if (
                    not failed
                    and Path(destination) == sidecar
                    and ".tmp-" in Path(source).name
                ):
                    failed = True
                    raise OSError("injected second publish failure")
                return real_replace(source, destination)

            with patch.object(
                bilayer_module.os, "replace", side_effect=fail_manifest_publish
            ):
                with self.assertRaisesRegex(OSError, "second publish"):
                    write_bilayer(
                        FIXTURE,
                        output,
                        interlayer_distance=3.3,
                        vacuum=20.0,
                        force=True,
                    )
            self.assertEqual((output.read_bytes(), sidecar.read_bytes()), old_pair)

            with patch.object(
                bilayer_module,
                "_write_bytes_fsync",
                side_effect=OSError("injected manifest staging failure"),
            ):
                with self.assertRaisesRegex(OSError, "manifest staging"):
                    write_bilayer(
                        FIXTURE,
                        output,
                        interlayer_distance=3.4,
                        vacuum=20.0,
                        force=True,
                    )
            self.assertEqual((output.read_bytes(), sidecar.read_bytes()), old_pair)
            self.assertEqual(list(root.glob(".*.tmp-*")), [])
            self.assertEqual(list(root.glob(".*.backup-*")), [])


if __name__ == "__main__":
    unittest.main()
