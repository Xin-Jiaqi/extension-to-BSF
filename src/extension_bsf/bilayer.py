"""Fail-closed construction of same-cell bilayers and stacking grids."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from io import StringIO
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Iterable
from uuid import uuid4

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.io import read, write

MANIFEST_SCHEMA = "extension-bsf/bilayer-manifest/v1"
_GEOMETRY_TOL = 1.0e-8
_DISTANCE_BLOCK_PAIRS = 65_536


class BilayerError(ValueError):
    """Raised when an input or generated bilayer violates the public contract."""


@dataclass(frozen=True)
class BilayerBuild:
    atoms: Atoms
    minimum_interlayer_distance: float
    internal_interlayer_distance: float
    outer_periodic_distance: float
    bottom_indices: tuple[int, ...]
    top_indices: tuple[int, ...]


def _positive_finite(value: float, name: str) -> float:
    parsed = float(value)
    if not np.isfinite(parsed) or parsed <= 0:
        raise BilayerError(f"{name} must be finite and greater than zero")
    return parsed


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest_bytes(value: bytes) -> str:
    return sha256(value).hexdigest()


def _resolved_path(path: Path, name: str) -> Path:
    try:
        return path.resolve(strict=False)
    except (OSError, RuntimeError) as exc:
        raise BilayerError(f"could not resolve {name} path {path}: {exc}") from exc


def _paths_alias(left: Path, right: Path) -> bool:
    if _resolved_path(left, "path") == _resolved_path(right, "path"):
        return True
    try:
        return left.exists() and right.exists() and left.samefile(right)
    except OSError as exc:
        raise BilayerError(
            f"could not compare paths {left} and {right}: {exc}"
        ) from exc


def _reject_output_aliases(
    outputs: tuple[tuple[Path, str], ...], protected: tuple[tuple[Path, str], ...]
) -> None:
    for output, output_name in outputs:
        for source, source_name in protected:
            if _paths_alias(output, source):
                raise BilayerError(
                    f"{output_name} path aliases protected {source_name} path: {output}"
                )
    for index, (left, left_name) in enumerate(outputs):
        for right, right_name in outputs[index + 1 :]:
            if _paths_alias(left, right):
                raise BilayerError(
                    f"{left_name} path aliases {right_name} path: {left}"
                )


def _temporary_path(parent: Path, name: str) -> Path:
    descriptor, value = tempfile.mkstemp(prefix=f".{name}.tmp-", dir=parent)
    os.close(descriptor)
    return Path(value)


def _fsync_file(path: Path) -> None:
    with path.open("rb") as stream:
        os.fsync(stream.fileno())


def _write_bytes_fsync(path: Path, value: bytes) -> None:
    with path.open("wb") as stream:
        stream.write(value)
        stream.flush()
        os.fsync(stream.fileno())


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish_pair(
    temporary_output: Path,
    destination: Path,
    temporary_manifest: Path,
    manifest_path: Path,
    *,
    force: bool,
) -> None:
    targets = (destination, manifest_path)
    temporaries = (temporary_output, temporary_manifest)
    published: list[Path] = []
    backups: dict[Path, Path] = {}
    original_exists = {path: os.path.lexists(path) for path in targets}
    try:
        if not force:
            try:
                for temporary, target in zip(temporaries, targets):
                    os.link(temporary, target)
                    published.append(target)
            except FileExistsError as exc:
                raise BilayerError(
                    f"refusing to overwrite existing output pair: {destination}"
                ) from exc
            _fsync_directory(destination.parent)
            return

        for target in targets:
            if original_exists[target]:
                backup = target.parent / f".{target.name}.backup-{uuid4().hex}"
                os.replace(target, backup)
                backups[target] = backup
        os.replace(temporary_output, destination)
        os.replace(temporary_manifest, manifest_path)
        _fsync_directory(destination.parent)
    except Exception as publish_error:
        try:
            if not force:
                for target in published:
                    if os.path.lexists(target):
                        target.unlink()
            else:
                for target in targets:
                    backup = backups.get(target)
                    if backup is not None and os.path.lexists(backup):
                        if os.path.lexists(target):
                            target.unlink()
                        os.replace(backup, target)
                    elif not original_exists[target] and os.path.lexists(target):
                        target.unlink()
            _fsync_directory(destination.parent)
        except Exception as rollback_error:
            retained = [str(path) for path in backups.values() if os.path.lexists(path)]
            raise BilayerError(
                "output-pair rollback failed; retained backups: "
                + (", ".join(retained) if retained else "none")
            ) from rollback_error
        raise publish_error
    else:
        for backup in backups.values():
            if os.path.lexists(backup):
                backup.unlink()
        _fsync_directory(destination.parent)
    finally:
        for temporary in temporaries:
            if os.path.lexists(temporary):
                temporary.unlink()


def _stage_and_publish_poscar_pair(
    destination: Path,
    manifest_path: Path,
    atoms: Atoms,
    manifest: dict[str, object],
    *,
    force: bool,
) -> None:
    temporary_output = _temporary_path(destination.parent, destination.name)
    temporary_manifest: Path | None = None
    try:
        temporary_manifest = _temporary_path(destination.parent, manifest_path.name)
        write(
            temporary_output,
            atoms,
            format="vasp",
            direct=True,
            sort=False,
            vasp5=True,
        )
        _fsync_file(temporary_output)
        output_digest = _sha256(temporary_output)
        manifest["output"] = {"path": str(destination), "sha256": output_digest}
        manifest_bytes = (
            json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        ).encode()
        _write_bytes_fsync(temporary_manifest, manifest_bytes)
        parsed_manifest = json.loads(temporary_manifest.read_bytes().decode("utf-8"))
        if parsed_manifest != json.loads(manifest_bytes.decode("utf-8")) or (
            output_digest != _sha256(temporary_output)
        ):
            raise BilayerError("temporary POSCAR/manifest pair failed verification")
        _publish_pair(
            temporary_output,
            destination,
            temporary_manifest,
            manifest_path,
            force=force,
        )
    finally:
        for temporary in (temporary_output, temporary_manifest):
            if temporary is not None and os.path.lexists(temporary):
                temporary.unlink()


def _stable_species_order(atoms: Atoms) -> tuple[Atoms, list[int]]:
    symbols = np.asarray(atoms.get_chemical_symbols(), dtype=str)
    order = np.argsort(symbols, kind="stable")
    ordered = atoms[order]
    return ordered, [int(value) for value in order]


def _positive_integer(value: int, name: str) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, np.integer))
        or value <= 0
    ):
        raise BilayerError(f"{name} must be a positive integer")
    return int(value)


def _surface_frame(
    atoms: Atoms,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cell = np.asarray(atoms.cell.array, dtype=float)
    if cell.shape != (3, 3) or not np.isfinite(cell).all():
        raise BilayerError("cell must be a finite 3 x 3 matrix")
    a, b, c = cell
    cross = np.cross(a, b)
    area = float(np.linalg.norm(cross))
    if area <= _GEOMETRY_TOL or abs(float(np.linalg.det(cell))) <= _GEOMETRY_TOL:
        raise BilayerError("cell must be invertible and have a non-zero in-plane area")
    normal = cross / area
    if float(np.dot(c, normal)) < 0:
        normal = -normal
    c_parallel = c - np.dot(c, normal) * normal
    if np.linalg.norm(c_parallel) > 1.0e-7 * max(np.linalg.norm(c), 1.0):
        raise BilayerError(
            "the third lattice vector must be parallel to the surface normal; "
            "sheared-vacuum cells are not supported"
        )
    positions = np.asarray(atoms.positions, dtype=float)
    if positions.shape != (len(atoms), 3) or not np.isfinite(positions).all():
        raise BilayerError("atomic positions must be finite N x 3 coordinates")
    return a, b, c, normal


def _read_monolayer_and_digest(path: str | Path) -> tuple[Atoms, str]:
    source = Path(path)
    if not source.is_file():
        raise BilayerError(f"input POSCAR does not exist: {source}")
    try:
        content = source.read_bytes()
        atoms = read(StringIO(content.decode("utf-8")), format="vasp")
    except Exception as exc:  # ASE exposes format-specific exception types.
        raise BilayerError(f"could not parse POSCAR {source}: {exc}") from exc
    if not isinstance(atoms, Atoms) or len(atoms) == 0:
        raise BilayerError("input must contain at least one atom")
    _surface_frame(atoms)
    atoms.pbc = (True, True, True)
    return atoms, _digest_bytes(content)


def read_monolayer(path: str | Path) -> Atoms:
    """Read a POSCAR with ASE and enforce the supported 2D-cell boundary."""
    return _read_monolayer_and_digest(path)[0]


def _unwrapped_normal_heights(atoms: Atoms, normal: np.ndarray) -> np.ndarray:
    """Unwrap a finite slab across z using its unique largest vacuum gap."""
    if len(atoms) == 1:
        return np.zeros(1, dtype=float)
    fractional_z = np.mod(atoms.get_scaled_positions(wrap=False)[:, 2], 1.0)
    order = np.argsort(fractional_z, kind="stable")
    sorted_z = fractional_z[order]
    gaps = np.diff(np.concatenate((sorted_z, [sorted_z[0] + 1.0])))
    largest = float(np.max(gaps))
    winners = np.flatnonzero(np.isclose(gaps, largest, rtol=0.0, atol=1.0e-10))
    if len(winners) != 1 or largest <= _GEOMETRY_TOL:
        raise BilayerError(
            "cannot uniquely unwrap the layer along z: the largest vacuum gap "
            "is absent or tied"
        )
    start = sorted_z[(int(winners[0]) + 1) % len(sorted_z)]
    unwrapped_fractional = np.mod(fractional_z - start, 1.0)
    normal_cell_length = float(np.dot(atoms.cell.array[2], normal))
    if normal_cell_length <= _GEOMETRY_TOL:
        raise BilayerError("third lattice vector has no positive normal length")
    return unwrapped_fractional * normal_cell_length


def _interlayer_minimum(
    atoms: Atoms,
    bottom: Iterable[int],
    top: Iterable[int],
    *,
    mode: str,
    max_cross_layer_pairs: int,
) -> float:
    bottom_positions = atoms.positions[np.fromiter(bottom, dtype=int)]
    top_positions = atoms.positions[np.fromiter(top, dtype=int)]
    max_cross_layer_pairs = _positive_integer(
        max_cross_layer_pairs, "max_cross_layer_pairs"
    )
    pair_count = len(bottom_positions) * len(top_positions)
    if pair_count > max_cross_layer_pairs:
        raise BilayerError(
            f"cross-layer distance check requires {pair_count} pairs, above "
            f"max_cross_layer_pairs={max_cross_layer_pairs}"
        )
    if pair_count == 0:
        raise BilayerError("could not compute finite interlayer distances")
    if mode == "internal":
        pbc = (True, True, False)
    elif mode == "full-periodic":
        pbc = (True, True, True)
    elif mode == "outer":
        pbc = (True, True, False)
    else:
        raise BilayerError(f"unknown interlayer-distance mode: {mode}")
    top_block_size = min(len(top_positions), _DISTANCE_BLOCK_PAIRS)
    bottom_block_size = max(1, _DISTANCE_BLOCK_PAIRS // top_block_size)
    minimum = np.inf
    for bottom_start in range(0, len(bottom_positions), bottom_block_size):
        bottom_block = bottom_positions[bottom_start : bottom_start + bottom_block_size]
        for top_start in range(0, len(top_positions), top_block_size):
            top_block = top_positions[top_start : top_start + top_block_size]
            vectors = top_block[None, :, :] - bottom_block[:, None, :]
            if mode == "outer":
                vectors = vectors - atoms.cell.array[2]
            _, distances = find_mic(vectors.reshape((-1, 3)), atoms.cell.array, pbc=pbc)
            if distances.size == 0 or not np.isfinite(distances).all():
                raise BilayerError("could not compute finite interlayer distances")
            minimum = min(minimum, float(np.min(distances)))
    return float(minimum)


def build_bilayer(
    monolayer: Atoms,
    *,
    interlayer_distance: float,
    vacuum: float,
    shift: tuple[float, float] = (0.0, 0.0),
    minimum_distance: float = 1.0,
    max_cross_layer_pairs: int = 4_000_000,
) -> BilayerBuild:
    """Create a two-copy bilayer with explicit envelope spacing and vacuum.

    ``shift`` is expressed in fractional coordinates of the first two cell
    vectors. Its components are reduced modulo one so periodic endpoints map
    to one canonical representation.
    """
    interlayer_distance = _positive_finite(interlayer_distance, "interlayer_distance")
    vacuum = _positive_finite(vacuum, "vacuum")
    minimum_distance = _positive_finite(minimum_distance, "minimum_distance")
    max_cross_layer_pairs = _positive_integer(
        max_cross_layer_pairs, "max_cross_layer_pairs"
    )
    if len(shift) != 2 or not np.isfinite(np.asarray(shift, dtype=float)).all():
        raise BilayerError("shift must contain two finite fractional coordinates")

    source = monolayer.copy()
    if len(source) == 0:
        raise BilayerError("monolayer must contain at least one atom")
    a, b, _, normal = _surface_frame(source)
    scaled = source.get_scaled_positions(wrap=False)
    in_plane = np.outer(np.mod(scaled[:, 0], 1.0), a) + np.outer(
        np.mod(scaled[:, 1], 1.0), b
    )
    heights = _unwrapped_normal_heights(source, normal)
    thickness = float(np.ptp(heights))
    normalized_heights = heights - float(np.min(heights))
    bottom_floor = vacuum / 2.0
    top_floor = bottom_floor + thickness + interlayer_distance
    canonical_shift = np.mod(np.asarray(shift, dtype=float), 1.0)
    shift_cart = canonical_shift[0] * a + canonical_shift[1] * b

    bottom_positions = in_plane + np.outer(bottom_floor + normalized_heights, normal)
    top_positions = (
        in_plane + shift_cart + np.outer(top_floor + normalized_heights, normal)
    )
    height = vacuum + 2.0 * thickness + interlayer_distance
    symbols = source.get_chemical_symbols() * 2
    result = Atoms(
        symbols=symbols,
        positions=np.vstack((bottom_positions, top_positions)),
        cell=np.vstack((a, b, normal * height)),
        pbc=(True, True, True),
    )
    count = len(source)
    pair_count = count * count
    if pair_count > max_cross_layer_pairs:
        raise BilayerError(
            f"cross-layer distance check requires {pair_count} pairs, above "
            f"max_cross_layer_pairs={max_cross_layer_pairs}"
        )
    repeated_constraints = []
    for constraint in source.constraints:
        try:
            repeated_constraints.append(constraint.copy().repeat((2, 1, 1), count))
        except Exception as exc:
            raise BilayerError(
                f"cannot safely duplicate constraint {type(constraint).__name__}"
            ) from exc
    if repeated_constraints:
        result.set_constraint(repeated_constraints)
    result.new_array(
        "layer_id",
        np.concatenate((np.zeros(count, dtype=int), np.ones(count, dtype=int))),
    )
    bottom_indices = tuple(range(count))
    top_indices = tuple(range(count, 2 * count))
    internal_distance = _interlayer_minimum(
        result,
        bottom_indices,
        top_indices,
        mode="internal",
        max_cross_layer_pairs=max_cross_layer_pairs,
    )
    outer_distance = _interlayer_minimum(
        result,
        bottom_indices,
        top_indices,
        mode="outer",
        max_cross_layer_pairs=max_cross_layer_pairs,
    )
    distance = _interlayer_minimum(
        result,
        bottom_indices,
        top_indices,
        mode="full-periodic",
        max_cross_layer_pairs=max_cross_layer_pairs,
    )
    if distance + 1.0e-10 < minimum_distance:
        raise BilayerError(
            f"minimum full-periodic cross-layer distance {distance:.6f} A is below the "
            f"configured limit {minimum_distance:.6f} A"
        )
    return BilayerBuild(
        result,
        distance,
        internal_distance,
        outer_distance,
        bottom_indices,
        top_indices,
    )


def write_bilayer(
    input_path: str | Path,
    output_path: str | Path,
    *,
    interlayer_distance: float,
    vacuum: float,
    shift: tuple[float, float] = (0.0, 0.0),
    minimum_distance: float = 1.0,
    max_cross_layer_pairs: int = 4_000_000,
    force: bool = False,
) -> dict[str, object]:
    source = Path(input_path)
    destination = Path(output_path)
    manifest_path = destination.with_name(destination.name + ".manifest.json")
    _reject_output_aliases(
        ((destination, "structure output"), (manifest_path, "bilayer manifest")),
        ((source, "monolayer input"),),
    )
    if not force and (os.path.lexists(destination) or os.path.lexists(manifest_path)):
        raise BilayerError(f"refusing to overwrite existing output pair: {destination}")
    monolayer, input_digest = _read_monolayer_and_digest(source)
    build = build_bilayer(
        monolayer,
        interlayer_distance=interlayer_distance,
        vacuum=vacuum,
        shift=shift,
        minimum_distance=minimum_distance,
        max_cross_layer_pairs=max_cross_layer_pairs,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    output_atoms, source_indices = _stable_species_order(build.atoms)
    canonical_shift = [float(v) for v in np.mod(np.asarray(shift, dtype=float), 1.0)]
    manifest: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "input": {"path": str(source), "sha256": input_digest},
        "parameters": {
            "interlayer_distance_angstrom": float(interlayer_distance),
            "vacuum_angstrom": float(vacuum),
            "fractional_shift": canonical_shift,
            "minimum_distance_angstrom": float(minimum_distance),
            "max_cross_layer_pairs": int(max_cross_layer_pairs),
        },
        "validation": {
            "atom_count_per_layer": len(monolayer),
            "total_atom_count": len(build.atoms),
            "minimum_interlayer_distance_angstrom": build.minimum_interlayer_distance,
            "internal_interlayer_distance_angstrom": build.internal_interlayer_distance,
            "outer_periodic_distance_angstrom": build.outer_periodic_distance,
            "coordinate_mode": "Direct",
            "species_sorted_on_write": True,
            "source_index_by_output_index": source_indices,
            "layer_id_by_output_index": [
                int(value) for value in output_atoms.arrays["layer_id"]
            ],
        },
    }
    _stage_and_publish_poscar_pair(
        destination,
        manifest_path,
        output_atoms,
        manifest,
        force=force,
    )
    return manifest


def periodic_shifts(nx: int, ny: int) -> tuple[tuple[float, float], ...]:
    """Return a canonical periodic grid on [0, 1) x [0, 1)."""
    nx = _positive_integer(nx, "grid x dimension")
    ny = _positive_integer(ny, "grid y dimension")
    shifts = tuple((i / nx, j / ny) for i in range(nx) for j in range(ny))
    if len(set(shifts)) != nx * ny:
        raise BilayerError("internal error: periodic scan contains duplicate endpoints")
    return shifts


def generate_scan(
    input_path: str | Path,
    output_directory: str | Path,
    *,
    grid: tuple[int, int],
    interlayer_distance: float,
    vacuum: float,
    minimum_distance: float = 1.0,
    max_scan_points: int = 10_000,
    max_total_atoms: int = 2_000_000,
    max_cross_layer_pairs: int = 4_000_000,
    force: bool = False,
) -> dict[str, object]:
    root = Path(output_directory)
    source = Path(input_path)
    source_lexical = Path(os.path.abspath(source))
    root_lexical = Path(os.path.abspath(root))
    source_resolved = _resolved_path(source, "monolayer input")
    root_resolved = _resolved_path(root, "scan output directory")
    if (
        source_lexical == root_lexical
        or source_lexical.is_relative_to(root_lexical)
        or source_resolved == root_resolved
        or source_resolved.is_relative_to(root_resolved)
    ):
        raise BilayerError(
            "scan output directory must not equal or contain the monolayer input"
        )
    scan_manifest = root / "scan.manifest.json"
    if root.exists() and not root.is_dir():
        raise BilayerError(f"scan output path is not a directory: {root}")
    if not force and (
        scan_manifest.exists() or (root.exists() and any(root.iterdir()))
    ):
        raise BilayerError(f"refusing to write into non-empty scan directory: {root}")
    try:
        nx, ny = grid
    except (TypeError, ValueError) as exc:
        raise BilayerError("grid must contain exactly two positive integers") from exc
    nx = _positive_integer(nx, "grid x dimension")
    ny = _positive_integer(ny, "grid y dimension")
    max_scan_points = _positive_integer(max_scan_points, "max_scan_points")
    max_total_atoms = _positive_integer(max_total_atoms, "max_total_atoms")
    max_cross_layer_pairs = _positive_integer(
        max_cross_layer_pairs, "max_cross_layer_pairs"
    )
    point_count = nx * ny
    if point_count > max_scan_points:
        raise BilayerError(
            f"scan has {point_count} points, above max_scan_points={max_scan_points}"
        )
    monolayer = read_monolayer(source)
    cross_layer_pair_count = len(monolayer) * len(monolayer)
    if cross_layer_pair_count > max_cross_layer_pairs:
        raise BilayerError(
            f"cross-layer distance check requires {cross_layer_pair_count} pairs, above "
            f"max_cross_layer_pairs={max_cross_layer_pairs}"
        )
    generated_atom_count = point_count * 2 * len(monolayer)
    if generated_atom_count > max_total_atoms:
        raise BilayerError(
            f"scan would generate {generated_atom_count} atoms, above "
            f"max_total_atoms={max_total_atoms}"
        )

    root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{root.name}.tmp-", dir=root.parent))
    backup: Path | None = None
    try:
        entries: list[dict[str, object]] = []
        for i, j_shift in enumerate(periodic_shifts(nx, ny)):
            ix, iy = divmod(i, ny)
            relative_poscar = Path(f"shift_{ix:03d}_{iy:03d}") / "POSCAR"
            poscar = staging / relative_poscar
            manifest = write_bilayer(
                input_path,
                poscar,
                interlayer_distance=interlayer_distance,
                vacuum=vacuum,
                shift=j_shift,
                minimum_distance=minimum_distance,
                max_cross_layer_pairs=max_cross_layer_pairs,
                force=False,
            )
            manifest["output"]["path"] = str(relative_poscar)  # type: ignore[index]
            poscar.with_name("POSCAR.manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n"
            )
            entries.append(
                {
                    "index": [ix, iy],
                    "fractional_shift": list(j_shift),
                    "poscar": str(relative_poscar),
                    "sha256": manifest["output"]["sha256"],  # type: ignore[index]
                }
            )
        payload: dict[str, object] = {
            "schema": "extension-bsf/scan-manifest/v1",
            "grid": [nx, ny],
            "periodic_endpoint_policy": "half-open-[0,1)",
            "entry_count": len(entries),
            "generated_atom_count": generated_atom_count,
            "max_cross_layer_pairs": max_cross_layer_pairs,
            "cross_layer_pair_count_per_structure": cross_layer_pair_count,
            "entries": entries,
        }
        (staging / "scan.manifest.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )

        if root.exists():
            backup = root.parent / f".{root.name}.backup-{uuid4().hex}"
            os.replace(root, backup)
        try:
            os.replace(staging, root)
        except Exception:
            if backup is not None and backup.exists() and not root.exists():
                os.replace(backup, root)
            raise
        if backup is not None and backup.exists():
            shutil.rmtree(backup)
        return payload
    finally:
        if staging.exists():
            shutil.rmtree(staging)
