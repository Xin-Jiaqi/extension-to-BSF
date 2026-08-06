# extension-to-BSF

`extension-to-BSF` is a Python 3 toolkit for generating auditable bilayer structures and periodic stacking grids from a validated monolayer POSCAR. The current engineering alpha focuses on structure preparation only: it does not run VASP, fit an energy surface, or assert a physical ground state.

## What is safe to use

The supported command-line path uses ASE for POSCAR parsing and writing, including the POSCAR scale factor, Direct/Cartesian coordinates, species ordering, and selective-dynamics constraints. Every generated structure is checked for a valid cell, finite coordinates, explicit vacuum and layer spacing, and a configurable minimum interlayer distance. A JSON manifest records the input digest, parameters, atom counts, minimum interlayer distance, and output digest.

```bash
python -m pip install -e .

# One bilayer at fractional in-plane shift (0.5, 0.0)
extension-bsf bilayer POSCAR build/POSCAR \
  --interlayer-distance 3.2 --vacuum 20 --shift 0.5 0.0 \
  --max-cross-layer-pairs 4000000

# A periodic 9 x 5 grid. Endpoints are not duplicated.
extension-bsf scan POSCAR scan --grid 9 5 \
  --interlayer-distance 3.2 --vacuum 20 \
  --max-scan-points 10000 --max-total-atoms 2000000 \
  --max-cross-layer-pairs 4000000
```

Commands refuse to overwrite outputs unless `--force` is supplied. Output and sidecar paths that alias the input through spelling, symbolic links, or existing hard links are always rejected, and a scan directory may not contain its input. The input bytes and digest are captured before any output write. Each POSCAR and sidecar is staged, fsynced, verified, and published as a rollback-capable pair; scan trees are likewise replaced only after every entry succeeds. Point-count, total-generated-atom, and cross-layer-pair ceilings are checked before output. `max_cross_layer_pairs` defaults to 4,000,000; accepted distance checks run in blocks of at most 65,536 pairs. Python-level failures restore the old pair or remove the incomplete new pair.

## Scope and limitations

- Inputs must be VASP 5+ POSCAR/CONTCAR files with an invertible three-dimensional cell.
- The first two lattice vectors define the periodic plane. The third vector must be parallel to its normal; sheared-vacuum cells fail closed.
- The two layers are copies of the same monolayer and share one in-plane cell. Heterostructural lattice matching belongs in [`heterojunction`](https://github.com/Xin-Jiaqi/heterojunction).
- `interlayer_distance` is the closest separation of the two layer envelopes along the surface normal, not a predicted equilibrium distance.
- Collision screening is geometric. It is not a bonding model and does not replace relaxation or convergence testing.
- Raising the atom, scan-point, or pair ceilings increases memory, storage, and runtime risk.
- The package does not yet consume `materials-structure-core` directly. The interoperability contract is documented in [STRUCTURE_CONTRACT.md](https://github.com/Xin-Jiaqi/extension-to-BSF/blob/main/docs/STRUCTURE_CONTRACT.md) so that a future adapter can be added without copying its implementation.

## Legacy research scripts

The original MATLAB, Python, Bash, and Slurm scripts are preserved in the [legacy directory](https://github.com/Xin-Jiaqi/extension-to-BSF/tree/main/legacy) as provenance. They are not supported entry points. They contain known errors involving POSCAR scale factors, Cartesian/Direct transformations, atom grouping, periodic wrapping, hard-coded paths, and inconsistent 9 x 5 versus 8 x 5 scan grids. See the [legacy notice](https://github.com/Xin-Jiaqi/extension-to-BSF/blob/main/legacy/README.md) before comparing historical results.

## Development and validation

```bash
python -m unittest discover -s tests -v
python -m build
```

The tests include scaled Cartesian POSCAR input, selective dynamics, round-trip writing, collision rejection, periodic-grid uniqueness, overwrite protection, and CLI smoke coverage. Passing tests establish the documented software contract only; independent numerical review is still required before using generated structures in a publication.

## Related tools

ASE supplies the maintained VASP parser/writer used here. Pymatgen offers broader interface construction and Zur--McGill lattice matching, Twister targets commensurate moire superlattices and structural relaxation, and atomate2 orchestrates provenance-rich first-principles workflows. This repository deliberately stays smaller: it prepares same-cell bilayer candidates and manifests that can later feed those workflow systems. Details and citations are in [ECOSYSTEM.md](https://github.com/Xin-Jiaqi/extension-to-BSF/blob/main/docs/ECOSYSTEM.md).

The package is available under the [BSD 3-Clause License](LICENSE). Public release remains **BLOCKED** pending the ownership and disclosure review described in [RELEASE_GATES.md](https://github.com/Xin-Jiaqi/extension-to-BSF/blob/main/RELEASE_GATES.md).
