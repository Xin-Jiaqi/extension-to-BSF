# extension-to-BSF

Research prototype for preparing two-dimensional structures, constructing transformed bilayers, and organizing VASP calculations of a stacking-energy surface.

> [!WARNING]
> This repository is not production-ready. The current checks cover source syntax only; they do not validate structures, symmetry operations, VASP inputs, energies, or scientific conclusions. Known coordinate and POSCAR correctness issues are listed below. Do not use generated structures for production calculations until those issues have been fixed and checked against trusted reference cases.

## Workflow

The scripts currently form a working-directory-based pipeline rather than an installed package:

1. `getCIFs.py` reads `c2db-2022-11-30.db`, groups records by layer-group number, and exports CIF and VASP structure files.
2. `lattice info--oblique.py` and `lattice info--rectangular.py` extract in-plane lattice metrics and write screening results to `vasp_results.xlsx`.
3. `convert_poscar_cartesian_to_direct.m`, `create_transformed_poscars.m`, and `create_bilayer_vasp_files.m` convert or transform coordinates and construct bilayer POSCAR files.
4. `scan_energy/swep.m` creates a lateral-displacement grid, `scan_energy/pos2opt.py` prepares calculation directories, the Slurm scripts submit VASP jobs, and `scan_energy/extract_energy_data.sh` extracts the last reported `TOTEN` values.

Most paths and parameters are currently hard-coded. Read each script before running it, use a disposable working directory, and retain the original structures separately.

## Requirements

Python dependencies are listed in `requirements.txt`:

- ASE
- NumPy
- openpyxl
- pandas

The Python scripts also use the standard-library `json`, `os`, `shutil`, `sqlite3`, and related modules. The shell utilities use Bash, `grep`, `awk`, and `tail`.

Some stages additionally require software that is not distributed by this repository:

- a local C2DB SQLite database;
- MATLAB for the `.m` files (no minimum version has been established);
- a licensed VASP installation, MPI, and a Slurm cluster for the calculation scripts.

No MATLAB, VASP, database, scheduler, or MPI version is implied by the static CI check.

## Current support boundary

- `getCIFs.py` exports layer groups LG8-LG48.
- `lattice info--oblique.py` reads LG1-LG7, while `lattice info--rectangular.py` reads generated `lgnum_*` directories and reports LG8-LG48 statistics.
- The repository does not currently provide a closed end-to-end path for LG1-LG7, and it does not cover LG49-LG80.
- The MATLAB POSCAR readers assume a simplified VASP 5 layout and do not support every valid POSCAR variant, including the optional `Selective dynamics` line.
- The stacking scan has a 9 x 5 generator but 8 x 5 directory preparation and extraction. The intended treatment of the periodic endpoint is not yet encoded consistently.

## Known correctness blockers

The following issues are known and intentionally remain unfixed in this publication-baseline branch:

1. `convert_poscar_cartesian_to_direct.m` changes the scale factor to `1.0` without scaling the lattice and divides the fractional coordinates by the original scale factor. `create_bilayer_vasp_files.m` uses the same extra division for Cartesian input.
2. `create_transformed_poscars.m` detects Cartesian coordinates but then applies transformations, wrapping, and z shifts as if they were fractional coordinates, while retaining the Cartesian label in its output.
3. Both lattice-screening scripts ignore the POSCAR scale factor when reporting lengths. Because the screening threshold is expressed as an absolute length, this can change classification results.
4. The custom VASP writer in `getCIFs.py` counts each element but preserves the incoming atom order. Structures whose atoms are not already grouped by element can therefore be assigned incorrectly. Files are named only by formula, so polymorphs with the same formula can overwrite one another.

Additional known risks include incomplete POSCAR validation, one-step rather than modulo-based coordinate wrapping, repeated element groups in generated bilayers, hard-coded layer identification by `z > 0.5`, and energy extraction without a convergence check.

## Data and cluster configuration

C2DB databases, VASP outputs, POTCAR files, generated structures, Excel results, and private cluster configuration are intentionally excluded by `.gitignore`. They must be supplied locally and used according to their respective licenses and institutional policies.

The tracked Slurm examples contain site-specific account, partition, task-count, and executable settings. Treat them as examples only. Before use, copy the relevant values into a private local configuration, remove personal or institutional identifiers, and verify the submission behavior on the target cluster. Never commit credentials, access tokens, private keys, licensed POTCAR data, or unpublished calculation results.

## Validation status

GitHub Actions compiles the Python sources and runs `bash -n` on the shell scripts. This detects syntax errors only. There are currently no reference structures, numerical regression tests, MATLAB tests, VASP integration tests, or scientific acceptance criteria.

No open-source license has been selected in this baseline; publication of the source does not by itself grant reuse rights.
