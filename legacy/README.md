# Legacy research scripts

These files are retained to make historical work reproducible and reviewable. They are not supported commands and are excluded from the package.

Known blockers include incorrect handling of non-unit POSCAR scale factors, applying fractional transformations to Cartesian coordinates, single-step wrapping instead of modulo wrapping, species counts that can disagree with atom order, formula-based filename collisions, hard-coded layer thresholds and cluster paths, copying directories without input validation, and inconsistent scan endpoints. The shell energy extractor also does not establish electronic or ionic convergence.

Do not patch these scripts into a production workflow. Compare any historically generated output against the Python 3 implementation and an independently reviewed reference structure.

`requirements.txt` records the unpinned Python packages used historically. MATLAB, VASP, Slurm, MPI, and the C2DB database are external and are not supplied here. Site-specific account, partition, executable, and user details have been removed; the retained Slurm file is only a placeholder template.
