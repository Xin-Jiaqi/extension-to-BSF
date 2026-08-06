# Ecosystem position and applications

## Maintained building blocks

- [ASE VASP I/O](https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#vasp) reads POSCAR/CONTCAR cells, positions, atom types, and constraints and writes either Direct or Cartesian coordinates. This project uses it instead of maintaining another POSCAR parser.
- [pymatgen interface analysis](https://pymatgen.org/pymatgen.analysis.interfaces.html) provides `ZSLGenerator` and `CoherentInterfaceBuilder` for broader lattice-matched film/substrate interfaces.
- [Twister](https://arxiv.org/abs/2102.07884) constructs commensurate moire superlattices using a coincidence-lattice method and supports relaxation through classical force fields.
- [atomate2](https://materialsproject.github.io/atomate2/user/index.html) composes first-principles jobs and stores structured task records at high-throughput scale.

## Deliberate niche

`extension-to-BSF` is a small pre-processing component for same-cell 2D bilayers. Its value is the auditable translation from one monolayer and an explicit stacking displacement to a collision-screened POSCAR and manifest. It can support stacking-energy data generation, registry comparison, ferroelectric-state candidate preparation, and regression datasets, but it does not decide which displacement is physically stable.

## Integration path

1. A structure source or database supplies a monolayer with provenance.
2. This package enumerates a half-open periodic shift grid and validates geometry.
3. A workflow engine such as atomate2 may consume each manifest/POSCAR pair.
4. Calculation convergence and relaxation status are added by the workflow layer.
5. Energy-surface analysis is performed only on complete, converged records.

Scientific validation still needs reference bilayers with independently checked spacing conventions, relaxed structures, energies, and periodic symmetry equivalences.
