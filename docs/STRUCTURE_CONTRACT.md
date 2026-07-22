# Structure contract and interoperability

The supported API uses ASE `Atoms` as its runtime adapter and follows this explicit boundary:

- lengths are in angstrom;
- lattice vectors are rows of a 3 x 3 matrix;
- fractional coordinates multiply the row-lattice matrix on the left (`fractional @ lattice`);
- the first two lattice vectors are periodic in-plane vectors;
- the third lattice vector is parallel to the surface normal;
- generated POSCAR files use Direct coordinates and group species on write;
- generated JSON contains only finite, versioned, machine-readable values.

Atoms split across the periodic z boundary are unwrapped at the unique largest vacuum gap. If the largest gap is tied, the physical layer image is ambiguous and construction fails closed. Scan directories are generated transactionally and replaced only after all structures and manifests pass validation.

Input bytes and their SHA-256 digest are captured before any output write. Output and sidecar paths may not alias the input through resolved paths, symbolic links, or existing hard links; a scan directory may not equal or contain its input. POSCAR/manifest pairs are staged in the target directory, fsynced, digest-verified, then published with rollback on Python exceptions. Cross-layer pair counts are checked before distance-vector allocation, recorded in manifests, and accepted distances are evaluated in blocks of at most 65,536 pairs. The default `max_cross_layer_pairs` ceiling is 4,000,000.

This is compatible with the public invariants of `materials-structure-core` but does not copy its model or hash implementation. A future adapter should construct a `StructureRecord`, preserve its result identity in the application manifest, and verify the output before calculation submission. Until such an adapter is released and tested, this repository's SHA-256 digests are file-integrity identifiers, not canonical structure identifiers.
