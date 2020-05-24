<!-- TITLE: Molecular Fingerprints -->
<!-- SUBTITLE: -->

# Molecular Fingerprints

Fingerprints are a very abstract representation of certain structural features 
of a molecule. 

Similarity measures, calculations that quantify the similarity of two molecules, and screening, 
a way of rapidly eliminating molecules as candidates in a substructure search, are both processes 
that use fingerprints.

## Available fingerprinters

| Fingerprinter | Description |
|--------|-------------|
| RDKFingerprint | RDKit topological fingerprint for a molecule |
| MACCSKeys | SMARTS-based implementation of the 166 public MACCS keys |
| AtomPair | Atom-Pair fingerprints of molecule, hashed |
| TopologicalTorsion | Topological torsion fingerprints of molecule, hashed |
| Morgan/Circular | Family of fingerprints, better known as circular fingerprints, is built by applying the Morgan algorithm to a set of user-supplied atom invariants |


See also:

  * [Fingerprints - Screening and Similarity](http://www.daylight.com/dayhtml/doc/theory/theory.finger.html)
  * [RDKFingerprint](https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity)
  * [MACCSKeys](https://pdfs.semanticscholar.org/ad40/b25e38314f39a82f193dc4806e6a1c2c6b69.pdf)
  * [AtomPair](https://pubs.acs.org/doi/abs/10.1021/ci00046a002)
  * [TopologicalTorsion](https://pubs.acs.org/doi/abs/10.1021/ci00054a008)
  * [Morgan/Circular](https://pubs.acs.org/doi/10.1021/ci100050t)
