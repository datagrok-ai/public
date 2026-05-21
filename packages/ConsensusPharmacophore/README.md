# Consensus Pharmacophore

[![Status](https://img.shields.io/badge/status-WIP-yellow)](#)

A Datagrok app that builds a **ligand-based 3D consensus pharmacophore** from a list of co-crystal PDB structures
of a shared target.

## What it does

Given a few PDB IDs of a target with bound ligands (e.g. five EGFR kinase structures), the pipeline:

1. **Enriches** each PDB with RCSB metadata, ligand identity, and quality filters (resolution, R-free, ligand RSCC).
2. **Aligns** all structures onto a global Cα frame (pass 1, Kabsch via MDAnalysis).
3. **Isolates** the binding pocket (atoms within 5 Å of any ligand; DBSCAN available for power users).
4. **Re-aligns** on the pocket-Cα subset (pass 2) for a tighter superposition.
5. **Extracts** pharmacophore features per ligand using the bundled 7-family SMARTS dictionary
   (Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative, Halogen).
6. **Clusters** features per family with k-means and emits a consensus model (centroids + frequency).
7. **Renders** the consensus as a fake-PDB block with B-factor encoding frequency — the docked Mol\* viewer
   colors it automatically.

A mid-pipeline cluster picker surfaces conformational sub-populations (e.g. active vs. inactive kinase
conformations) before pocket-level alignment runs.

## Quick start

- **App:** `Apps → Bio → Consensus Pharmacophore`
- **Demo:** the input panel offers "Use demo input" — loads 5 active-state EGFR kinase PDBs (1XKK / 3W2S / 4WKQ / 5HG5 / 5HG8). All five sit in the DFG-in active conformation, so the conformational cluster picker auto-skips. To exercise the picker, paste a mixed-state set, e.g. add `2HZ0` (DFG-out, 2.10 Å) and `2RGP` (αC-out, 2.00 Å) to the textarea before clicking Build.
- **Custom run:** paste a list of PDB IDs (one per line or comma-separated), pick QC thresholds, hit **Build**.

## Status

Phase 0 — package scaffolded. Implementation in progress; see `blueprint/ligand-based-3d-pharmacophore-datagrok/blueprint.md`.

## License

The TypeScript and Python sources are MIT-licensed under the Datagrok public repository terms. Stage 2/2b
align scripts invoke MDAnalysis (GPL-2.0) at the process boundary inside the Python script handler — no source
linkage.
