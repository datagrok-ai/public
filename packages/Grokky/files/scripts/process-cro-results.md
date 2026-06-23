# Process CRO Results

Triage a freshly returned CRO assay: standardize compounds, separate known from novel hits, flag liabilities, score the promising ones, and hand off a shared project for review.

Work through the steps below in order. After each one, post a one-line status update to chat.

1. Open the dataset containing molecular structures.
2. Standardize structures — canonicalize SMILES and strip salts.
3. Flag structural alerts on all rows: PAINS, BMS, SureChEMBL. Add a `Flags` column.
4. For unflagged rows, score toxicity risks and drug-likeness (Lipinski).
5. Build a chemical space map (t-SNE/UMAP on Morgan fingerprints), colored by Lipinski violations.
6. Create a project named "CRO Triage {today}".
7. Share the project with the `toxicology-review` group at View access.
8. Share the project in the chat.