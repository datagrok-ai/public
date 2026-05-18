# Process CRO Results

Triage a freshly returned CRO assay: standardize compounds, separate known from novel hits, flag liabilities, score the promising ones, and hand off a shared project for review.

Work through the steps below in order. After each one, post a one-line status update to chat.

1. Load the attached CSV (columns: Compound_ID, SMILES, Assay_pIC50).
2. Standardize structures — canonicalize SMILES and strip salts.
3. Deduplicate against the corporate compound catalog by canonical SMILES. Tag each row as `Known` or `Novel`.
4. Flag structural alerts on all rows: PAINS, BMS, SureChEMBL. Add a `Flags` column.
5. For Novel + unflagged rows, score toxicity risks and drug-likeness (Lipinski).
6. Create a project named "CRO Triage {today}" with the current table view.
7. Share the project with the `toxicology-review` group at View access.
8. Post the project URL to chat.