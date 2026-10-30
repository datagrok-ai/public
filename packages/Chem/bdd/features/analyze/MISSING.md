# Chem analyze features — what is still missing

## Chem alerts

Features: `structural-alerts.feature`, `pharmacophores.feature` (30 Oct 2026, local stand, RDKit engine).

### Questions (Chem alerts)

- **No help for the Structural Alerts dialog.** `help/datagrok/solutions/domains/chem/info-panels/structural-alerts.md`
  describes the rule sets and the info pane only; `chem.md` has a "How to use" for the Pharmacophores dialog but
  none for Chem > Analyze > Structural Alerts... (which rule set is ticked by default, the "<rule set> (<column>)" names of
  the columns it adds). The feature claims what `structuralAlertsTopMenu` does; should the help page get a
  "How to use" like the Pharmacophores one?
- **Run twice.** By the code (`runStructuralAlerts`, `getUnusedName`; not measured), a second OK with the same
  rule set adds a second column under a new name rather than refreshing the first. Neither help page says which
  is intended; not claimed.
