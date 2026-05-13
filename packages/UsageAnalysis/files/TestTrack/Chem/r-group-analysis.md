---
feature: chem
sub_features_covered: [chem.analyze.r-groups, chem.analyze.r-groups.top-menu, chem.analyze.r-groups.decomposition, chem.sketcher]
target_layer: playwright
coverage_type: edge
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/r-group-analysis.md
migration_date: 2026-05-11
migration_report: r-group-analysis-migration-report.md
related_bugs: [GROK-16329]
---

# R-Groups Analysis

Bug-focused R-Group Analysis scenario locking in two invariants:

- **Block A — GROK-16329 empty-result repro.** On `smiles.csv`, the R-Groups dialog
  with MCS + **Visual analysis** must produce a graceful `No R-Groups were found`
  balloon — NOT a null-reference crash. Block A is the literal reproduction path for
  bug-library `GROK-16329` (fixed in 1.20.0, status `fixed`, `test_coverage: needed`).
- **Block B — Replace Latest matrix on `sar_small.csv`.** Three sequential R-Group
  runs with **Replace latest** toggled OFF then ON, plus an OK-without-MCS run that
  must surface the `No core was provided` balloon (negative-path invariant guarding
  the empty-core branch in `rGroupDecomp`).

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario
(`depends_on: []`), `classification: medium`, `pyramid_layer: bug-focused`,
`target_layer: playwright`, strategy `simple`. UI coverage owned
(`ui_coverage_delegated_to: null`) over the R-Groups editor surface: MCS button,
Visual analysis checkbox, Replace latest checkbox, No-core balloon, No-rgroups
balloon, trellis plot rendering. Overlaps atlas critical path
`chem.cp.r-groups-analysis` (p1; "MCS-strategy empty case is the sharp edge").

## Setup

1. **Provision linked datasets.** Both files are bundled in the platform's System
   file shares (no external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — Block A dataset; structurally diverse
     SMILES set on which MCS cannot decompose meaningful R-groups (the empty-result
     trigger for GROK-16329).
   - `System:DemoFiles/chem/sar_small.csv` — Block B dataset; small SAR matrix
     suitable for MCS decomposition with visible R-group columns.
2. **Login** — standard test user. No special privileges required (no project
   save, no share).

## Scenarios

### Block A — `smiles.csv` empty-result balloon (GROK-16329 repro)

Locks the invariant: MCS on diverse SMILES → graceful empty-result message, never
null crash.

1. Open `System:DemoFiles/chem/smiles.csv`. Wait for the table view to render with
   the molecule column populated.
2. Run **Chem > Analyze > R-Groups Analysis...** from the top menu. The
   **R-Groups Analysis** dialog opens with the sketcher, the MCS button, and the
   form (column input, column prefix, **Visual analysis** checkbox).
3. Click the **MCS** button. The platform computes the Most Common Substructure
   and populates the sketcher with the MCS molfile.
4. Confirm the **Visual analysis** checkbox is checked (it defaults to `true`).
5. Click **OK** to run decomposition.
6. **Expected:** a balloon `No R-Groups were found` is shown. The dialog closes
   without throwing a null-reference error in the console. No trellis plot is
   created. The grid acquires no R-group columns.

### Block B — `sar_small.csv` Replace Latest matrix

Walks the three-pass Replace Latest matrix and the empty-core negative-path
balloon.

1. Open `System:DemoFiles/chem/sar_small.csv`. Wait for the table view to render.
2. Run **Chem > Analyze > R-Groups Analysis...**. The dialog opens.
3. Click the **MCS** button. Sketcher populates with the MCS molfile.
4. Click **OK**. **Expected:** a trellis plot with R-group results is added to the
   view; R-group columns are appended to the grid.
5. Run **R-Groups Analysis** once more (top menu re-invocation). The dialog opens
   again. The **Replace latest** checkbox appears because a prior analysis exists
   for this dataframe.
6. Click **MCS** to re-populate the sketcher with the MCS molfile.
7. **Uncheck** the **Replace latest** checkbox, then click **OK**. **Expected:** a
   second trellis plot is displayed in addition to the first. The grid now
   contains two sets of R-Group columns (the first run's columns are preserved;
   the second run's columns are appended with a suffix).
8. Run **R-Groups Analysis** a third time.
9. Click **MCS** to re-populate the sketcher.
10. **Check** the **Replace latest** checkbox, then click **OK**. **Expected:** the
    latest results (R-Group columns and the second trellis plot) are removed; the
    third run's columns and trellis plot replace them. The first run's results
    remain intact.
11. Run **R-Groups Analysis** a fourth time. In the opened dialog, do NOT click
    **MCS** (leave the sketcher empty), then click **OK**. **Expected:** a balloon
    `No core was provided` is shown. No new R-Group columns are added. No trellis
    plot is created.

## Notes

- "Visual analysis" checkbox label verified against
  `public/packages/Chem/src/analysis/r-group-analysis.ts:97`
  (`ui.input.bool('Visual analysis', {value: true})`) — current canonical label.
- "Replace latest" checkbox label verified against
  `r-group-analysis.ts:99` — current canonical label.
- "No core was provided" balloon text verified against
  `r-group-analysis.ts:241` (`grok.shell.error('No core was provided')`).
- "No R-Groups were found" balloon text verified against
  `r-group-analysis.ts:192` (`grok.shell.error('No R-Groups were found')`) —
  note the slight casing/hyphenation difference from the original scenario's
  "No R Groups were found"; the canonical product wording is `No R-Groups were
  found`.
- The **Replace latest** checkbox is rendered conditionally — it only appears
  when `latestAnalysisCols[col.dataFrame.name]?.length` is truthy
  (`r-group-analysis.ts:162`). Block B Step 5 implicitly relies on this: the
  re-invocation after Step 4 is what surfaces the checkbox.
- Existing sibling spec at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/r-group-analysis-spec.ts`
  (per `existing-test-index.yaml`) — Automator should regenerate or extend the
  spec from this migrated scenario.
- No helpers from `helpers-registry.yaml` apply to this scenario (R-Groups dialog
  driving has no registered helper).
