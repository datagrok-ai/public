---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p0
realizes: [bio.cp.convert-notation-roundtrip]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/convert.md
migration_date: 2026-05-31
realized_as:
  - convert-spec.ts
source_text_fixes:
  - missing-step-4-renumbered-1-2-3-5-to-1-2-3-4
  - implicit-four-action-by-three-notation-matrix-made-explicit
  - check-a-new-column-expanded-to-per-cell-result-column-assertions
  - get-region-menu-label-aligned-to-atlas-extract-region
  - set-options-clarified-as-accept-defaults-or-prefilled-dialog
candidate_helpers: []
unresolved_ambiguities:
  - polytool-vs-transform-convert-submenu-path
scope_reductions: []
related_bugs:
  - GROK-15176
  - GROK-12164
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T12:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T13:00:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T15:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T22:47:00Z
    spec_runs:
      - spec: convert-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 207
        failure_keys: []
---

# Bio | Transform / Calculate — convert-notation matrix

Covers the four `Bio` transform / calculate top-menu actions that
produce a new column from a Macromolecule column, exercised across
three canonical Macromolecule notations. Verifies that each action
dispatches its dialog and produces a new column with the right shape.

Matrix shape: **4 actions x 3 notations = 12 cells**.

Per-function deeper scenarios (Sequence Space, Activity Cliffs,
Composition) live elsewhere; this scenario is the convert-notation
matrix runner.

## Setup

Datasets (canonical Bio test fixtures, one per source notation):

- `System.AppData/Bio/tests/filter_FASTA.csv` — FASTA notation.
- `System.AppData/Bio/tests/filter_HELM.csv` — HELM notation.
- `System.AppData/Bio/tests/filter_MSA.csv` — separator-with-MSA-flag.

Actions (each invoked from the `Bio` top-menu):

- **Calculate > Extract Region...** (menu label is "Extract
  Region..."; the source scenario spelling "Get region" is
  reconciled here to this label).
- **PolyTool > Convert** (source-scenario spelling). The platform
  registers the convert-notation surface under
  `Bio | Transform | Convert Sequence Notation...`; the
  "PolyTool > Convert" entry in the source scenario is preserved as
  the current PolyTool-namespaced placement.
- **Transform > To Atomic Level...**
- **Transform > Split to Monomers...**

Per-cell expected output column (guidance from the 2026-04-23
reference run in `convert-run.md`):

| Action                        | FASTA result column            | HELM result column           | MSA result column         |
|-------------------------------|--------------------------------|------------------------------|---------------------------|
| Calculate > Extract Region    | `Sequence: (1-39)` (Macromol.) | `HELM: (1-17)` (Macromol.)   | `MSA: (1-17)` (Macromol.) |
| PolyTool > Convert            | `molfile(Sequence)` (Molecule) | `molfile(HELM)` (Molecule)   | `molfile(MSA)` (Molecule) |
| Transform > To Atomic Level   | `molfile(Sequence) (2)` (Mol.) | `molfile(HELM) (2)` (Mol.)   | second molfile column (Mol.) |
| Transform > Split to Monomers | 39 `Monomer` columns (1..39)   | 17 `Monomer` columns (1..17) | 17 `Monomer` columns (1..17) |

These are guidance values — per-cell assertions verify column
appearance + semType, not the exact name string when the column
suffix is environment-dependent.

## Scenarios

### Scenario 1 — Open dataset and verify Macromolecule detection

For each dataset in
`{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The
   Macromolecule detector classifies the sequence column
   synchronously; the table view opens with the sequence column
   tagged `quality=Macromolecule, units=<fasta|helm|separator>`.

### Scenario 2 — Dispatch each transform action and verify a new column appears

For each of the four actions x each of the three datasets (12 cells):

2. On the menu ribbon, open **Bio** and dispatch the action:
   - **Calculate > Extract Region...** -> opens `dialog-Get-Region`
     (`GetRegionEditor`).
   - **PolyTool > Convert** -> opens the convert dialog. Current
     deployment uses a shared widget with the To Atomic Level dialog;
     `name=` attribute observed as `dialog-To-Atomic-Level` per
     `convert-run.md` 2026-04-23.
   - **Transform > To Atomic Level...** -> opens
     `dialog-To-Atomic-Level`.
   - **Transform > Split to Monomers...** -> opens
     `dialog-Split-to-Monomers` (`SplitToMonomersEditor`).

3. Set options:
   - When the dialog opens prefilled with the active table's
     Macromolecule column and any default parameters, leave the
     prefilled values and click **OK**.
   - When the dialog requires a column input that is not auto-bound
     (rare on these fixtures because the detector already tagged the
     sequence column), select the table's sequence column from the
     dropdown, then click **OK**.

4. Verify a new column appears on the dataframe matching the expected
   per-cell semantics in the Setup table:
   - **Extract Region** -> a sub-region Macromolecule column appears
     (name pattern `<source>: (start-end)`, `units` matches source
     notation).
   - **PolyTool > Convert** -> a molfile column appears (semType
     `Molecule`, units `molblock`). Per `convert-run.md`, this entry
     currently reuses the To-Atomic-Level widget — the "new column"
     intent is satisfied by molfile-column appearance, not by widget
     identity.
   - **Transform > To Atomic Level** -> a molfile column appears
     (semType `Molecule`, units `molblock`); when PolyTool > Convert
     ran first in the same dataset, the new column lands as
     `molfile(<source>) (2)`.
   - **Transform > Split to Monomers** -> N per-position `Monomer`
     columns appear (N = alignment width of the source column;
     `monomerCellRenderer` dispatches on the new columns).

   Between actions on the same dataset, wait for the previously
   opened dialog to detach before triggering the next so that two
   actions sharing a `name=` attribute (PolyTool > Convert and
   Transform > To Atomic Level both expose
   `[name="dialog-To-Atomic-Level"]`) do not collide.

## Notes

- This scenario checks that the convert action produces a new column,
  but not that the cell renderer switches correctly afterward —
  that's the dedicated regression guard in `bio-renderer-dispatch.md`
  (GROK-12164). Similarly, the To Atomic Level molfile-validity
  contract for GROK-15176 is exercised more thoroughly in
  `bio-transform-atomic-level.md`.
- Open question: the source scenario lists Convert under "PolyTool >
  Convert", but the Transform menu also registers a "Convert Sequence
  Notation..." entry under a different submenu path. Either the menu
  hierarchy changed over time, or this scenario predates the current
  Transform-submenu placement — the PolyTool path is confirmed
  working on the live dev server as of 2026-04-23, so this scenario
  keeps that spelling until the menu path is confirmed.

---
{
  "order": 5
}
