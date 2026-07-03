---
feature: bio
sub_features_covered:
  - bio.calculate.get-region
  - bio.calculate.get-region.top-menu
  - bio.transform.convert-notation
  - bio.transform.convert-notation.top-menu
  - bio.transform.to-atomic-level
  - bio.transform.split-to-monomers
  - bio.detector
target_layer: playwright
coverage_type: regression
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

Source-matrix scenario for the four `Bio` transform / calculate
top-menu actions that produce a new column from a Macromolecule
column, exercised across three canonical Macromolecule notations.
Verifies that each action dispatches its dialog and produces a new
column whose semantics match the atlas registration.

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

- **Calculate > Extract Region...** (atlas
  `bio.calculate.get-region.top-menu`, `package.ts#L483`; menu label
  is "Extract Region..." per atlas — the source scenario spelling
  "Get region" is reconciled here to the atlas label).
- **PolyTool > Convert** (source-scenario spelling). Atlas registers
  the convert-notation surface under
  `Bio | Transform | Convert Sequence Notation...`
  (`bio.transform.convert-notation.top-menu`, `package.ts#L1131`); the
  "PolyTool > Convert" entry in the source scenario is preserved as
  the current PolyTool-namespaced placement (see
  `unresolved_ambiguities: polytool-vs-transform-convert-submenu-path`
  for the menu-path reconciliation TODO).
- **Transform > To Atomic Level...** (atlas
  `bio.transform.to-atomic-level`, `package.ts#L863`).
- **Transform > Split to Monomers...** (atlas
  `bio.transform.split-to-monomers`, `package.ts#L1225`).

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
   Macromolecule detector classifies the sequence column synchronously
   (atlas `bio.detector`, `detectors.js#L1`); the table view opens
   with the sequence column tagged
   `quality=Macromolecule, units=<fasta|helm|separator>`.

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

- **Source-text fixes silently applied** during migration:
  - The source scenario's step numbering ran `1, 2, 3, 5` (step 4
    absent). Renumbered to `1, 2, 3, 4` in this migration. See
    frontmatter `source_text_fixes`.
  - The implicit 4-action x 3-notation matrix (Step 1 lists three
    datasets; Step 2 lists four sub-menu actions) was not explicit in
    the source; surfaced under Setup as a 12-cell matrix.
  - The source "Check a new column" assertion was expanded into the
    per-cell expected-output table in Setup so the verification
    surface is decidable across the 12 cells (which column, which
    semType / units / shape).
  - Source-scenario menu label "Get region" reconciled to the atlas
    label "Extract Region..." (`bio.calculate.get-region.top-menu`,
    atlas L496).
  - Source-scenario step "Set options" clarified into two cases:
    accept defaults on a prefilled dialog, or select the Macromolecule
    column manually when not auto-bound.
- **Sub-features covered:**
  - `bio.calculate.get-region` + `.top-menu` — Extract Region action
    plus its top-menu entry-point.
  - `bio.transform.convert-notation` + `.top-menu` — Convert Sequence
    Notation surface (current PolyTool > Convert dispatch).
  - `bio.transform.to-atomic-level` — To Atomic Level transform.
  - `bio.transform.split-to-monomers` — Split to Monomers transform.
  - `bio.detector` — Macromolecule classification on dataset open
    (gate condition for the transforms to bind to the correct
    column).
- **Bug-context (`related_bugs`)** surfaced for downstream awareness;
  per-bug coverage is delegated to chain-level
  `bug_focused_candidates`:
  - **GROK-15176** (MSA -> To Atomic Level molfile carries illegal
    isotope=1 on heavy-atom oxygens; downstream PubChem
    standardization rejects). Cross-cutting with `msa.md` /
    `pepsea.md`; chain proposes `bio-grok-15176-spec.ts`. This
    scenario exercises To Atomic Level on canonical sources but does
    NOT chain into Chem renderer / PubChem standardization — the
    molfile-validity contract is covered cross-scenario.
  - **GROK-12164** (cell renderer dispatch failure after HELM ->
    SEPARATOR convert with "-" separator — renderer stays on
    `bilnSequenceCellRenderer` instead of dispatching to
    `separatorSequenceCellRenderer`). This scenario verifies the
    convert action produces a new column but does NOT assert
    post-convert renderer dispatch on the same column. Chain marks
    this bug below the cross-scenario trigger threshold
    (`bug_match_attempts_skipped: GROK-12164`); no dedicated spec
    proposed at chain level. Surfaced here for awareness; renderer-
    dispatch regression is the responsibility of
    `bio.cp.convert-notation-roundtrip` follow-up.
- **Unresolved ambiguity** (carried in frontmatter
  `unresolved_ambiguities`): the source scenario lists Convert under
  "PolyTool > Convert" submenu, but atlas registers
  `bio.transform.convert-notation.top-menu` under
  "Bio | Transform | Convert Sequence Notation..."
  (`package.ts#L1131`). Either the menu hierarchy migrated (current
  PolyTool nesting) or the scenario predates the Transform-submenu
  placement. The reference run `convert-run.md` 2026-04-23 ran
  successfully via the PolyTool path on the live dev server.
  Operator confirmation of current menu path needed before Migrator
  rewrites the action description; this scenario preserves the source
  spelling "PolyTool > Convert" so the assertion can be re-evaluated
  against the live server.
- **See atlas** `bio.cp.convert-notation-roundtrip` (p0),
  `bio.cp.to-atomic-level` (p1), `bio.cp.split-to-monomers` (p1) for
  the canonical critical-path framing; this scenario realizes their
  matrix runner.

---
{
  "order": 5
}
