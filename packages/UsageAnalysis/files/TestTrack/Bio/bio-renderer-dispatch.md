---
feature: bio
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [GROK-12164]
produced_from: atlas-driven
realized_as:
  - bio-renderer-dispatch-spec.ts
related_bugs:
  - GROK-12164
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T07:15:19Z
    spec_runs:
      - spec: bio-renderer-dispatch-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 191
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T08:00:00Z
    failure_keys: []
    review_round: 1
---

# Bio | Rendering — non-FASTA cell-renderer dispatch matrix

Covers the Macromolecule cell-renderer family beyond the FASTA path
already exercised by `convert.md` and `composition-analysis.md`. Bio
registers four sequence cell renderers keyed off the `units` column
tag the detector sets: `fastaSequenceCellRenderer` (fasta),
`separatorSequenceCellRenderer` (separator),
`bilnSequenceCellRenderer` (biln/HELM), and
`customSequenceCellRenderer` (custom, for pluggable notation
providers). A fifth renderer, `monomerCellRenderer`, paints the
per-position `Monomer` columns produced by the Split-to-Monomers
transform.

This scenario checks that renderer dispatch follows the detected
`units` tag for each of the three remaining notation classes
(separator, HELM, custom), and that `monomerCellRenderer` paints the
per-position columns from a transform. It's the non-FASTA complement
to the FASTA dispatch covered in `convert.md`.

Scenario 2 is a regression guard for GROK-12164: converting a HELM
column to SEPARATOR notation used to leave the old BILN renderer
dispatched on the new column instead of switching to the separator
renderer. The expected outcome here (post-convert dispatch switches
to `separatorSequenceCellRenderer`) is the fixed behavior.

## Setup

Datasets (canonical Bio test fixtures, one per source notation):

- `System.AppData/Bio/tests/filter_HELM.csv` — HELM notation,
  detector sets `units=helm` (bilnSequenceCellRenderer dispatch).
- `System.AppData/Bio/tests/filter_MSA.csv` — separator-with-
  MSA-flag, detector sets `units=separator`
  (separatorSequenceCellRenderer dispatch).

The Macromolecule cell renderer family is registered through
`grok.functions` with cell-type keyed by the `units` column tag
(`public/packages/Bio/src/package.ts#L302..L344`); the detector
sets `units` synchronously on column open
(`detectMacromolecule`, `public/packages/Bio/detectors.js#L1`).

Helper handle: top-menu Bio entries used in Scenario 2 reuse
the `bio.flow.convert-matrix` helper-registry entry; the
detector classification used in every scenario reuses the
standard `bio.flow` dataset-open dataset path.

## Scenarios

### Scenario 1 — Detector + renderer dispatch for HELM and SEPARATOR

Steps:

1. Open `System.AppData/Bio/tests/filter_HELM.csv` from the
   Files browser. The Macromolecule detector classifies the
   sequence column synchronously and sets column tags
   `quality=Macromolecule, units=helm`.
2. Confirm the sequence column renders via the HELM dispatch
   path. The `bilnSequenceCellRenderer` registration is keyed on
   `units=helm` (HELM strings are BILN-shape under the bio
   library notation hierarchy); cells paint monomer letters
   with the bio-library monomer-color palette.
3. Open `System.AppData/Bio/tests/filter_MSA.csv` from the
   Files browser in a second table view. The detector
   classifies the column and sets
   `units=separator, .separator='-', aligned=SEQ.MSA`.
4. Confirm the separator dispatch path. The
   `separatorSequenceCellRenderer` registration is keyed on
   `units=separator`; cells paint per-monomer with the
   separator character honoured.

Expected:

- After step 1, the HELM table is open and the sequence column
  carries tags `quality=Macromolecule, units=helm`. The
  `bilnSequenceCellRenderer` is the renderer for the sequence
  column (visible via the per-monomer painted cells; no plain-
  string fallback rendering).
- After step 3, the MSA table is open and the sequence column
  carries tags `units=separator, .separator='-'`. The
  `separatorSequenceCellRenderer` is the renderer for the
  sequence column (cells show separator-divided monomers, not
  the BILN-style continuous letter run).
- No error balloon appears (`isErrorBallon` returns false on
  both opens).
- The `bio.rendering` base surface is exercised (top-level cell-
  renderer registration list is reachable; the test runner can
  enumerate Macromolecule renderers via the grok cell-renderer
  registry).

### Scenario 2 — Convert HELM -> SEPARATOR re-dispatches renderer (GROK-12164 guard)

Steps:

1. Open `System.AppData/Bio/tests/filter_HELM.csv` from the
   Files browser. Detector sets `units=helm`; cells render via
   `bilnSequenceCellRenderer`.
2. Right-click the sequence column header. On the column
   context menu, run
   **Bio > Transform > Convert Sequence Notation...**. When the
   convert dialog opens, set target notation to **Separator**
   with separator string `-` (the GROK-12164 reproducer
   parameters), then click **OK**.
3. The convert action produces a new sequence column whose tags
   are `units=separator, .separator='-'`.
4. Confirm the new column renders via the
   `separatorSequenceCellRenderer`, not the
   `bilnSequenceCellRenderer` that was active on the source
   HELM column. This is exactly the failure mode reported in
   GROK-12164: pre-fix, the cell renderer dispatch did not
   follow the new `units` tag and the new column was painted
   by the bilnSequenceCellRenderer.

Expected:

- After step 2 the convert dialog closes cleanly and a new
  column appears on the dataframe with `units=separator`.
- The new column's painted cells show the separator character
  between monomers (separator dispatch). The original HELM
  column continues to render via the BILN dispatch.
- The post-convert renderer dispatch follows the new column's
  `units` tag, not the source column's tag — this is the
  regression-guard contract of GROK-12164.

### Scenario 3 — Split-to-Monomers produces Monomer columns rendered by monomerCellRenderer

Steps:

1. Open `System.AppData/Bio/tests/filter_MSA.csv` from the
   Files browser (separator-with-MSA-flag dataset; detector
   sets `units=separator, aligned=SEQ.MSA`).
2. On the menu ribbon, open
   **Bio > Transform > Split to Monomers...**. When the dialog
   opens prefilled with the active table's Macromolecule column,
   click **OK**.
3. The transform produces N per-position columns named
   `Monomer (1)`..`Monomer (N)` (N = alignment width of the
   source column). Each Monomer column carries
   `quality=Monomer`.
4. Confirm the per-position Monomer columns render via
   `monomerCellRenderer`. The renderer paints a single monomer
   per cell with the monomer-library color and label.

Expected:

- After step 2 the dialog closes cleanly and N new Monomer
  columns appear on the dataframe (N = alignment width;
  visible as adjacent grid columns).
- Each new column has `quality=Monomer` and paints a single
  monomer letter per cell with the monomer-library color
  palette (`monomerCellRenderer`, not the parent sequence
  renderer family).
- No error balloon appears.

### Scenario 4 — Custom-notation column dispatches to customSequenceCellRenderer

Steps:

1. Open `System.AppData/Bio/tests/filter_HELM.csv` from the
   Files browser. Detector sets `units=helm`.
2. Programmatically set `units=custom` on the sequence column
   to simulate a pluggable notation provider (this registration
   is keyed on the units value `custom`, and exists for the
   pluggable-notation extension point). The runner uses
   `column.setTag('units', 'custom')` and triggers a grid
   invalidate via `view.grid.invalidate()`.
3. Confirm the `customSequenceCellRenderer` becomes the
   dispatch target. Because the pluggable-notation provider is
   not registered on the default Bio runtime, the renderer
   falls back to a string-painted display of the underlying
   sequence string rather than the monomer-segmented BILN
   layout used pre-step-2.

Expected:

- After step 2 the sequence column tag set is
  `units=custom` and `view.grid.invalidate()` re-dispatches
  the renderer.
- The painted display differs from the pre-step-2 BILN-style
  paint (the customSequenceCellRenderer registration is
  exercised even when no pluggable provider is wired in;
  fallback render is the contract).
- No error balloon appears.

## Notes

- This scenario doesn't cover the macromolecule-difference (diff)
  renderer — that one's correctness is a manual/visual check, out
  of scope here.

---
{
  "order": 20
}
