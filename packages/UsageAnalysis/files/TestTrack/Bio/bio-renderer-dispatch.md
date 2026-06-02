---
feature: bio
sub_features_covered:
  - bio.rendering.biln
  - bio.rendering.custom
  - bio.rendering.monomer
  - bio.rendering.separator
  - bio.rendering
  - bio.detector
target_layer: playwright
coverage_type: regression
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

Regression scenario for the Macromolecule cell-renderer family
beyond the FASTA path already exercised by `convert.md` and
`composition-analysis.md`. Bio registers four sequence cell
renderers keyed off the `units` column tag set by the detector:
`fastaSequenceCellRenderer` (units=fasta),
`separatorSequenceCellRenderer` (units=separator),
`bilnSequenceCellRenderer` (units=biln), and
`customSequenceCellRenderer` (units=custom) for pluggable
notation providers. A fifth renderer, `monomerCellRenderer`,
renders the per-position `Monomer` columns produced by the
Split-to-Monomers transform.

This scenario verifies that for each of the three remaining
units classes (separator, biln/HELM, custom-via-monomer) the
renderer dispatch follows the column tags set by the
`detectMacromolecule` synchronous classifier, and that the
`monomerCellRenderer` paints the per-position columns produced
downstream of a transform. It is the non-FASTA complement to
the FASTA dispatch already covered by `bio.cp.detect-open-fasta`
(realized in `convert.md` and other scenarios under
`bio.rendering.fasta`) and the HELM dispatch sketched by
`bio.cp.detect-open-helm` (atlas critical_path, `priority: p0`).

GROK-12164 (HELM -> SEPARATOR convert leaves the bilnSequence
renderer dispatched on a now-separator column) is a regression
guard for Scenario 2 — the convert step in that scenario is
exactly the failure trigger from the bug, and the scenario's
expected outcome (post-convert renderer dispatches to
`separatorSequenceCellRenderer`) is the fixed behaviour.

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
   Files browser. The Macromolecule detector
   (atlas `bio.detector`) classifies the sequence column
   synchronously and sets column tags
   `quality=Macromolecule, units=helm`.
2. Confirm the sequence column renders via the HELM dispatch
   path. The `bilnSequenceCellRenderer` registration
   (atlas `bio.rendering.biln`, `package.ts#L344`) is keyed on
   `units=helm` (HELM strings are BILN-shape under the bio
   library notation hierarchy); cells paint monomer letters
   with the bio-library monomer-color palette.
3. Open `System.AppData/Bio/tests/filter_MSA.csv` from the
   Files browser in a second table view. The detector
   classifies the column and sets
   `units=separator, .separator='-', aligned=SEQ.MSA`.
4. Confirm the separator dispatch path. The
   `separatorSequenceCellRenderer` registration
   (atlas `bio.rendering.separator`, `package.ts#L330`) is
   keyed on `units=separator`; cells paint per-monomer with
   the separator character honoured.

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
   **Bio > Transform > Convert Sequence Notation...**
   (atlas `bio.transform.convert-notation.top-menu`,
   `package.ts#L1131`). When the convert dialog opens, set
   target notation to **Separator** with separator string `-`
   (the GROK-12164 reproducer parameters), then click **OK**.
3. The convert action (`bio.transform.convert-notation.action`)
   produces a new sequence column whose tags are
   `units=separator, .separator='-'`.
4. Confirm the new column renders via the
   `separatorSequenceCellRenderer` (atlas
   `bio.rendering.separator`), not the
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
   **Bio > Transform > Split to Monomers...** (atlas
   `bio.transform.split-to-monomers`, `package.ts#L1225`).
   When the dialog opens prefilled with the active table's
   Macromolecule column, click **OK**.
3. The transform produces N per-position columns named
   `Monomer (1)`..`Monomer (N)` (N = alignment width of the
   source column). Each Monomer column carries
   `quality=Monomer`.
4. Confirm the per-position Monomer columns render via
   `monomerCellRenderer` (atlas `bio.rendering.monomer`,
   `package.ts#L218`). The renderer paints a single monomer
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
   to simulate a pluggable notation provider (the
   `bio.rendering.custom` registration is keyed on units value
   `custom` and exists for the pluggable-notation extension
   point — see atlas
   `bio.rendering.custom`, `package.ts#L302`). The runner uses
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

- atlas entry derived from sub_features
  `bio.rendering.{biln, separator, monomer, custom}` and
  parent `bio.rendering` + `bio.detector`; the convert step in
  Scenario 2 derives from atlas critical_path
  `bio.cp.convert-notation-roundtrip` (`priority: p0`) and
  atlas edge_case GROK-12164.
- target_layer rationale: `playwright` — the cell-renderer
  dispatch contract is observable only by inspecting painted
  cells in the grid (the renderer registry is keyed off
  column tags); apitest can confirm the tag transitions on
  the column but cannot confirm which renderer paints the
  cells. The convert top-menu in Scenario 2 and the split-to-
  monomers top-menu in Scenario 3 are ribbon-driven; the
  custom-notation tag injection in Scenario 4 can be staged
  via `column.setTag` from a Playwright `evaluate`. All four
  scenarios are UI-driven and assert against the painted
  grid surface.
- coverage_type rationale: `regression` — the renderer family
  is the long-standing dispatch contract under which the
  GROK-12164 regression fired (renderer dispatch missing the
  post-convert tag change). Three of the four scenarios are
  positive-path coverage of dispatch correctness; Scenario 2
  is the explicit regression guard for GROK-12164. Atlas
  edge_case GROK-12164 carries `coverage_type: regression`;
  the scenario frontmatter mirrors it.
- Sub-features covered:
  - `bio.rendering.biln` (`package.ts#L344`) — HELM/BILN cell
    renderer; covered by Scenarios 1 + 2 pre-convert.
  - `bio.rendering.separator` (`package.ts#L330`) — separator
    cell renderer; covered by Scenarios 1 + 2 post-convert.
  - `bio.rendering.monomer` (`package.ts#L218`) — per-position
    monomer cell renderer; covered by Scenario 3.
  - `bio.rendering.custom` (`package.ts#L302`) — pluggable
    notation cell renderer; covered by Scenario 4.
  - `bio.rendering` (`utils/cell-renderer.ts#L1`) — parent
    cell-renderer surface; exercised across all four scenarios
    via the dispatch enumeration in Scenario 1 step 2.
  - `bio.detector` (`detectors.js#L1`) — Macromolecule
    classification on dataset open; gate condition for every
    scenario's renderer dispatch.
- Bug-context (`related_bugs`):
  - **GROK-12164** (HELM -> SEPARATOR convert leaves the
    bilnSequenceCellRenderer dispatched on the new
    `units=separator` column). Atlas edge_case
    `coverage_type: regression`; Scenario 2 is the explicit
    regression guard. This scenario is the bug-aware-tie-break
    pick from the breadth-extension net-new pool: each of the
    four `bio.rendering.*` sub_features carries net-new
    coverage, and `bio.rendering.separator` is GROK-12164's
    primary affected surface. Per STEP C bug-aware tie-break
    (this skill's net-new refusal section), the renderer
    dispatch scenario yields double value (breadth + bug-
    repro adjacency) over a same-net-new annotate-only
    scenario.
- Net-new ids vs `live_covered_union` (this cycle's
  authoritative covered set):
  - `bio.rendering.biln` — net-new.
  - `bio.rendering.custom` — net-new.
  - `bio.rendering.monomer` — net-new.
  - `bio.rendering.separator` — net-new.
  - `bio.rendering`, `bio.detector` — already in
    `live_covered_union` (anchors, not net-new).
  - net_new = 4 ids; satisfies the SR-loop progress-sensitive
    bound (delta > 0) and the STEP C net-new refusal.
- Manual-only subset: none of the six covered sub_features
  appear in atlas `manual_only[]` (verified against atlas rev
  3 `manual_only[]` list — `bio.viewers.web-logo`,
  `bio.viewers.vd-regions`, `bio.rendering.column-header`,
  `bio.rendering.macromolecule-difference`, the demo entries,
  `bio.panels.{structure-3d, atomic-level, tooltip}`). In
  particular `bio.rendering.macromolecule-difference` IS
  manual_only and is therefore excluded from the renderer
  dispatch scope of this scenario; the difference renderer is
  out of scope here.
- Deferrals: none. All four scenarios are UI-driven and
  observable via Playwright.
- The section has no Bio/grok-browser ref-doc verb-form H2
  matching the citation regex, so `## Notes` citation
  pointers reference atlas entries only — no
  `See: <path>#<heading>` citation form applies here.
- This scenario covers 6 sub_features
  (`F-STRUCT-DENSITY-01` floor: 2;
  `F-STRUCT-INTERACTION-01` floor: 3 in a multi-sub_feature
  scenario — satisfied).

---
{
  "order": 20
}
