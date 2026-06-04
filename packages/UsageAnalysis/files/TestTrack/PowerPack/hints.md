---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-function-hint-tooltip
ui_coverage_delegated_to: add-new-column.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/hints.md
migration_date: 2026-05-20
source_text_fixes:
  - fucntion-typo-corrected-to-function
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
realized_as:
  - hints-spec.ts
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T09:49:00Z
    spec_runs:
      - spec: hints-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 23
        failure_keys: []
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-hints
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-hints
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
    review_round: 1
---

# Add New Column dialog — function-hint tooltip

Focused UI scenario verifying that hovering over an inserted function
name in the Add New Column dialog's formula editor surfaces a tooltip
containing the function's signature.

`pyramid_layer: ui-smoke` per `scenario-chains/powerpack.yaml` rev 1
(Rule 1 fallback — single ultra-short UI flow over one source, `simple`
classification, no related bugs). NOT the chain's elected smoke witness;
`coverage_type: regression` is used here (the chain reserves
`coverage_type: smoke` for the elected smoke scenario at top-level
`add-new-column.md`, which exercises the dialog UI holistically).

`ui_coverage_delegated_to: add-new-column.md` — this scenario owns the
specialty function-hint-tooltip flow on its responsibility list and
delegates the basic dialog-open-and-add flow to the smoke witness.

`related_bugs: []` — no curated bug in
`bug-library/powerpack.yaml` intersects this scenario's
`sub_features_covered`. (GROK-17109 and GROK-17004 touch the same
sub-features but their reproduction surfaces — formula-recalc-on-rename
and complex-paste-handler crash — are not exercised by the hover/tooltip
flow.)

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining (`depends_on: []` per chain rev 1).

## Scenarios

### Hover over inserted function name → tooltip with signature appears

1. Open the demog dataset (`System:DemoFiles/demog.csv`) — for example
   via **Browse** > **Files**, **File** > **Open**, or the equivalent
   JS-API loader. **Verify:** the demog grid renders.
2. Open the **Add New Column** dialog — click the `Add New Column`
   icon on the table view toolbar (or **Edit** > **Add New Column**
   from the top menu). **Verify:** the Add New Column dialog opens
   with its formula editor (CodeMirror), functions panel, and preview
   grid visible.
3. Insert any function into the formula text field using any
   supported approach (typing the function name and confirming via
   autocomplete, clicking a function entry in the functions panel, or
   dragging a function onto the formula field). For example, type
   `a` and select `Abs` from the autocomplete tooltip; the function
   is inserted as `Abs(num)`.
4. Hover the mouse pointer over the inserted function name (e.g.
   over `Abs` in `Abs(num)`). **Verify:** a tooltip appears
   containing the function's signature (function name plus its input
   parameter list, e.g. `Abs(num)`).

## Notes

- **target_layer: playwright** — chosen because a sibling
  `hints-spec.ts` already exists at the playwright layer (per
  `existing-test-index.yaml` line 33484, layer `playwright`). All
  other Add-New-Column siblings (`autocomplete-spec.ts`,
  `highlight-spec.ts`, `functions-sorting-spec.ts`,
  `input-functions-spec.ts`, `formula-refreshing-spec.ts`,
  `add-new-column-spec.ts`) are also Playwright specs, so the
  house style is consistent.
- **coverage_type: regression** — per the chain's smoke-scenario
  election (`ui_coverage_plan.smoke_scenario: add-new-column.md`),
  the top-level `add-new-column.md` owns the `smoke` slot. This
  narrower per-flow scenario lands in `regression`. (Per the
  migration prompt's `coverage_type` heuristic: `smoke` for one
  fast happy-path test per feature; `regression` is the default
  for the rest.)
- **Source-text fix applied:** the original scenario contains
  the typo `fucntion` twice in Step 4 ("Hover over fucntion name -
  a tooltip with fucntion signature should appear"). The migrated
  body uses the correct spelling `function`. Tracked in
  `source_text_fixes: [fucntion-typo-corrected-to-function]`.
- **Helpers (already in registry, available for downstream
  Automator):** `softStep`, `loginToDatagrok`, `specTestOptions`,
  `stepErrors` from
  `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts` —
  used by the existing `hints-spec.ts`.
- **Bug-library status:** consulted —
  `bug-library/powerpack.yaml` exists and was scanned; no curated
  bug intersects this scenario's `sub_features_covered` (cross-
  cutting candidates GROK-17109 and GROK-17004 surface in
  sibling scenarios `AddNewColumn/add-new-column.md` /
  `AddNewColumn/formula-refreshing.md` and `AddNewColumn/highlight.md`
  respectively, per chain rev 1 `bug_focused_candidates[]`).
- **Decision log status:** queried — no
  `failed_attempts WHERE feature == powerpack` entries that touch
  the `AddNewColumn/hints.md` migration. No retry-skip approaches
  apply.
- **Atlas linkage (`derived_from:` provenance):**
  - `powerpack.dialogs.add-new-column` interaction "type expression
    with caret on a function name → hint div shows function
    signature" is derived from
    `public/packages/PowerPack/src/tests/add-new-column.ts#L68`
    (atlas `interactions[]` entry, test-mined).
  - `powerpack.dialogs.add-new-column-func` interactions "click
    Add New Column toolbar icon → opens AddNewColumnDialog" and
    "Edit | Add New Column top-menu → opens AddNewColumnDialog"
    are code-derived from
    `public/packages/PowerPack/src/package.ts#L405`.

---
{
  "order": 3
}
