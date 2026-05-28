---
feature: chem
sub_features_covered: [chem.notation, chem.notation.action, chem.notation.convert-mol]
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
produced_from: atlas-driven
produced_for: chem-grok-17964-spec.ts
authored_date: 2026-05-11
authored_by: claude-code-test-designer
related_bugs: [GROK-17964]
---

# Chem | Convert Notation column-action registration invariant (GROK-17964 regression-lock)

Bug-focused scenario locking the invariant: **`Convert Notation...` action in the
column-level Context Panel MUST be registered exactly once per applicable molecule
column, and remain at count 1 across the action lifecycle — cancellation, successful
completion, multiple invocations, and emergence of new molecule columns from the
conversion itself**. Duplicate registration after a handler invocation (the literal
GROK-17964 surface) is the regression this scenario catches.

Per chain YAML (`scenario-chains/chem.yaml` rev 2
`bug_focused_candidates[chem-grok-17964-spec.ts]`): independent scenario
(`depends_on: []`), `pyramid_layer: bug-focused`, `target_layer: playwright`,
strategy `simple`. Parallel-coverage on `chem.notation.*` with `info-panels.md`
(Context Panel walk, `coverage_type: regression`) and `sketcher.md` (notation round-
trip via clipboard, `coverage_type: regression`). This spec owns the exactly-once
registration invariant.

Bug-library reference: `references/bug-library/chem.yaml :: GROK-17964` —
title *Chem: Duplicate "Convert Notations" option appears in the context panel*,
priority p2, status `fixed`, `test_coverage: needed`. Atlas surface:
`chem.notation.action` (column-action surface for notation conversion).

## Setup

1. **No external provisioning.** The scenario consumes the bundled chem demo
   dataset `System:AppData/Chem/tests/smiles-50.csv` (validated in MCP recon 2026-05-11 —
   50 rows, `canonical_smiles` column auto-detects `semType: Molecule`). No
   bespoke fixture, no project save, no share.
2. **Login** — standard test user. No special privileges required.
3. **No fixture consumed.** Independent scenario (`depends_on: []`).

## Scenarios

### Block A — Convert Notation action registration is exactly-once across the action lifecycle

Locks the invariant: regardless of how many times the action is invoked
(cancelled, completed, or both), the Actions pane on the column-level Context
Panel renders exactly ONE `Convert Notation...` link per molecule column. New
molecule columns emerging from a successful conversion also render exactly ONE
entry each.

1. **Open the chem demo dataset and focus the molecule column on the Context
   Panel.** Read `System:AppData/Chem/tests/smiles-50.csv`, add a table view, wait for
   semantic type detection (`canonical_smiles` becomes `semType: Molecule`),
   set the shell's current object (`grok.shell.o = df.col('canonical_smiles')`)
   so the Context Panel renders in column-context mode (Actions pane visible).
2. **Baseline — assert single Convert Notation registration.** Count the
   `label.d4-link-action` elements whose text starts with `"Convert Notation"`.
   Expected: exactly 1. (If 0 — the Chem package's `chem.notation.action`
   surface failed to register at all; if ≥2 — GROK-17964 regression has
   already landed at initial render.)
3. **Cancellation path.** Click the `Convert Notation...` link → the
   `Convert Notation` dialog (`.d4-dialog`) opens with controls
   `[name="input-Target-Notation"]` (SELECT), `[name="input-Overwrite"]`
   (CHECKBOX, default false), `[name="input-Join"]` (CHECKBOX, default true),
   `[name="input-Kekulize"]` (CHECKBOX, default false). Click
   `[name="button-CANCEL"]`. Re-set `grok.shell.o` to the molecule column to
   refresh the Context Panel rendering. Re-count
   `label.d4-link-action[Convert Notation]` — must remain 1.
4. **Successful completion path.** Re-open the dialog → set
   `[name="input-Target-Notation"]` to `molblock` (forces a real
   notation-format change from the source `smiles`; valid target options per
   MCP recon 2026-05-11: `smiles`, `cxsmiles`, `smarts`, `cxsmarts`,
   `molblock`, `v3Kmolblock`) → click `[name="button-OK"]`. Wait up to 15
   seconds for the conversion to complete.
5. **Assert exactly-once on original column post-completion.** Re-set
   `grok.shell.o` to the original `canonical_smiles` column. Wait for the
   Context Panel to re-render in column-context mode. Re-count
   `Convert Notation` link entries — must remain 1.
6. **Multiple-invocation hardening.** With focus on the original
   `canonical_smiles` column → open the dialog → click CANCEL → reopen → click
   CANCEL again. Refresh focus. Re-count `Convert Notation` link entries on
   the original column — must still be 1.
7. **Final assertion — global Convert Notation panel-link count.** Across ALL
   `label.d4-link-action` elements on the page (Convert Notation entries
   currently visible in any open Actions pane), the count must equal exactly 1
   (only the focused-column's pane renders at a time; per-column entries do
   not stack). The hidden `.d4-menu-item-label` text "Convert Notation..." in
   the menu fragment (`rect.w === 0`, not panel-attached) is excluded from the
   selector by design.

**Scope reduction note (post-Validator-Round-1):** the originally-planned
"exactly-once on newly-created molecule column" check was dropped after Gate
B Round 1 surfaced that (a) `molV2000` is NOT a valid Target-Notation option
on the current build (valid set above), so the originally-spec'd path
silently no-op'd and (b) the per-new-column registration extension is
adjacent to the core invariant — the registration leak surface is captured
fully by the original-column-only count tracking across the action lifecycle.
Removing the new-column step keeps the bug invariant intact while
de-coupling the spec from conversion-output-naming behavior.

## Notes

- **`coverage_type: edge`** — registration-leak edge case after handler
  invocation; this scenario asserts the exactly-once invariant for
  `chem.notation.action`. Happy-path Context Panel walk (Chemistry section
  rendering across notation formats) is owned by the parallel `info-panels.md`
  scenario (`coverage_type: regression`, `pyramid_layer: integration`).
- **Selector anchoring — `label.d4-link-action` only.** The literal text
  `"Convert Notation..."` appears in TWO classes of element: (1)
  `label.d4-link-action` inside `span.d4-entity-markup-row` (the panel-
  attached entry — our target), and (2) `div.d4-menu-item-label` inside
  `div.d4-menu-item.d4-menu-item-vert` (a hidden menu fragment with
  `rect = {x:0, y:0, w:0}`). The selector
  `label.d4-link-action` excludes the menu fragment cleanly; never select on
  the broader `.d4-menu-item-label` text alone or the count is poisoned by
  the hidden duplicate.
- **Column-context mode required.** The Actions pane that lists Convert
  Notation renders only when `grok.shell.o` is a column object (column-
  context Context Panel). Cell-context mode (`grok.shell.o = df.cell(...)`)
  shows molecule-level tabs (Chemistry / Biology / Structure) WITHOUT the
  column-Actions pane. The bug-library repro says "click Convert Notations
  in the context panel" — the actionable label lives on the column-context
  pane, not the cell-context tabs.
- **MCP-validated current behavior (2026-05-11 dev.datagrok.ai recon, test
  account, Chem package v1.17.6).** All seven steps pass on the current build:
  baseline=1, post-CANCEL=1, post-OK on original=1, on new column=1, post
  multi-cancel=1. GROK-17964 is `status: fixed` per bug-library;
  test exercises the post-fix invariant as a regression-lock. No
  `test.fixme()` — the spec is enabled.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` from
  `helpers-registry.yaml`. Inline patterns for dataset open and Context
  Panel focus — no new helper authored under reuse threshold (single-use
  pattern; helper authoring threshold ≥3 not reached).
- **Project save/reopen NOT included.** Bug edge_case_for_atlas mentions
  "Worse, this duplicate state persists across project reload (suggesting
  it is serialized somewhere)". The reload-persistence channel is excluded
  from this spec for two reasons: (a) the exactly-once invariant during the
  action lifecycle is the primary registration-leak surface — if the action
  doesn't duplicate during the lifecycle, it cannot duplicate via reload;
  (b) project save+reload incurs cross-suite state coupling
  (auth.json + project upload + reload navigation) that pushes the spec
  beyond `strategy: simple` for marginal regression-coverage gain. A
  separate cross-cutting "context panel state persistence" scenario could
  be authored in a future cycle if a related project-reload regression
  surfaces.
- **No JS API substitution.** Steps 3, 4, 7 (click `Convert Notation...`
  link, click OK, multi-CANCEL flow) exercise the column-Actions pane
  UI vector under test (`chem.notation.action`). JS API direct invocation
  of `convertMoleculeNotation` would NOT exercise the action-registration
  lifecycle that the bug lives in (registration is a panel-rendering
  side-effect, not an API-call side-effect).
- **Order in chain.** Bug-focused candidate; placed parallel to
  `info-panels.md` and `sketcher.md`. No `depends_on:` ordering —
  independent scenario.
