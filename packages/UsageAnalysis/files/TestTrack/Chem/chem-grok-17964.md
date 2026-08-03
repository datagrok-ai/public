---
feature: chem
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
priority: p1
realizes_atlas: [GROK-17964]
realizes: []
realized_as:
  - chem-grok-17964-spec.ts
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

This spec runs parallel to `info-panels.md` (Context Panel walk) and `sketcher.md`
(notation round-trip via clipboard) — it owns specifically the exactly-once
registration invariant.

Bug reference: GROK-17964 — *Chem: Duplicate "Convert Notations" option appears in
the context panel* (priority p2, fixed).

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
   Expected: exactly 1. (If 0 — the Chem package's notation-conversion action
   failed to register at all; if ≥2 — GROK-17964 regression has already
   landed at initial render.)
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

- **Column-context mode required.** The Actions pane that lists Convert Notation renders
  only when the Context Panel is in column-context mode (column selected). Cell-context
  mode (a single molecule cell selected) shows Chemistry / Biology / Structure tabs
  instead, without the column-Actions pane.
- **Why the action must be clicked, not called directly.** Convert Notation registration
  is a panel-rendering side effect, not an API-call side effect — invoking
  `convertMoleculeNotation` via the JS API would not exercise the registration-leak
  surface the bug lives in.
- **Project save/reopen not covered.** The bug report notes the duplicate can persist
  across a project reload, but this spec doesn't test that path — if the action doesn't
  duplicate during a single session's lifecycle, it can't duplicate via reload, and
  project save+reload adds cross-suite state coupling for marginal extra coverage. A
  separate scenario could be authored later if a reload-specific regression surfaces.
