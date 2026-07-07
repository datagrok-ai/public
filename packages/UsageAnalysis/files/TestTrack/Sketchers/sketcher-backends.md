---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [chem.cp.sketcher-open-set-readback, chem.cp.sketcher-backend-switch]
produced_from: atlas-driven
pyramid_layer: bug-focused
migration_date: 2026-06-01
original_path: public/packages/UsageAnalysis/files/TestTrack/Sketchers/sketcher-backends.md
authored_date: 2026-06-01
realized_as: [ sketcher-backends-spec.ts ]
related_bugs:
  [
    GROK-16340,
    GROK-12685,
    GROK-12391,
    GROK-12297,
    GROK-12966,
    GROK-14028,
    GROK-12826,
    GROK-12581,
    GROK-12905,
    GROK-12758,
    CLAUDE-5
  ]
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-03-sketchers-migrate-03
    timestamp: 2026-06-03T23:15:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-sketchers-migrate-03
    timestamp: 2026-06-03T22:43:00Z
    spec_runs:
      - spec: sketcher-backends-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 212
        failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-03-sketchers-migrate-03
    timestamp: 2026-06-03T20:15:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-03-sketchers-migrate-03
    timestamp: 2026-06-03T16:00:00Z
    failure_keys: []
---

# Chem | Sketcher — battery across OpenChemLib → Ketcher → ChemDraw (UI backend switch)

Runs the same battery of functional checks (C1-C8 below) against each molecule sketcher backend —
OpenChemLib, Ketcher, and ChemDraw — on a single sketcher panel, switching backends manually
through the hamburger ≡ menu the way a user would. Each check in the battery guards a specific,
previously reported bug, so every backend gets identical coverage.

A sketcher backend is one of the molecule-drawing widgets available from the sketcher's hamburger
menu (`OpenChemLib / Ketcher / Marvin / ChemDraw`). This scenario covers OpenChemLib, Ketcher, and
ChemDraw (Marvin is excluded — it returns an empty molecule on dev without a license).

## The battery (per backend)

| # | Check | Bug it guards |
|---|-------|---------------|
| C1 | `setMolecule('c1ccccc1')` → `getSmiles()` is benzene; input field is not the literal string `'undefined'` | GROK-12685 (copy-as-SMILES → 'undefined') |
| C2 | After set, the sketcher widget persists ~2s: inner widget alive, root still in DOM, molecule still present (does **not** vanish) | GROK-16340 (sketcher disappears after creation) |
| C3 | MOLBLOCK (V2000) read back then re-applied → still a molecule | baseline round-trip |
| C4 | A large **V3000** molblock loads and is not empty | GROK-12391 (OCL doesn't open molV3000) |
| C5 | `setMolecule(SMARTS, substructure=true)` → `getSmarts()` round-trips | GROK-12297 / GROK-12966 (SMARTS) |
| C6 | A malformed molecule string does **not** throw | CLAUDE-5 (removeDuplicates throws on malformed) |
| C7 | `resize()` does **not** throw | GROK-12826 (incorrect behavior on resize) |
| C8 | `setMolecule('')` clears: `isEmpty()` and the molecular input field is empty | GROK-14028 (Reset doesn't clear input) |
| — | After the menu switch the active backend == the selected one | GROK-12581 / GROK-12905 (switch doesn't take / sketcher breaks after switch) |
| — | Zero sketcher `console.error` during each backend's battery | GROK-12758 (console errors) |

## Setup

1. **Chem package loaded** (registers the sketcher backends, enables `DG.chem.Sketcher`). No
   dataset needed — the battery drives the widget directly (molecule strings in/out). The C4
   V3000 data is an **embedded molblock constant** (no file, no server `convert` call — an
   unbounded `convert` was observed to hang the whole test under load).
2. **No fixture consumed.** Independent scenario.
3. **Backends present on the target env:** OpenChemLib (always), Ketcher, ChemDraw. The scenario
   selects each by name through the menu.

## Scenario

1. **Open the sketcher.** A `DG.chem.Sketcher` is shown in a dialog (a handle is kept so the
   battery can read molecule state back). The initial backend is whatever the user last selected
   (Chem init restores it from settings), so it is **not** assumed.
2. **Select OpenChemLib via the hamburger menu → battery.** Click ≡, click **OpenChemLib**, wait
   for the switch, run C1-C8. All assert **hard**.
3. **Select Ketcher via the hamburger menu → battery.** Click ≡, click **Ketcher**, wait for the
   external widget to load, run the same C1-C8 on the same panel. All assert **hard**.
4. **Select ChemDraw via the hamburger menu → battery.** Click ≡, click **ChemDraw**, wait, run
   the same C1-C8. All assert **hard**.

Each step also asserts the active backend equals the one just selected (GROK-12581 / GROK-12905)
and zero sketcher console errors (GROK-12758).

**Pass:** all three batteries (OpenChemLib, Ketcher, ChemDraw) pass in full. Verified stable
across 3 consecutive Playwright runs (Gate B: 3/3).

## Notes

- **Manual switch, not JS-API switch.** Backend selection is a real DOM click on the hamburger
  menu radio item (`.d4-menu-item-label` with the backend name), matching the manual gesture —
  not a programmatic `sketcherType` assignment. Molecule state is *read back* via the kept widget
  handle (verification, not action).
- **One persistent widget is the key to stability.** Creating a fresh sketcher per backend
  flakily **crashed the headless renderer** (~50% of runs) — an uncatchable Playwright fatal.
  Switching backend **in place on a single panel** via the menu is stable (3/3 green headless).
- **Why JS-API for the checks.** Neither the OCL canvas nor the external iframes (Ketcher /
  ChemDraw) can be reliably "drawn into" with Playwright; `setMolecule` / `getSmiles` /
  `getMolFile` / `getSmarts` / `resize` is the robust functional probe (same reasoning as
  `sketcher-spec.ts`).
- **Console-error scope.** The zero-console-error assertion (GROK-12758) counts only
  **sketcher-relevant** errors (filtered by sketcher/search-pattern/copy-as/molfile patterns);
  unrelated platform/login noise is excluded, mirroring `sketcher-spec.ts`.
- **Starting-backend determinism.** Chem init restores the last-selected sketcher from
  userSettings (`package.ts:188-197`), so the dialog may open on any backend. The OCL step
  therefore selects OpenChemLib via the menu like the others; the spec restores
  `userSettings 'sketcher'/'selected' = OpenChemLib` on teardown.
- **Out of scope here (cross-context UI — Filter Panel / Context Pane).** These bugs need the
  Filter Panel or Context Pane wired to the sketcher, not the standalone widget — candidates for
  a dedicated filter-panel spec, flagged not dropped:
  - GROK-12903 — filter synchronization (filter panel ↔ hamburger after backend switch)
  - GROK-12905 — sketcher won't reopen from the **filter panel** after switching back (the
    in-place re-switch is covered here; the filter-panel reopen path is not)
  - GROK-12581 — switching on the Filter Panel doesn't update the hamburger / Context Pane
  - GROK-12505 — `.structure-filter-type = Categorical` should yield a categorical filter
  - GROK-12803 — editing a DB query parameter via the Sketcher on the Context Pane
- **Sibling specs.** `sketcher.md` / `sketcher-spec.ts` cover the OCL cell-editor walk
  (Favorites / Recent / Copy / Paste); `sketcher-ui.md` carries the #2448 stereo-on-highlight
  manual invariant. This scenario is the cross-backend functional battery.
