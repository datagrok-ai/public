# Migration Report — sketcher.md (split into sketcher.md + sketcher-ui.md)

First-cycle migration of `sketcher.md` (TestTrack/Chem section) per chain YAML
`scenario-chains/chem.yaml` rev 2 (`output_plan.sketcher.md`: `target_layer: playwright`,
strategy `simple`, classification `medium`, `pyramid_layer: bug-focused`,
`depends_on: []`). Chain rev 2 footer note (d) (Olena 2026-05-11) directs the Migrator to
**split** the original scenario into two outputs:

1. **`sketcher.md`** (this report's primary output) — canonical sketcher cell-editor walk
   over Favorites / Recent / Copy / Paste / backend enumeration. `coverage_type: regression`,
   `target_layer: playwright`. The original's `pyramid_layer: bug-focused` chain annotation
   becomes inert after the source-text transformations (#1608 drop + #2448 split) — what
   remains is multi-action regression-of-the-set.
2. **`sketcher-ui.md`** (new, `produced_from: decomposed`) — ui-only manual smoke for the
   #2448 stereochemistry-preserved-on-SMILES-highlight visual invariant.
   `coverage_type: smoke`, `target_layer: manual`, `pyramid_layer: ui-smoke`.

## Step mapping

The original is a 31-line, 13-numbered-step scenario (with a step-numbering glitch
`1.` / `2.` / `3.` / `4.` / `5.` / `6.` / `7.` / `8.` / `10.` / `11.` / `9.` / `12.` / `13.`
— gap at 9, doubled-out-of-order 11/9) plus a 2-paragraph "Check:" block carrying two
GitHub-issue regressions (#1608 + #2448) plus a JSON footer
`{ "order": 5, "datasets": ["System:DemoFiles/chem/smiles.csv"] }`. Body defects per chain
rev 2 directive (Olena 2026-05-11, chain footer note (c) for typo + step renumbering and
footer note (d) for the #1608 drop + #2448 split).

| Original step / block | Migrated output | Migrated location | Decision |
|-----------------------|-----------------|-------------------|----------|
| Step 1 — "Open smiles. csv (Browse > Files > Demo Files/chem)" | sketcher.md | Setup step 1 (provisioning — `System:DemoFiles/chem/smiles.csv`) + Scenarios > "Sketcher cell-editor walk on OCL backend" step 1 (dataset open) | preserved (typo "smiles. csv" normalized to `smiles.csv`; split for clarity between provisioning declaration in Setup and the open action in scenario body) |
| Step 2 — "Double-click a molecule." | sketcher.md | Scenarios step 2 (double-click any molecule cell → sketcher modal opens) | preserved |
| Step 3 — "In the hamburger menu of the sketcher, click Favorites > Add to Favorites." | sketcher.md | Scenarios step 3 (Favorites add) | preserved |
| Step 4 — "Enter C1CCCCC1 (not using ctrl+V), press Enter, click OK, sketcher closes" | sketcher.md | Scenarios step 4 (type-directly SMILES C1CCCCC1 → Enter → OK; sketcher closes; underlying cell updated) | preserved (explicit "not via Ctrl+V paste" guard kept verbatim) |
| Step 5 — "Open sketcher again by double clicking the structure. In the hamburger menu, check Recent and Favorites." | sketcher.md | Scenarios step 5 (reopen + inspect Recent + Favorites lists; verify Recent contains cyclohexane from step 4 and Favorites contains the original molecule from step 3) | preserved (explicit verification of both lists made concrete) |
| Step 6 — "In the hamburger menu, click Copy as SMILES." | sketcher.md | Scenarios step 6 (Copy as SMILES → clipboard receives the SMILES) | preserved |
| Step 7 — "Change the molecule" | sketcher.md | Scenarios step 7 (modify the sketch canvas so it differs from the clipboard) | preserved |
| Step 8 — "Go to molecular input field and press CTRL+V. Press Enter. Copied molecule should be displayed" | sketcher.md | Scenarios step 8 (Ctrl+V into molecular input field → Enter → canvas shows the SMILES copied at step 6) | preserved (CTRL+V normalized to canonical Ctrl+V; paste round-trip verification made explicit) |
| Original step "10." (out-of-order; should be 9) — "In the hamburger menu, click Copy as MOLBLOCK." | sketcher.md | Scenarios step 9 (Copy as MOLBLOCK → clipboard receives the MOLBLOCK) | preserved (renumbered cleanly from out-of-order `10.` → `9.`) |
| Original step "11. 9." (collision + out-of-order) — "Change the molecule" | sketcher.md | Scenarios step 10 (modify the sketch canvas again) | preserved (renumbered from collided `11. 9.` → `10.`) |
| Step 12 — "Go to molecular input field and press CTRL+V. Copied molecule should be displayed" | sketcher.md | Scenarios step 11 (Ctrl+V → MOLBLOCK paste round-trip verification → OK; sketcher closes; cell updated) | preserved (CTRL+V → Ctrl+V; renumbered `12.` → `11.`; OK click added as the natural closing step per the parallel earlier Step 4 OK click) |
| Step 13 — "Repeat steps 2-9 for all available sketcher types for selection from hamburger menu" | sketcher.md | Scenarios > "Backend enumeration (opportunistic)" (OCL required as the full block 1-11; Ketcher / Marvin / ChemDraw opportunistic per qa-pw env availability) | preserved as enumeration loop (hard-coded OCL as required; others opportunistic per chain rev 2 directive footer note (d) "hard-code OCL as the only required backend") |
| "Check:" block paragraph #1608 — large-substructure tooltip in Filter Panel sketch box (with verbatim peptide-like SMILES) | (dropped entirely; not propagated) | n/a | **dropped** per chain rev 2 directive footer note (d): "#1608 → IGNORE (drop from scenario body during migration; do not propagate, no bug-library entry)". No #1608 entry added to `bug-library/chem.yaml`. No related-bug frontmatter citation. The verbatim peptide-like SMILES is not preserved anywhere. |
| "Check:" block paragraph #2448 — stereochemistry-preserved-on-SMILES-highlight invariant on `SMILES_highlighted.csv` via Context panel > Cheminformatics > Highlights > Sketch + scaffold tree / `isosmiles` structure filter | sketcher-ui.md | Scenario > "#2448 — Stereochemistry preserved on SMILES highlight (ui-only)" (entire block) | **split** per chain rev 2 directive footer note (d): "#2448 → SPLIT into a new ui-only scenario sketcher-ui.md (target_layer: manual, pyramid_layer: ui-smoke)". `coverage_type: smoke` (visual-fidelity invariant; not edge — though stereochemistry is arguably a fidelity edge case, Olena's directive judgment-call sets `smoke` for the manual-layer ui-smoke placement). |
| JSON footer `{ "order": 5, "datasets": ["System:DemoFiles/chem/smiles.csv"] }` | sketcher.md (Setup step 1) | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration). The `SMILES_highlighted.csv` dataset that was inlined in the original "Check:" #2448 paragraph is provisioned via `sketcher-ui.md` Setup step 1 (preferred path `System:AppData/Chem/tests/SMILES_highlighted.csv` + fallback provisioning instruction). |

No original step is silently dropped without an explicit directive citation. The #1608
paragraph is dropped per explicit chain rev 2 directive footer note (d); the #2448
paragraph is split into `sketcher-ui.md` per the same directive. Step-numbering
collision + out-of-order (original had `1.` / `2.` / `3.` / `4.` / `5.` / `6.` / `7.` /
`8.` / `10.` / `11. 9.` / `12.` / `13.`) renumbered cleanly into the migrated
`sketcher.md` body's `Setup` + `Scenarios > "Sketcher cell-editor walk on OCL backend"` +
`Scenarios > "Backend enumeration (opportunistic)"` linear sequence.

## Decisions

- **Two-output split per chain rev 2 directive (Olena 2026-05-11 chain footer note (d)).**
  The original `sketcher.md` carried two distinct concerns: (a) a canonical sketcher cell-
  editor walk (steps 1-12) plus a backend-enumeration loop (step 13), AND (b) two
  GitHub-issue regression reproductions in a trailing "Check:" block (#1608 + #2448).
  Chain rev 2 directive resolves these via:
  - **#1608 drop**: rationale per Olena 2026-05-11 — "#1608 (large-substructure tooltip):
    IGNORE (drop from scenario body during migration; do not propagate, no bug-library
    entry)". The bug surface (Filter Panel sketch box tooltip rendering for large
    structures) is a Filter Panel concern, not strictly a sketcher concern; the existing
    `filter-panel.md` is the natural home if a Filter Panel tooltip invariant is later
    deemed test-worthy. Per Olena's directive, no propagation: no bug-library entry, no
    related_bugs citation, no atlas curator candidate flag generated by this migration.
  - **#2448 split**: rationale per Olena 2026-05-11 — "#2448 (stereochemistry on SMILES
    highlight): SPLIT into a new ui-only scenario `sketcher-ui.md` (target_layer:
    manual, pyramid_layer: ui-smoke)". The invariant is a visual stereochemistry-
    rendering fidelity check; no DOM assertion library can reasonably verify "wedge bond
    remains visible" or "E/Z marking remains visible" without a pixel-diff or vector-
    render comparison harness that the section does not currently have. Manual-layer
    ui-smoke is the natural fit. `coverage_type: smoke` selected (judgment call per
    prompt: this is a stereochemistry-visual-fidelity invariant; could arguably be
    `edge` but Olena's directive footer suggests smoke or edge — manual-layer ui-smoke
    placement defaults to `smoke` per section ui-smoke convention; the chain-level
    A-STRUCT-02 (edge/perf at chain) is satisfied by the chain's
    `bug_focused_candidates[]` and the section's `r-group-analysis.md` already, so
    `sketcher-ui.md` doesn't need to carry edge on its own).
- **Why this `target_layer` (sketcher.md):** chose `playwright` per
  `scenario-chains/chem.yaml` `output_plan.sketcher.md.target_layer = playwright`. The
  scenario requires DOM-level UI driving: double-click open + sketcher modal +
  hamburger-menu walks + molecular input field typing + clipboard keyboard shortcuts +
  paste round-trip. Chain YAML reason field: "Bug-focused scenario (#1608 + #2448
  reproduction) layered on canonical sketcher modal walk … Sibling sketcher-spec.ts
  already exists at playwright per existing-test-index. Sketcher modal + hamburger menu
  + clipboard round-trip + Highlight panel large-structure tooltip rendering require
  real DOM that playwright handles natively." After the #1608 drop + #2448 split the
  "bug-focused" framing in the chain reason field is partly obsolete, but the playwright
  target_layer remains correct for the residual canonical sketcher walk.
- **Why this `target_layer` (sketcher-ui.md):** chose `manual` per chain rev 2 directive
  footer note (d). Visual stereochemistry-fidelity invariant has no straightforward
  playwright assertion path; manual review is the section convention for ui-smoke visual
  invariants (cf. section sibling `WideSmokeTest/Chem/*-ui.md` files for the manual-layer
  ui-smoke idiom).
- **Why this `coverage_type` (sketcher.md):** chose `regression` per chain rev 2 directive
  footer note (d) and the post-split scenario shape. After dropping #1608 and splitting
  #2448 out, what remains is multi-action regression-of-the-set (Favorites + Recent +
  Copy + Paste + backend enumeration). Not `smoke` (section's smoke is
  `Advanced/scaffold-tree-functions.md`). Not `edge` / `perf` (no specific failure-mode
  invariant or threshold). Not `bug-focused` at the migrated-scenario layer (the
  chain-level `pyramid_layer: bug-focused` annotation is now inert because the bug-
  focused content was extracted).
- **Why this `coverage_type` (sketcher-ui.md):** chose `smoke` per chain rev 2 directive
  footer note (d) and the manual-layer ui-smoke convention. Could be `edge` (the
  stereochemistry-on-highlight invariant is technically a visual-fidelity edge case),
  but Olena's directive flagged this as a judgment call and the manual-layer
  `pyramid_layer: ui-smoke` placement defaults to `coverage_type: smoke`. Chain-level
  A-STRUCT-02 (`coverage_type: edge` or `perf`) is already satisfied chain-wide by 10
  bug-focused candidates and by `r-group-analysis.md`.
- **Why this `strategy`:** `simple` per chain YAML `output_plan.sketcher.md.strategy =
  simple`. Single scenario per output file, no cross-file fixture; chain analyzer
  classified `simple`. Pattern 1 (decomposition) IS triggered here per chain rev 2
  directive footer note (d) — `sketcher-ui.md` is a `produced_from: decomposed` child,
  but the decomposition is **directive-driven**, not classifier-driven. Pattern 2
  (bug-focused slice) does NOT apply at the migrated-scenario layer; chain-level
  `bug_focused_candidates[]` reference sketcher.md as a span pointer but no
  bug-focused-slice migration is performed by this Migrator step.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher-spec.ts` — existing
    playwright-layer test ("Chem: Sketcher", category Chem, `layer: playwright`,
    `covers: full` for `chem.sketcher` per coverage-map). Existing spec covers the
    sketcher modal; Automator will extend / align to the migrated body's Favorites +
    Recent + Copy + Paste + backend-enumeration walk. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis.md` +
    `elemental-analysis-migration-report.md` (same-cycle predecessor) — shape anchor
    for migrated `.md` + report structure (frontmatter, Setup / Scenarios / Notes
    order, per-variant data-driven walk, SR-01 carryforward template). Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels.md` +
    `info-panels-migration-report.md` (same-cycle predecessor) — multi-block / multi-
    surface walk house-style anchor. Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate.md` +
    `calculate-migration-report.md` (same-cycle predecessor) — calculator-style
    section-mate; pattern reference for the step-renumbering silent-fix decision.
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/WideSmokeTest/Chem/info-panels-ui.md`
    + `activity-cliffs-ui.md` + `r-group-analysis-ui.md` — the section's existing
    ui-only manual-layer smoke notes. Used as shape anchor for `sketcher-ui.md`
    (single-block manual scenario with a few numbered steps and a clear visual
    expected-result clause). Note that `sketcher-ui.md` is placed alongside its
    parent `sketcher.md` (per chain rev 2 directive target path
    `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher-ui.md`) rather than
    under WideSmokeTest. Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator uses in spec `beforeEach`.
    - `softStep`
      (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`)
      — per-step soft-fail wrapper; Automator wraps each of the 11 sketcher walk
      steps in `softStep` for granular failure reporting.
    - `closeAllViews`
      (`public/packages/UITests/playwright/scripts/helpers.ts:129`)
      — cleanup at end of scenario.
  - **Candidate helpers (NOT yet in registry — flagged for addition via helpers-
    registry curator; per migration-prompt Helpers discipline §):**
    - `helpers.playwright.chem.openSketcherFromCell(page, rowIndex, columnName)` —
      double-clicks a molecule cell and waits for the sketcher modal to open.
      Surfaced by Scenarios step 2 (sketcher.md).
    - `helpers.playwright.chem.sketcherHamburgerMenuClick(page, itemPath)` — opens
      the sketcher's hamburger menu and clicks the specified menu item path (e.g.
      `["Favorites", "Add to Favorites"]` or `["Copy as SMILES"]`). Surfaced by
      Scenarios steps 3 / 6 / 9 (sketcher.md).
    - `helpers.playwright.chem.sketcherTypeSmilesInInput(page, smiles, pressEnter)` —
      types a SMILES into the sketcher's molecular input field (typed, not pasted)
      and optionally presses Enter. Surfaced by Scenarios step 4 (sketcher.md).
    - `helpers.playwright.chem.sketcherPasteFromClipboard(page, pressEnter)` —
      focuses the molecular input field and triggers Ctrl+V then optionally Enter.
      Surfaced by Scenarios steps 8 / 11 (sketcher.md).
    - `helpers.playwright.chem.sketcherBackendSwitcherSelect(page, backendLabel)` —
      opens the sketcher hamburger > "Sketcher type" (or equivalent) and selects
      the named backend, returning a boolean indicating whether the backend was
      available. Surfaced by Scenarios > "Backend enumeration" (sketcher.md).
- **Bug library consulted:** yes — `bug-library/chem.yaml` revision 1 (12 curated_bugs,
  generated 2026-05-05). Grepped `curated_bugs[]` for "sketcher" / "1608" / "2448" —
  zero matches. Neither #1608 nor #2448 is in the curated bug-library, consistent with
  chain rev 1 notes ("they're embedded as scenario-internal regression checks") and the
  chain rev 2 directive footer note (d) explicitly instructs **not** to add #1608 /
  #2448 to the bug-library on migration. `related_bugs: []` in both output frontmatters.
- **Decision log queried:** yes — `decision-log.yaml` grepped for `feature: chem` —
  **zero matches** (consistent with the same-cycle sibling migrations:
  elemental-analysis / info-panels / calculate / scaffold-tree / r-group-analysis).
  This is the FIRST CYCLE for the chem section. No prior `migration_decisions`,
  `layer_decisions`, `manual_only`, or `failed_attempts` entries for chem apply. No
  "approaches off the table" constraints. Additionally grepped for "sketcher" — zero
  matches; no prior sketcher migration decision applies.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).** Per
  `migration-prompt.md` "Cross-cutting bug citations from chain YAML" §: this scenario
  IS a span pointer in three chain `bug_focused_candidates[]` entries:
  - `chem-grok-12758-spec.ts` spans `Advanced/scaffold-tree.md:Step 1` +
    `filter-panel.md:Step 1` + `sketcher.md:Step 2` (cross-cutting Scaffold Tree ×
    Sketcher × Substructure Search). The post-rev-2 sketcher.md Step 2 (double-click
    molecule cell, opens sketcher) is the pointer surface — invariant
    "scaffold-node-open-in-sketcher must not corrupt subsequent searchSubstructure
    state" is NOT exercised by this scenario directly.
  - `chem-grok-14028-spec.ts` spans `filter-panel.md:Step 2` +
    `Advanced/structure-filter.md:Step 1` + `sketcher.md:Step 6` (Filter Panel Reset
    cleanup gap). The post-rev-2 sketcher.md Step 6 (Copy as SMILES — clipboard
    receives SMILES) is the pointer surface; the Reset-clears-BitSet-but-sketcher-
    UI-retains-SMILES invariant is NOT exercised by this scenario directly.
  - `chem-grok-17964-spec.ts` spans `info-panels.md:Step 4` + `sketcher.md:Step 6`
    (Convert Notations duplicate after handler error). The post-rev-2 sketcher.md
    Step 6 (Copy as SMILES — touches Convert Notations action surface) is the
    pointer surface; the action-registration-leak invariant is NOT exercised here
    directly.
  Per migration-prompt these cross-cutting citations are RECOMMENDED, not mandatory;
  F-BUG-COVERAGE-01 at section-complete is authoritative. The dedicated bug-focused
  specs will land later per the chain `bug_focused_candidates[]` plan; none of them is
  the responsibility of this migration step.
- **UI delegation status.** Per chain YAML `ui_coverage_plan.smoke_scenario:
  Advanced/scaffold-tree-functions.md` (sketcher.md is NOT the section smoke) and the
  chain entry `ui_coverage_plan.delegations[].scenario: sketcher.md` /
  `delegated_to: null`, this scenario owns its UI coverage directly (no upstream
  delegation). Chain `ui_coverage_responsibility` lists 10 entries; the two trailing
  entries `chem-info-panel-rendering-highlight-large-tooltip` and
  `chem-info-panel-highlight-scaffold-tree-or-filter` were attached because of the
  #1608 + #2448 sub-blocks. With #1608 dropped, the `-large-tooltip` entry
  effectively transfers to `filter-panel.md`'s ownership if and when the tooltip
  rendering is later test-worthy (no migration action here). With #2448 split,
  `chem-info-panel-highlight-scaffold-tree-or-filter` transfers to `sketcher-ui.md`
  for manual coverage. The remaining 8 sketcher-cell-editor entries are exercised in
  full by the migrated sketcher.md body. No SCOPE_REDUCTION proposal substitutes JS
  API for any UI flow; reaffirmed in the migrated Notes ("No JS API substitution").
- **Scenario constraint extraction (per migration-prompt § Scenario constraint extraction).**
  - **(a) FORBIDDEN substitutions:** "Enter C1CCCCC1 (not using ctrl+V)" in the
    original step 4 is a FORBIDDEN-substitution sentinel: typing the SMILES directly
    is the assertion target; pasting via Ctrl+V would invalidate the molecular-input-
    field typing behavior. Preserved verbatim in sketcher.md Scenarios step 4
    ("type the SMILES directly, **not** via Ctrl+V paste"). `pyramid_layer:
    bug-focused` (chain annotation) implicitly forbids JS API substitution for the
    8 sketcher-cell-editor `ui_coverage_responsibility` flows — surfaced explicitly
    in the migrated Notes (`No JS API substitution`).
  - **(b) REQUIRED actions:** All 8 sketcher-cell-editor entries of
    `ui_coverage_responsibility` are exercised via UI driving in the migrated
    sketcher.md body (double-click open + hamburger Favorites add + Recent inspect
    + Copy as SMILES + paste round-trip + Copy as MOLBLOCK + paste round-trip +
    backend switcher). No flow is deferred at the migrated-scenario layer.
  - **(c) Missing-selector escalation:** Migrator does NOT have explicit selector
    definitions for the sketcher modal's hamburger menu structure or the molecular
    input field in the section's `grok-browser/references/` corpus
    (`references/chem.md` is absent per chain rev 2 notes lines 56-58 / 80-85 —
    Edit 8 ui_consolidation_proposals synthesis cannot run). Spec-time selector
    discovery is the Automator's responsibility (per existing sibling
    `sketcher-spec.ts` patterns + `page.evaluate` introspection); a reference-file
    approval-required proposal to add Sketcher-modal-specific selectors to
    `grok-browser/references/widgets/chem-sketcher.md` is a Phase 2 deliverable —
    out of this migration's scope.
  - **(d) Reference templates:** Code-style anchors for the Automator (cited
    verbatim, not modified here):
    `public/packages/UsageAnalysis/files/TestTrack/Chem/sketcher-spec.ts` (existing
    section-mate spec) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis-spec.ts`
    (section-mate Calculate-menu-style flow as a dialog-walk parallel) +
    `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts`
    (multi-block playwright walk reference).
  - **(e) Scenario authority clause:** Prompt-vs-scenario conflict surfaced and
    resolved by chain rev 2 directive: the original body's #1608 + #2448 blocks
    are removed per Olena's directive (chain footer note (d) is the authoritative
    transformation). No `prompt_scenario_conflict:` flag for orchestrator beyond
    what the chain already encodes.
- **Source-text fixes (silently applied per chain rev 2 directive, Olena 2026-05-11
  chain footer notes (c) + (d)).**
  - **Typo "smiles. csv" → `smiles.csv`** (sketcher.md Setup step 1).
  - **Keyboard-shortcut "CTRL+V" → "Ctrl+V"** (sketcher.md Scenarios steps 8 / 11) —
    canonical form per section playwright-test convention.
  - **Step-numbering glitch `1. 2. 3. 4. 5. 6. 7. 8. 10. 11. 9. 12. 13.`** renumbered
    cleanly to `1.` / `2.` / `3.` … `11.` in sketcher.md's `Scenarios > "Sketcher
    cell-editor walk on OCL backend"` block, with the original step 13 (backend
    enumeration loop) extracted to its own `Scenarios > "Backend enumeration
    (opportunistic)"` block per the chain rev 2 directive's backend-enumeration
    scope clarification.
  - **#1608 paragraph dropped** per chain rev 2 directive footer note (d). No
    propagation; no bug-library entry; no related_bugs citation; no atlas-curator
    candidate flag. The verbatim peptide-like SMILES is not preserved.
  - **#2448 paragraph split** into `sketcher-ui.md` per chain rev 2 directive footer
    note (d) — `target_layer: manual`, `pyramid_layer: ui-smoke`, `coverage_type:
    smoke`. The `SMILES_highlighted.csv` dataset that was inlined in the original
    paragraph is provisioned in `sketcher-ui.md` Setup step 1 (preferred path
    `System:AppData/Chem/tests/SMILES_highlighted.csv` + fallback as a Setup-time
    provisioning requirement). The Context panel > Cheminformatics > Highlights >
    Sketch walk + the scaffold-tree-or-`isosmiles`-filter alternate entry-point is
    preserved verbatim in `sketcher-ui.md` Scenario.
  - **Step 13 backend-enumeration scope hard-coded.** Original step 13 says "Repeat
    steps 2-9 for all available sketcher types". Per chain rev 2 directive footer
    note (d) "hard-code OCL as the only required backend (the fallback per atlas
    chem.sketcher.ocl). … other backends (Ketcher / Marvin / ChemDraw)
    opportunistic per qa-pw env. Each is in a separate package." The migrated
    sketcher.md Scenarios > "Backend enumeration (opportunistic)" reflects this
    scope: OCL block is the pass-blocking required exercise (steps 1-11); other
    backends are pass-extending opportunistic re-runs. Flagged for Automator at
    spec design as `test.describe.parallel`-with-`test.skip` shape per env-
    detected backend availability.
- **Backend enumeration scope (which backends in migrated body).** Per chain rev 2
  directive: **OCL only is required**. Ketcher / Marvin / ChemDraw are opportunistic
  per qa-pw environment availability; each is in a separate plugin package. The
  migrated sketcher.md Scenarios > "Backend enumeration (opportunistic)" encodes
  this scope explicitly (OCL block is the canonical pass-blocking exercise; alternate-
  backend re-runs are pass-extending only). No atlas-curator candidate flag generated
  for non-OCL backends — they are environment-availability concerns, not test-coverage
  gaps.
- **No invented sub_features / helpers.** Every id in `sub_features_covered`
  for both outputs (sketcher.md: `chem.sketcher`, `chem.sketcher.ocl`,
  `chem.sketcher.cell-editor`, `chem.actions.copy-smiles`,
  `chem.actions.copy-molfile-v2000`, `chem.actions.copy-as`; sketcher-ui.md:
  `chem.sketcher`, `chem.sketcher.cell-editor`) is present in `feature-atlas/chem.yaml`
  rev 2 (verified per atlas keys at lines 26909 / 29151 / 31393 / 302675 / 304917 /
  300433). Every reused helper is in `helpers-registry.yaml` (verified per lines
  3543 / 3597 / 3603); candidate helpers are flagged as candidates, NOT referenced
  by name in the migrated bodies.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage) — applies to sketcher.md only

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is edge case or
negative path — i.e. its frontmatter `coverage_type` is `edge` or `perf`") is a
**section-level structural invariant**. After the #1608 drop + #2448 split, the
migrated `sketcher.md` is a canonical sketcher-walk-plus-backend-enumeration regression
scenario — `coverage_type: regression` is the natural fit per chain rev 2
(`pyramid_layer: bug-focused` becomes inert after the bug-focused content is removed;
multi-action regression-of-the-set). Forcing `coverage_type: edge` would mis-classify
(the scenario asserts the canonical hamburger-menu + clipboard round-trip walk across
the OCL backend, not a specific failure-mode invariant such as backend-package-missing
fallback, clipboard-permission-denied behavior, or empty-clipboard paste).

**Resolution path:** A-STRUCT-02 satisfaction is chain-wide via 10
`bug_focused_candidates[]` in `scenario-chains/chem.yaml` rev 2 (three of which span
sketcher.md as a pointer: GROK-12758, GROK-14028, GROK-17964 — each chem-bug-focused
spec naturally carries `coverage_type: edge` once authored). Additionally
`r-group-analysis.md` (same-cycle, same section, `coverage_type: edge`) already
organically satisfies the section-level invariant. Decision-log line 8528-8530 records
the radar.md 2026-05-07 carryforward precedent. Same SR-01 carried by sibling-cycle
reports `elemental-analysis-migration-report.md`, `calculate-migration-report.md`,
`scaffold-tree-functions-migration-report.md`, `scaffold-tree-migration-report.md`,
`info-panels-migration-report.md`, `Advanced/structure-filter-migration-report.md`,
`Advanced/similarity-search-migration-report.md`.

**Scope of waiver:** A-STRUCT-02 verdict on `sketcher.md` is deferred to chain-level
evaluation. Critic A should return SR (carryforward), not FAIL.

### No SR-01 on sketcher-ui.md

`sketcher-ui.md` (`target_layer: manual`, `coverage_type: smoke`) does not require an
A-STRUCT-02 carryforward — A-STRUCT-02 is a structural invariant on the playwright-layer
content authored by the migration cycle; manual-layer scenarios are evaluated against
A-STRUCT-02 differently per the section's manual-layer convention. Additionally chain
rev 2's `r-group-analysis.md` already satisfies the section's edge-coverage need at
chain level; the manual-layer `sketcher-ui.md` does not add to or subtract from that
satisfaction.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through A-STRUCT-06, A-COVERAGE-*,
A-MERIT-*) are expected to PASS on both `sketcher.md` and `sketcher-ui.md` without any
other SR proposals. No step is opted out for effort. The SR-01 entry above cites a
real section-level structural property whose satisfaction path is owned by the chain,
not by this per-scenario migration (D-MERIT-01 compliant)._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet. The migrated
sketcher.md walk (11 steps + opportunistic backend-enumeration loop) is realizable
against the current playwright + helpers-registry + atlas state, given the existing
`sketcher-spec.ts` sibling. The candidate helpers surfaced in Decisions are convenience
abstractions that the Automator may inline pending registration — they do NOT block
spec realization. The `sketcher-ui.md` manual scenario is realizable as a manual run
once the `SMILES_highlighted.csv` dataset is provisioned (Setup step 1 fallback
provisioning instruction). The reference-file write proposal (Decisions § Scenario
constraint extraction (c)) is a downstream Phase 2 enhancement, not a blocker on this
scenario's automation. The #1608 (Filter Panel large-substructure tooltip) surface is
not deferred — it is **dropped** per Olena's directive; if a later cycle revisits it,
the natural home is `filter-panel.md` not `sketcher.md`._

## Edge cases

The original lists no explicit "edge case" keyword (the trailing "Check:" block carried
two GitHub-issue regressions, but those are dropped/split per chain rev 2 directive,
not preserved as edge cases on this scenario). Implicit edge cases derivable from the
post-split scenario content + atlas:

- **Sketcher backend availability variance across qa-pw environments.** OCL is the
  platform fallback that must always be present; Ketcher / Marvin / ChemDraw are
  installed per package availability. The migrated sketcher.md Scenarios > "Backend
  enumeration (opportunistic)" handles this via opportunistic re-runs (pass-extending
  only). PRESERVED as scenario block.
- **Clipboard-permission denied in test browser context.** The Copy as SMILES + paste
  round-trip and Copy as MOLBLOCK + paste round-trip require clipboard read+write
  permissions in the playwright context. PRESERVED as Setup step 3 ("Confirm
  clipboard access is granted"). The failure mode if clipboard access is denied is
  out of scope (atlas-curator candidate: "Sketcher Copy-as cycle when clipboard
  permission is denied: expected behavior is a clear user-visible error, not a
  silent paste failure" — flagged, not in this migration's scope).
- **Empty Recent / empty Favorites on first open.** The migrated body assumes the
  Favorites add at Scenarios step 3 takes effect before the reopen at step 5;
  similarly Recent is populated only after the SMILES edit at step 4. If a previous
  test run polluted the per-user Favorites store, step 5's Favorites assertion may
  see additional entries. PRESERVED as scenario step (Scenarios step 5 — verification
  is "Favorites contains the molecule added at step 3", not "Favorites contains only
  that molecule"). Automator at spec time should clear the per-user Favorites store
  between runs (fixture cleanup); flagged for spec design at Automator stage.
- **Sketcher type-vs-paste behavior for SMILES.** Original step 4 explicitly directs
  "Enter C1CCCCC1 (not using ctrl + V)" — typing the SMILES into the molecular input
  field, not pasting it. The "type, not paste" guard exists because the typing path
  exercises the molecular-input-field SMILES-parsing handler, whereas pasting would
  exercise the clipboard-paste handler (which is exercised separately at steps 8 / 11
  with the Copy-as round-trip). PRESERVED as scenario step (Scenarios step 4 explicit
  "type directly, not via Ctrl+V paste" guard).
- **MOLBLOCK paste behavior when the input field expects SMILES.** Original step 12
  pastes MOLBLOCK (multi-line text) into the molecular input field. The migrated
  Scenarios step 11 preserves this — the input field must accept the MOLBLOCK paste
  and parse it into the same molecule as the pre-edit canvas. If the input field
  rejects MOLBLOCK, the assertion fails. PRESERVED as scenario step.
- **Console errors throughout.** Implicit across the sketcher walk: each hamburger-
  menu click + clipboard interaction + backend switch must complete without console
  errors. The migrated sketcher.md does NOT add an explicit per-step "no console
  errors" verification (sibling `elemental-analysis.md` adds this; here the action
  density is higher and explicit per-step console-error verification would clutter
  the walk). Automator at spec time should add a section-wide `expect(consoleErrors).
  toHaveLength(0)` after the scenario block. Flagged for spec design at Automator
  stage.

For sketcher-ui.md (manual layer):

- **Curated SMILES_highlighted.csv content variance.** The dataset is not bundled
  in the platform's System file shares with verified content per atlas curator
  flag. If the file's exact stereo-feature row distribution differs from the
  invariant's expectation (chiral C atoms, E/Z double bonds, explicit wedge bonds),
  the visual comparison may degrade in coverage. PRESERVED as Setup step 1 fallback
  provisioning instruction.
- **Alternate "make these structures highlighted" entry point.** Original step 5
  of the #2448 paragraph says "either use scaffold tree to highlight or use
  structure filter on 'isosmiles' column". Both are preserved verbatim in
  `sketcher-ui.md` Scenario step 6 as an OR-choice. The manual QA chooses one path;
  flagged as Unresolved ambiguity (the original does not specify which path is
  canonical for the invariant).

No edge case is moved to atlas, manual_only, deferred, or a separate scenario
silently. The optional clipboard-permission-denied / empty-Favorites edge cases are
flagged for atlas curator (no atlas write here).

## Unresolved ambiguities

- **`SMILES_highlighted.csv` dataset path.** Original body included a TODO: "Add to
  linked datasets" — the file is not in the JSON `datasets` list and is not a known
  System path at migration time. `sketcher-ui.md` Setup step 1 documents the
  preferred path `System:AppData/Chem/tests/SMILES_highlighted.csv` (section
  convention for Chem test datasets) and a fallback provisioning instruction (QA
  uploads the file before running the scenario). Flag for QA pair review at first
  manual run; if the preferred path resolves, this becomes a Setup confirmation
  rather than a provisioning step. If the file requires curation (selecting rows
  that exhibit stereochemistry), a follow-up atlas-curator workflow is needed to
  formalize the dataset's specification.
- **Backend-switcher hamburger-menu label.** The migrated sketcher.md Scenarios >
  "Backend enumeration (opportunistic)" refers to the backend switcher as
  "(Sketcher type / equivalent label)" — the exact menu item label is not in the
  Migrator's surface knowledge. Automator at spec time resolves via DOM
  introspection; flag for spec design.
- **Sketcher hamburger-menu Favorites + Recent list selectors.** Spec-time selector
  discovery responsibility (Decisions § Scenario constraint extraction (c)). Flag
  for Automator at spec design; a reference-file write proposal to add Sketcher-
  modal-specific selectors to `grok-browser/references/widgets/chem-sketcher.md`
  is a Phase 2 deliverable.
- **Alternate "make these structures highlighted" entry point in #2448.** Original
  text offers an OR-choice (scaffold tree OR structure filter on `isosmiles` column).
  `sketcher-ui.md` preserves the OR-choice verbatim. Flag for QA pair review on
  first manual run — if one path is canonically preferred for exercising the
  stereo-preservation invariant, the scenario can be tightened to a single path.
- **#2448 `coverage_type` smoke-vs-edge.** Olena's directive flagged this as a
  judgment call. The migration selected `smoke` for the manual-layer ui-smoke
  placement per section convention; could be `edge` if the visual-fidelity
  stereo-preservation invariant is later judged to be an edge-case rather than a
  smoke. No code/text consequence on this cycle (manual scenario); flag for
  re-classification at section-complete Critic F.
- **`pyramid_layer: bug-focused` chain annotation inert after migration.** The
  chain YAML annotates `sketcher.md` as `pyramid_layer: bug-focused` based on the
  pre-rev-2 source body that included #1608 + #2448 sub-blocks. After the chain
  rev 2 directive transformations (drop #1608 + split #2448), the migrated
  `sketcher.md` no longer carries bug-focused content; the chain annotation is
  inert. No action required on this migration (the migrated frontmatter uses
  `coverage_type: regression` directly; chain-layer `pyramid_layer` annotation is
  metadata-not-frontmatter). Flag for chain author at rev 3 if the chain wants to
  re-annotate sketcher.md as `pyramid_layer: integration` post-split.
- **Cross-cutting bug citations not propagated as `related_bugs`.** The three
  chain `bug_focused_candidates[]` entries that span `sketcher.md` (GROK-12758
  Step 2, GROK-14028 Step 6, GROK-17964 Step 6) reference the scenario as a
  pointer surface, not as a covered bug. Per migration-prompt's
  cross-cutting-citation discipline these are not added to `related_bugs` in
  frontmatter (which is reserved for bugs whose reproduction is performed by the
  scenario). Flag noted; no action this cycle.
