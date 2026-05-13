# Migration Report — upload-project.md

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open demog.csv" | Setup step 1 | preserved (dataset path resolved to `System:DemoFiles/demog.csv` from the original's trailing `datasets` JSON metadata) |
| 2. "Add scatter plot, bar chart and line chart viewers" | Setup steps 2–4 | preserved (split for clarity — one viewer per step) |
| 3. "Save a project. In the opened Save Project dialog keep the data sync on" | Scenarios > "Save with Data Sync on" steps 1–2 | preserved |
| 4. "Click OK - wait for project tot be saved" | Scenarios > "Save with Data Sync on" step 3 | preserved (typo "tot" → "to" silently corrected; original wait-for-save now expressed as an explicit verification of project existence on the server) |
| 5. "The Share demog dialog opens - click Cancel" | Scenarios > "Save with Data Sync on" step 4 | preserved as verification (the auto-open of the Share dialog AND its dismissal are both verified) |
| 6. "Close All" | Scenarios > "Save with Data Sync on" step 5 | preserved as cleanup |

## Decisions

- **Why this `target_layer`:** chose `playwright` because the existing
  sibling specs in `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts`
  (Cases 1, 3, 4, 9 — all "Save project from ..." flows) are Playwright,
  and this scenario is a smaller/simpler variant of the same workflow.
  Aligning the layer keeps the section's house style consistent. The
  scenario also exercises an OS-style modal handoff (Save dialog →
  Share dialog → Cancel) that is more naturally driven by Playwright
  than by a JS API package test.
- **Why this `strategy`:** `simple` — single standalone scenario, no
  matrix axes, no chained tests in this file. (Step 1 chain analysis,
  if run on the Projects section, will likely identify this scenario
  as `must_run_first` and producer of a `demog-project-with-viewers`
  fixture for downstream scenarios; that's a Step 1 concern.)
- **Sibling tests consulted:**
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts`
    — confirms Playwright is the sectional convention; Cases 1/3/4/9
    show the save-with-sync UI flow shape.
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading.md`
    — adjacent manual scenario covering the broader "save from N
    sources" matrix; this `upload-project.md` is the smaller golden
    path for it.
- **Helpers reused:**
  - None directly. The `uploadProject(projectName, tableInfo, view, df)`
    helper at `public/packages/UITests/src/gui/gui-utils.ts:100`
    matches semantically but is registered as `grok_test_layer` (a/b
    style only). It is NOT Playwright-compatible, so it cannot be
    reused at the chosen target layer.
  - **Candidate helper (not invented; flagged for registry):**
    `helpers.playwright.projects.uploadDemogWithViewers(page)` — a
    Playwright-layer counterpart to `uploadProject` that drives the
    Save Project dialog and dismisses the Share dialog. Propose to
    Andrew for addition to `helpers-registry.yaml :: playwright_layer`
    if Step 2 of the orchestrator chain produces ≥2 scenarios that
    would reuse it.
- **Bug library consulted:** yes — `bug-library/projects.yaml` is
  present (8 curated bugs as of 2026-04-28). Three intersect this
  scenario's coverage (`projects.upload`, `projects.api.save`):
  - `GROK-19750` — Save Copy with table link mode drops original's
    viewers. Reproduction is for save-copy, NOT this golden-path
    save; belongs in `projects-copy-clone.md` migration.
  - `GROK-19212` — Project fails to open after rename + save with
    data sync. Reproduction is rename-then-save, NOT this scenario;
    belongs in a dedicated rename+sync scenario (atlas curation).
  - `GROK-19103` — Join result silently saved as separate project
    that fails to open. Reproduction is join-in-project, NOT this
    scenario; belongs in `complex.md` or a dedicated join+save
    scenario.
  These three are listed in the migrated frontmatter's `related_bugs`
  for traceability but are intentionally NOT exercised by this smoke
  — they have their own scenarios.
- **Decision log queried:** skipped — `decision-log.yaml` does not
  exist (Block 10 not yet run; expected per plan 01).

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every step of the original was preserved at the same
target layer; no SCOPE_REDUCTION was proposed for this migration.

## Deferred items (NOT opt-outs)

(none) — every step of the original is in the migrated body.

## Edge cases

The original scenario lists no explicit edge cases and contains no
implicit ones beyond the verifications already preserved as steps
(project saved successfully → step 3; Share dialog auto-opens →
step 4; Share dialog dismissable via Cancel → step 4).

(none additional)

Note: feature-wide edge cases for `projects.api.save` documented in
`bug-library/projects.yaml` (GROK-19750, GROK-19212, GROK-19103,
github-3550, etc.) are **not** in the original of this scenario;
they belong in separate scenarios specifically targeting
save-copy-link, rename+sync, join+save, and entity-rename
invalidation flows. They are listed in this report's Decisions
section under "Bug library consulted" for traceability.

## Unresolved ambiguities

- **Step 3 — "keep the data sync on":** the original's phrasing
  implies the Save Project dialog defaults `Data Sync` to ON and the
  user simply leaves it; it does not explicitly state whether the
  user must confirm the toggle's current state or actively verify it.
  Migrated as "leave **Data Sync** toggled ON" assuming default-on
  (consistent with `uploading-spec.ts` Cases 1/3/4/9 which all save
  with sync enabled by default). Flag for retro if Gate B reveals
  the default has shifted.
- **Step 4 — "wait for project to be saved":** the original does not
  state HOW to assert the save completed. Migrated as "verify the
  project exists on the server" with two suggested signals (Browse >
  Dashboards entry, or the underlying `POST /projects` response).
  Automator should pick one; if both are flaky in practice, this is
  a candidate for a `helpers.playwright.projects.waitForSaved(page)`
  registry entry.
- **Step 5 — "The Share demog dialog opens":** assumes the dialog
  title contains the saved project's name (`demog`). If project
  naming changes (e.g. timestamp suffix to avoid collisions in CI),
  the title-based assertion may need to switch to a structural
  selector. Flag for Gate B if observed flake on dialog identification.
- **Original's typo "tot be saved" → "to be saved":** silently
  corrected in the migrated body. No semantic change.
