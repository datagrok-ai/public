# Migration Report — uploading.md

Re-migrated 2026-05-04 against `scenario-chains/projects.yaml` rev 3
(supersedes the 2026-04-30 migration). Chain entry classifies this
scenario as `classification: matrix`, `pyramid_layer: source-matrix`,
with `ui_coverage_delegated_to: upload-project.md`. All 16 logical
paths (8 cases × 2 sync states) are preserved at the .md level per
D-STRUCT-02 — honoring the prior `mig-2026-04-30-uploading-migration`
SCOPE_REDUCTION decision (4-of-9 spec coverage was an Automator-stage
decision, not a Migrator opt-out). 18→16 path reduction is the
Phase B cleanup 2026-05-05 (Case 7 split out — see addendum below).

## Step mapping

The original was a 9-case matrix with each case running twice (Data
Sync ON, Data Sync OFF) — 18 logical paths total. Phase B cleanup
2026-05-05 split Case 7 (Get Top 100 / Get All + DB table
double-click + Join Tables) out to
`../Queries/get-all-get-top-100.md`, leaving an 8-case matrix
(numbered 1-6, 8, 9 — Case 7 missing) × 2 = 16 logical paths.
Step numbering in the original is per-case (each Case starts at
step 1). The migrated body preserves all 8 remaining cases under a
single `## Scenarios` H2 with each Case as an `### Case N` H3, and
preserves both Sync ON and Sync OFF variants within each Case body.

| Original location                                                | Migrated location                                              | Decision                                                                                                                                                                                                          |
|------------------------------------------------------------------|----------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `### Matrix: Table source combinations` (header + 9-row table)   | `## Scenarios > ### Matrix: Table source combinations` (table) | preserved (matrix table moved under the `## Scenarios` H2; row content unchanged)                                                                                                                                 |
| Spaces prelude (Setup section, JS API)                           | `## Setup > Spaces prelude` (JS API block)                     | preserved (sa-2026-05-03-spaces-inline-prelude-pattern citation added)                                                                                                                                            |
| Spaces postlude (Setup section, JS API)                          | `## Setup > Spaces postlude` (JS API block)                    | preserved                                                                                                                                                                                                         |
| Case 1, original steps 1–6                                       | Case 1, migrated steps 1–6                                     | preserved (Browse > Files navigation, dataset double-click, Link Tables, verify linking)                                                                                                                          |
| Case 1, original step 7 ("Save with Data Sync ON")               | Case 1, migrated step 7                                        | preserved                                                                                                                                                                                                         |
| Case 1, original steps 8–9                                       | Case 1, migrated steps 8–9                                     | preserved (close, reopen, verify on reopen)                                                                                                                                                                       |
| Case 1, original step 10 ("Save with Data Sync OFF")             | Case 1, migrated step 10                                       | preserved                                                                                                                                                                                                         |
| Case 1, original steps 11–12                                     | Case 1, migrated steps 11–12                                   | preserved                                                                                                                                                                                                         |
| Case 2, original steps 1–12                                      | Case 2, migrated steps 1–12                                    | preserved                                                                                                                                                                                                         |
| Case 3, original steps 1–12                                      | Case 3, migrated steps 1–12                                    | preserved                                                                                                                                                                                                         |
| Case 4, original steps 1–12                                      | Case 4, migrated steps 1–12                                    | preserved                                                                                                                                                                                                         |
| Case 5, original steps 1–10                                      | Case 5, migrated steps 1–10                                    | preserved                                                                                                                                                                                                         |
| Case 6, original steps 1–10                                      | Case 6, migrated steps 1–10                                    | preserved                                                                                                                                                                                                         |
| Case 7, original steps 1–14                                      | extracted to `../Queries/get-all-get-top-100.md`               | Phase B cleanup 2026-05-05 — Get Top 100 / Get All + Join Tables flow moved to its dedicated home; case-numbering left at 1-6, 8, 9 (Case 7 missing) so existing spec mappings stay valid                         |
| Case 7 — implicit GROK-19212 / GROK-19103 / GROK-18345 framing   | replaced by `>` summary block in uploading.md                  | Phase B cleanup 2026-05-05 — uploading.md now carries a brief callout pointing at the extracted Case 7 + redistributed bug anchors (GROK-19103 → Cases 8-9; GROK-18345 → Cases 4-6; GROK-19212 → complex.md only) |
| Case 8, original steps 1–12                                      | Case 8, migrated steps 1–12                                    | preserved                                                                                                                                                                                                         |
| Case 9, original steps 1–12                                      | Case 9, migrated steps 1–12                                    | preserved                                                                                                                                                                                                         |
| Notes section (helpers, ordering, env, source-class substitutions) | `## Notes` section                                            | preserved (re-organised under explicit headings; decision-log citations added for the 3 Session A substitutions and the 1 Spaces inline prelude pattern)                                                          |

No original numbered step is silently dropped; every step has an
explicit migrated counterpart (Case 7's 14 steps relocated to
`../Queries/get-all-get-top-100.md` per Phase B cleanup 2026-05-05).
Numbering hygiene matches the prior 2026-04-30 migration; this
re-migration adds the chain rev 3 metadata (`pyramid_layer`,
`ui_coverage_responsibility`, `ui_coverage_delegated_to`).

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 3
  `output_plan.uploading.md.target_layer = playwright`, consistent
  with the existing sibling spec `uploading-spec.ts` (already
  Playwright). Save → reload → verify across navigation requires
  real persistence, ruling out `api-contract`.
- **Why this `pyramid_layer: source-matrix`:** chain rev 3
  classifies this scenario via Rule 2 (matrix discriminator). Matrix
  axes enumerate over data sources (Files / Query / Spaces / DB /
  Pivot / Aggregate), not button states; the `source-matrix` label
  disambiguates from `ui-smoke` (owned by upload-project.md per
  chain `ui_coverage_plan.smoke_scenario`) and `bug-focused` (no
  single dominant bug repro).
- **Why `ui_coverage_delegated_to: upload-project.md`:** Save
  Project dialog and Data Sync toggle smoke is owned by
  upload-project.md per chain rev 3 `ui_coverage_plan.smoke_covers`.
  This scenario owns the residual source-specific dialog flows
  (Link Tables, Join Tables, Pivot Table > Add to workspace,
  Aggregate Rows > Add to workspace) — those 5 flows appear in
  `ui_coverage_responsibility` but the Save-Project / Data-Sync
  smoke surface delegates to upload-project.md. Get Top 100 / Get
  All flows are owned by `../Queries/get-all-get-top-100.md` (split
  out 2026-05-05 — see Phase B cleanup addendum below). UI coverage citation per migration-prompt.md
  "UI delegation in SCOPE_REDUCTION proposals" §: **UI coverage for
  Save Project dialog + Data Sync toggle: upload-project.md
  (declared in chain `ui_coverage_plan.smoke_covers`).**
- **Why `coverage_type: regression`:** matrix coverage of the
  upload+save+reopen lifecycle across 8 source combinations and 2
  sync states is regression-shaped (NOT a single golden-path smoke;
  golden-path lives in upload-project.md).
- **Why this `strategy` (consumed by chain, not in frontmatter):**
  `data_driven` per chain rev 3 `output_plan.uploading.md.strategy`.
  The 8×2 matrix is the canonical data-driven shape — same step
  skeleton parameterised over 8 table-source combinations and 2
  sync states. (Was 9×2=18 before Phase B cleanup 2026-05-05; Case 7
  extracted to `../Queries/get-all-get-top-100.md`.)
- **D-STRUCT-02 compliance — all 16 paths preserved:** Per
  migration-prompt.md "matrix scenarios" §, the 8×2=16 paths are
  preserved at the .md level. The prior decision-log
  `mig-2026-04-30-uploading-migration` (4-of-9 spec coverage) is an
  **Automator-stage** decision, not a Migrator-stage opt-out — the
  migrated .md remains the full source of truth. Re-litigation
  against rev 3: not warranted; chain rev 3 still classifies as
  matrix and prior reduction rationale (env dependencies on Spaces /
  Postgres queries; spec-time helper coverage) holds. Honoring prior
  decision per migration-prompt.md "Decision log query (D10)" §.
- **Sibling tests consulted (per migration-prompt.md "Sibling test
  consultation" §):**
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts` —
    Playwright; covers Cases 1, 3, 4, 9 (Sync ON only); uses local
    `saveProject(page, name)` helper, `evalJs(page, ...)` for
    JS-API replay of table-open, `softStep(...)` for soft-error
    accumulation. Pattern reused at .md level (prose Save Project
    flow; helper invention deferred to Automator).
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/upload-project.md` —
    smoke owner per chain `ui_coverage_plan.smoke_scenario`; this
    scenario delegates Save-Project / Data-Sync UI coverage to it.
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-run.md` —
    2026-03-09 run summary (8/14 step-equivalents passed; Cases
    4–8 skipped due to missing Spaces — closed by
    sa-2026-05-03-spaces-inline-prelude-pattern).
- **Helpers reused (per migration-prompt.md "Helpers discipline"
  §):**
  - At Playwright layer: `saveProject(page, name)` — local helper
    in `uploading-spec.ts` (NOT in helpers-registry.yaml). Migrator
    uses prose ("**File** > **Save Project**") at the .md level;
    helper invocation is Automator's call.
  - At grok_test_layer: `uploadProject(projectName, tableInfo,
    view, df)` at `helpers-registry.yaml :: grok_test_layer` — NOT
    Playwright-compatible; not the chosen target layer for this
    scenario. Same flag as upload-project.md.
  - **Candidate helper (not invented; flagged for registry):**
    `helpers.playwright.projects.saveAndReopen(page, name,
    syncOn)` — Playwright-layer counterpart that drives the Save
    Project dialog with Data Sync toggle, saves under given name,
    closes, and reopens. Would be reused by all 16 paths in this
    scenario plus residual paths in scenarios 3, 5, 6. **Propose
    add to helpers-registry.yaml externally** if Automator
    materialises ≥4 paths reusing it.
- **Bug library consulted (per migration-prompt.md "Bug library
  query" §):** yes — `bug-library/projects.yaml` rev 2 read. Three
  bugs cited via cross-cutting candidate spans in chain rev 3
  `bug_focused_candidates`:
  - `GROK-19103` (Join result silently saved as separate project
    that fails to open). Chain rev 3 spans (post Phase B cleanup):
    `uploading.md:Cases 8 & 9` (Pivot/Aggregate Add to workspace) +
    `upload-project.md:Step 1` + `projects-copy-clone.md:Step 4` +
    `complex.md:Step 1`. Cross-cutting invariant: derivation lands
    in active project, not as stray new project. **Bug GROK-19103
    has cross-cutting candidate spec — see chain
    `bug_focused_candidates` entry, proposed_spec:
    projects-grok-19103-spec.ts.**
  - `GROK-19212` (Projects fail to open with "Could not resolve
    table" after referenced table is renamed). Chain rev 3 spans
    pre-Phase-B: `uploading.md:Case 7 Step 7` + `complex.md:Step 7`.
    Post Phase B cleanup 2026-05-05: Case 7 split out to
    `../Queries/get-all-get-top-100.md`; sole anchor for this
    scenario span is now `complex.md:Step 7`. Cross-cutting
    invariant unchanged. **Bug GROK-19212 has cross-cutting candidate
    spec — see chain `bug_focused_candidates` entry, proposed_spec:
    projects-grok-19212-spec.ts.**
  - `GROK-18345` (Recipient cannot open shared project that uses a
    Spaces dataset saved with data sync). Chain rev 3 spans (post
    Phase B cleanup): `uploading.md:Cases 4-6` (Spaces sources) +
    `complex.md:Step 12`. Cross-cutting invariant: Spaces +
    datasync + share → recipient open succeeds. **Bug GROK-18345
    has cross-cutting candidate spec — see chain
    `bug_focused_candidates` entry, proposed_spec:
    projects-grok-18345-spec.ts.**
  - `related_bugs` frontmatter post Phase B cleanup is
    `[GROK-19103, GROK-18345]` (GROK-19212 removed because Case 7
    extraction left no anchor in this scenario; rev-2 shape
    restored). Earlier rev-3 expansion to
    `[GROK-19103, GROK-19212, GROK-18345]` reverted 2026-05-05.
  - Other bugs do not intersect: GROK-19750 (save-copy, not save
    here), GROK-19403 (share-with-unshared-deps, not in scope),
    GROK-19728 (view-and-use-failure-state, sharing flow), and
    github-3550 (query rename, lives in complex.md only).
- **Decision log queried (per migration-prompt.md "Decision log
  query (D10)" §):** yes — `decision-log.yaml` read. Relevant
  entries honored:
  - `mig-2026-04-30-uploading-migration` — D-STRUCT-02 18-path
    preservation honored verbatim; 4-of-9 reduction remains
    Automator-stage.
  - `sa-2026-05-03-postgres-queries-public-data-substitution` —
    Cases 2/3/6 use `Samples:PostgresCustomers` (substituting
    `PostgresAll`); same-query-twice degradation noted in Notes.
  - `sa-2026-05-03-spgi-file-public-data-substitution` — Cases
    1/5/8 use `spgi-100.csv` from `System:AppData/Chem/tests`
    (substituting `SPGI_v2_infinity.csv` from `System:Demo`).
  - `sa-2026-05-03-spaces-inline-prelude-pattern` — Cases 4–6
    Spaces prelude uses inline JS API
    (`grok.dapi.spaces.createRoot` → `addEntity` → postlude
    `delete`); env dependency on pre-provisioned Space removed.
  - No prior `failed_attempts` entry for `feature: projects` AND
    `scenario: uploading` is on the table; the prior
    `mig-2026-04-30-uploading-migration` is a successful
    migration entry (decision_type: scenario-migration), not a
    failed_attempt.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every numbered step of every case is preserved at the
same target layer; no SCOPE_REDUCTION is proposed for the .md
migration. The fact that `uploading-spec.ts` only covers Cases 1,
3, 4, 9 with Sync ON is a property of the existing spec
(Automator-stage decision per
`mig-2026-04-30-uploading-migration`), not a reduction made by
this migration.

UI delegation citation (per migration-prompt.md "UI delegation in
SCOPE_REDUCTION proposals" §): although no SR is applied, this
scenario delegates UI coverage of the Save Project dialog + Data
Sync toggle smoke to `upload-project.md` per chain rev 3
`ui_coverage_plan.smoke_covers`. The 5 source-specific flows in
this scenario's `ui_coverage_responsibility` (Link Tables, Join
Tables, Pivot Table > Add, Aggregate Rows > Add,
save-project-dialog) are owned here and not delegated. Get Top 100
/ Get All flows split out to `../Queries/get-all-get-top-100.md`
2026-05-05 (Phase B cleanup).

## Deferred items (NOT opt-outs)

- **Case 9 DB-table reference not yet substituted.** Case 9 still
  references `Postgres > Northwind > public > orders` directly via
  plain double-click; the Session A
  `sa-2026-05-03-postgres-queries-public-data-substitution`
  mapping does not cover bare DB-table references.
  **Prerequisite:** decision on whether to map to
  `Samples:PostgresOrders` (Session A scope expansion) or accept
  env-dependent skip on dev. Tracked under decision-log
  `sa-2026-05-03-postgres-queries-public-data-substitution
  related_pending`. (Case 7 — also a DB-table reference
  pre-cleanup — split out 2026-05-05 to
  `../Queries/get-all-get-top-100.md`.)
- **`saveAndReopen` Playwright helper registration.**
  `helpers.playwright.projects.saveAndReopen(page, name, syncOn)`
  candidate is flagged but not invented at migration. **Prerequisite:**
  Automator-stage decision based on whether ≥4 paths reuse it.
  Tracked here for cross-cycle visibility.
- **Concurrent-CI-collision risk on deterministic project names.**
  All 16 saved projects use deterministic names
  (`Test_Case<N>_Sync` / `Test_Case<N>_NoSync`). Concurrent runs
  collide at the server level. **Prerequisite:** Automator
  applies `Date.now()` suffix per existing `uploading-spec.ts`
  pattern (`AutoTest-Upload-Case<N>-<timestamp>`). Tracked at
  Automator-stage; deletion delegated to `deleting.md`
  (must_run_last per chain rev 3).

## Phase B cleanup addendum (2026-05-05)

Case 7 (Get Top 100 / Get All + DB table double-click + Join
Tables) was extracted from this scenario to
`../Queries/get-all-get-top-100.md`. Rationale: that dedicated
scenario is the canonical home for Get Top 100 / Get All UI
coverage; the Projects scenario should not duplicate it as a
data-source step.

Changes synced:
- `uploading.md`: Case 7 block + matrix row + `get-top-100` /
  `get-all` frontmatter entries removed; counts 9→8 cases,
  18→16 paths.
- `feature-atlas/projects.yaml`: Get All example dropped from
  the db_table source-class entry.
- `scenario-chains/projects.yaml`: get-top-100 / get-all
  references attached to uploading.md removed (rev 4 reconciled).
- `coverage-gaps/projects.md`: cross-cuts pointing back to
  uploading.md for `dataops.getTopN` removed; helper-candidate
  itself stays for the dedicated Queries spec.

Case-numbering preserved (1–6, 8, 9 with Case 7 missing) so
existing spec mappings (`uploading-spec.ts` covers 1, 3, 4, 9)
and run artefacts stay valid.

## Edge cases

For each edge case in the original (explicit OR implicit), record
the outcome:

- **Pivot Table preserve-on-reopen** (Case 8) — preserved as
  scenario step (Case 8 step 9 verifies "both the source table
  and the pivot table are loaded").
- **Aggregate Rows preserve-on-reopen** (Case 9) — preserved as
  scenario step (Case 9 step 9 verifies the aggregated table loads
  alongside the source).
- **4 join-types preserve-on-reopen** — extracted with Case 7 to
  `../Queries/get-all-get-top-100.md` (Phase B cleanup 2026-05-05).
  Originally Inner / Outer / Left / Right joins each produced a
  separate result table; Case 7 step 11 verified "all tables
  (source + 4 joined) are loaded".
- **Data Sync toggle persistence** (every case) — implicit;
  verified indirectly by the reopen-and-verify-no-console-errors
  steps. Preserved as scenario step (12-of-16 paths assert this
  via reopen verification post Phase B cleanup).
- **Selection-to-filter linking-direction asymmetry.** Original
  verifies one direction only (selection in Table 1 → filter in
  Table 2). Migrated preserves the original's one-direction
  verification. **Decision: preserved as scenario step.** Reverse
  direction is flagged for atlas curator if a regression surfaces.
- **Project-naming collision risk.** See Deferred items —
  Automator-stage handling.
- **"Add to workspace" semantics for Pivot/Aggregate (Cases 8,
  9).** Migrated preserves "**Verify:** the pivot table opens as a
  separate table in a new tab" (and the aggregate equivalent in
  Case 9). Selector / DOM hook is unspecified at .md level (per
  house style). **Decision: preserved as scenario step.**
  Selector resolution is Automator's call.

Feature-wide edge cases for `projects.api.save`,
`projects.api.files.sync`, and `projects.add-relation` documented
in `bug-library/projects.yaml` rev 2 (GROK-19103, GROK-19212,
GROK-18345) are listed in this scenario's `related_bugs` for
traceability and cited via cross-cutting candidate spans (see
Decisions § "Bug library consulted").

## Unresolved ambiguities

- **"Selection to filter" linking-direction asymmetry.** See Edge
  cases — preserved as one-direction verification matching
  original; flag if reverse-direction regression surfaces.
- **`SPGI_v2_infinity.csv` source ambiguity.** Closed by
  `sa-2026-05-03-spgi-file-public-data-substitution` — substituted
  with `spgi-100.csv` from `System:AppData/Chem/tests`. No longer
  ambiguous.
- **Project-naming collision risk in CI.** See Deferred items —
  Automator-stage.
- **"Add to workspace" semantics** (Cases 8, 9) — selector/DOM
  hook unspecified at .md level. Automator to cross-reference
  `grok-browser/references/viewers/pivot_table.md` at spec time.
- **Original step renumbering across Sync-ON / Sync-OFF blocks.**
  Closed by prior 2026-04-30 migration's normalisation; this
  re-migration preserves the normalised numbering.
- **`order: 1` collision with `upload-project.md`.** Recorded in
  `scenario-chains/projects.yaml` rev 3 `unresolved_ambiguities`
  (priority: low). Migrator does not modify the chain artifact.
