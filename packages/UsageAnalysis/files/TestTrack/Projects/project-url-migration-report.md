# Migration Report — project-url.md

## Step mapping

The original is 4 numbered steps (Markdown source uses `1.` for all
items, auto-numbered as 1-4 on render) + an "Expected" line + trailing
JSON metadata `{ "order": 4 }`. All 4 steps are preserved; the
"Expected" line becomes an explicit verification step.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Go to Browse > Dashboards" | Scenarios > "URL deep-link reopen for each variant" step 1 | preserved |
| 2. "Click the project (from the previous step): original / copied with the link / copied with clone / saved with personal view customizations" | Scenarios > step 2 (with 4 variant sub-bullets) + parameterisation note | preserved as parameterised over 4 variants. Source-text-correction context: "with layout" → "saved with personal view customizations" was already applied to source on disk per `mig-2026-04-29-source-text-correction`; migrated body preserves the corrected 4-variant list verbatim. |
| 3. "Go to the Context Panel > Links and copy URL" | Scenarios > step 3 | preserved |
| 4. "Open a new tab in browser and insert the URL" | Scenarios > step 4 | preserved (split: open tab + paste/navigate) |
| "Expected: the corresponding project should be opened" | Scenarios > step 5 (verification) | preserved as verification. D-STEP-02 axis-clarity: original's bare "Expected" line becomes a numbered verification step listing what is verified (tables, viewers, view layout render; no console errors). |
| Trailing JSON `{ "order": 4 }` | (dropped from body) | metadata-not-step (per orchestrator chain analysis convention; original `order: 4` captured in `scenario-chains/projects.yaml` rev 3 `order_from_files`) |

No original numbered step or expected-result is silently dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 3
  `output_plan.project-url.md.target_layer = playwright`. The
  scenario exercises new-tab navigation, browser session continuity,
  URL deep-link parsing, and full project re-render — all
  Playwright-driven. `api-contract` cannot exercise the new-tab UI
  flow.
- **Why this `coverage_type`:** chose `regression` per
  A-STRUCT-MECH-06 enum (`smoke | regression | edge | perf`). The
  URL-deep-link reopen flow is regression-prone (URL serialisation /
  deserialisation, query-param round-trip, view-layout reconstruction
  across sessions) and the 4-variant matrix verifies uniform behavior
  — NOT a single golden-path smoke. Per-cycle Invariant 1 honored:
  Migrator SKILL.md severity enum (`p0|p1|p2|p3`) NOT consulted.
- **Why this `pyramid_layer`:** `integration` per chain rev 3
  classification (Rule 4 Variant C — source-agnostic operation).
  The Context-Panel-Links URL-copy + new-tab URL-apply path works
  identically across all six atlas source classes (files, query,
  script, spaces, db_table, derived); the 4-variant enumeration
  exercises uniformity across variant modes (link / clone /
  personal-customizations / original) but NOT across source classes
  per se.
- **Source-agnostic representative source (Rule 4 Variant C):**
  **file-share (`demog`)**. Selected because
  (a) `upload-project.md` is the cheapest, smallest upstream
  producer in the chain — its `demog-project-with-viewers` fixture
  (js-api-replay reusable_in: project-url.md) is available with
  minimal setup cost; (b) the three copy-mode variants
  (`<original>-link`, `<original>-clone`,
  `<original>-personal-view-customizations`) are derived by
  `projects-copy-clone.md` from the same `demog` baseline, so all
  4 variants share file-share as the underlying source class;
  (c) the existing `project-url-spec.ts` (Wave 1a B70 follow-up)
  already uses `System:DemoFiles/demog.csv` as the test fixture,
  aligning .md and spec on the same representative source. Other
  source classes are covered by atlas-driven
  `proactive_lifecycle_specs` (one spec per source_class × dep_op
  cell) per chain rev 3, NOT by this scenario.
- **Why this `strategy`:** `end_to_end_fixtures` per
  `scenario-chains/projects.yaml` rev 3
  `output_plan.project-url.md.strategy = end_to_end_fixtures`. The
  scenario has cross-fixture dependencies (upload + uploading +
  copy-clone-customizations-variants) — the Automator's `beforeAll`
  block coordinates a composite fixture build before exercising the
  4-variant URL-reopen test body. `chained_tests` would not capture
  the multi-producer fixture nature; `data_driven` would imply matrix
  axes which this scenario does not have (it has a fixed 4-variant
  enumeration, not a multi-axis cross-product).
- **Sibling tests consulted (READ-ONLY per Invariant 2 of the
  per-cycle override):**
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/project-url-spec.ts`
    (Wave 1a B70 follow-up — exists; SCOPE_REDUCTION to single
    `demog` variant via `page.goto({BASE}/p/{nqName})`). Read-only
    inspection; the rev-3 .md migration aligns with the spec's
    existing `demog`-as-representative-source choice.
  - Adjacent specs (`opening-spec.ts`, `uploading-spec.ts`,
    `projects-copy-clone-spec.ts`) all use the `loginToDatagrok`,
    `softStep`, `evalJs`, `closeAll` convention — that convention
    applies here too. Read-only inspection confirms `Date.now()`
    suffix pattern is the section-wide naming convention.
- **Helpers reused / candidate:** the migration body identifies a
  candidate Playwright-layer helper:
  `helpers.playwright.projects.buildVariantsComposite(page, originalSource, copyOptions)`
  to coordinate fixture setup across two producers (the original
  source from upload/uploading, and the link/clone/personal-view
  variants from `projects-copy-clone.md`). Flagged for B14 propose-
  only flow; not applied autonomously. Already noted in
  `decision-log.yaml :: mig-2026-04-30-project-url-migration ::
  candidates_flagged_for_followup`.
- **Bug library consulted:** yes — `bug-library/projects.yaml`
  rev 2. No bug intersects this scenario's flow:
  - GROK-19750 (save-copy with Link drops viewers) — concerns the
    save-copy operation, not the subsequent URL-deep-link reopen.
    The "copied-with-link" variant in this scenario IS a
    save-with-link product, so a regression of GROK-19750 would
    manifest as the "copied-with-link" variant having lost viewers
    when reopened via URL — but the visible failure mode is
    different from this scenario's verification step (which checks
    that the URL opens A project, not that the project is viewer-
    complete). Not listed in `related_bugs` because the
    reproduction path differs; coverage is incidental at best.
  - Other bugs (GROK-19212, GROK-19103, GROK-19403, GROK-18345,
    GROK-19728, github-3550) do not intersect this URL-deep-link
    flow.
  `related_bugs: []` is the correct call.
- **Cross-cutting bug citation (chain rev 3
  `bug_focused_candidates`):** chain rev 3 lists no
  `bug_focused_candidates[]` entry whose `spans[]` references
  `project-url.md` — this scenario is not a cross-cutting span for
  any curated bug. Citation N/A. (F-BUG-COVERAGE-01 at
  section-complete will verify chain YAML coverage authoritatively.)
- **Decision log queried:** yes — `decision-log.yaml` (rev current)
  read; entry `mig-2026-04-30-project-url-migration` honored as the
  rev-2 baseline for this re-migration. The
  `mig-2026-04-29-source-text-correction` entry directly applies:
  source step 2 was edited to use the corrected "personal view
  customizations" terminology before this migration ran. Migrated
  body honors the correction verbatim. The
  `mig-2026-04-29-fixture-synthesized-inline` entry also applies:
  it confirms that "save personal view customizations" is a
  save-flow flag produced by `projects-copy-clone.md`, not an
  external fixture — so the Setup section's fixture prerequisite
  list correctly attributes it to `projects-copy-clone.md` as the
  producer. No `failed_attempts` for `feature: projects` overlap
  with this scenario's chosen approach.
- **UI delegation under SCOPE_REDUCTION:** N/A — no SCOPE_REDUCTION
  is proposed at the .md level for this migration. The existing
  `project-url-spec.ts` Wave 1a B70 follow-up applies a
  spec-level SCOPE_REDUCTION (1 variant instead of 4) but that
  SR is owned by the spec, not by this migration; it is documented
  in the spec's header comment. The .md preserves all 4 variants
  per D-STRUCT-02.
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `coverage_type: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — `project-url-spec.ts` exists, was read for
    convention, and was NOT modified by this migration.
  - Invariant 3 (atlas-aware sub_features_covered for share):
    correctly does NOT apply — step 3's "Context Panel > Links"
    is a clipboard-copy of a deep-link URL, NOT the right-click
    Share dialog. The URL-deep-link mechanism is covered by
    `projects.url-params.build-share-link` (already in the
    sub_features_covered list); `projects.shell.share-via-
    context-menu` is correctly OMITTED.
- **Rev 3 schema fields populated:** `pyramid_layer: integration`,
  `ui_coverage_responsibility: [context-panel-links-url-copy,
  new-tab-open-url]`, `ui_coverage_delegated_to: null` — sourced
  directly from chain rev 3 `dependency_graph.project-url.md`.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every numbered step + the Expected verification of the
original is preserved at the `.md` level; the parameterisation over
4 variants is preserved per D-STRUCT-02.

## Deferred items (NOT opt-outs)

- **Composite fixture build.** The 4 variants come from 2 different
  producers (`upload-project.md` / `uploading.md` for "original";
  `projects-copy-clone.md` for the 3 copy-mode variants). The
  Automator's `beforeAll` block needs to either (a) run the
  upstream scenarios sequentially in the chain, OR (b) build the
  fixture state via direct API replay. Real prerequisite: chain
  ordering enforced by orchestrator + composite fixture helper not
  yet registered (`helpers.playwright.projects.buildVariantsComposite`
  flagged as candidate).
- **Browser session continuity across new tabs.** Step 4 opens the
  URL in a "new tab in browser". In Playwright, this can be
  implemented as a new `BrowserContext.newPage()` (preserves
  authentication cookies) OR a new `Browser.newContext()` (fresh
  session, would re-authenticate). The original scenario implies
  the SAME authenticated session continues — the new tab inherits
  cookies. Real prerequisite: Playwright session-handling pattern
  decision (BrowserContext.newPage vs newContext) not yet codified
  in `helpers-registry.yaml :: playwright_layer.session`.
- **URL clipboard-copy mechanism.** Step 3 says "copy URL" from
  Context Panel > Links. In Playwright, clipboard access requires
  permissions setup; an alternative is to read the URL value
  directly from the DOM (the displayed link's `href` attribute or
  inner-text) without going through the OS clipboard. Real
  prerequisite: clipboard-permissions helper or DOM-read fallback
  helper not yet registered.
- **"Original project" identity disambiguation.** The migrated
  body's Setup section says the "original" can be either the
  `demog` project from `upload-project.md` (Variant C
  representative source) OR any `Test_Case<N>_Sync` from
  `uploading.md`. Real prerequisite: chain-run state introspection
  helper or runtime fixture-availability probe not yet registered;
  Automator decides which specific project to use as "original"
  based on chain-run state and fixture availability.

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable from the steps:

- **URL with stale view-layout reference.** If the original
  project's view layout is later modified, the copied URL still
  encodes the OLD layout state. New-tab open should either render
  the encoded state OR fetch the latest server state — behavior
  depends on URL parameter design (snapshot vs. live). Flagged
  for atlas curator to mark `manual_only` or capture as a separate
  scenario; not in scope here.
- **URL with session-bound permissions.** If the URL-copy session
  has different permissions than the new-tab session (unlikely in
  same-browser scenario but possible if cookies are invalidated
  mid-test), the new-tab open may fail with a permission error.
  Flagged for atlas curator; cross-session permission boundary is
  a separate scenario shape.
- **Variant-specific render differences.** The 4 variants may
  render differently on URL reopen — e.g. "copied-with-link"
  retains the link to the original tables, while "copied-with-
  clone" has its own table copies. The verification step says
  "the corresponding project should be opened" which is satisfied
  if the project loads with its variant-specific content. The
  exact render-correctness assertion per variant is left to
  Automator (preserved as scenario step at the parameterised
  verification level).
- **URL with embedded URL params (cross-variant).** If a variant
  was saved with URL parameters (e.g. a filter or row-selection),
  those params should be encoded in the deep-link URL and applied
  on reopen. Implicit; verified by the project rendering with the
  same state on URL load (preserved in the Step 5 verification).
- **Source-class variation across the URL-deep-link path.** The
  Variant C representative-source claim (file-share `demog`) is
  itself a source-agnostic invariant that the URL flow holds for
  query / script / spaces / db_table / derived sources too. NOT
  exercised here (out of scope per Variant C); covered by atlas-
  driven `proactive_lifecycle_specs` (one spec per source_class
  × dep_op cell) per chain rev 3.

## Unresolved ambiguities

- **"Click the project" — single-click vs. double-click.** The
  original step 2 says "Click the project" without specifying
  single- or double-click. Single-click selects-for-Context-Panel-
  inspection; double-click opens. Step 3 ("Go to Context Panel >
  Links") implies the project must be SELECTED but NOT yet opened
  (Context Panel reflects the selection of the Dashboards card).
  The migrated body inherits this ambiguity — Automator needs to
  use single-click to keep the project in selected state for
  Context Panel inspection in step 3.
- **"Context Panel > Links" location.** The original assumes the
  Context Panel has a "Links" section/tab containing the project's
  deep-link URL. The exact UI control name and DOM hook are not
  specified. Automator must cross-reference
  `grok-browser/references/projects.md:81-101` for Context Panel
  structure and the Links section selector at spec time. The
  existing spec sidesteps the DOM-read by computing the URL via
  `nqName` from JS API; that approach satisfies the URL-apply
  invariant but bypasses the Context-Panel-Links UI surface that
  rev 3 `ui_coverage_responsibility` explicitly attributes to
  this scenario. Flag for retro: spec-level SR may need to be
  re-examined when the Context-Panel-Links DOM read pattern lands.
- **URL format & query parameters.** The original is silent on
  the URL format. The atlas sub_feature `projects.url-params.
  build-share-link` describes URL construction from project nqName
  + param mapping + FuncCall values. Automator should NOT assume a
  specific URL format; instead read whatever URL the Context
  Panel > Links section displays. If the URL has expiry or
  session-bound tokens, the new-tab test may fail intermittently.
- **"New tab" vs. "incognito tab".** The original says "Open a
  new tab in browser" — this implies the same browser session
  (cookies preserved). An incognito tab would force re-auth.
  Migrated body uses "new tab" in the same browser session per
  the original's likely intent. If Automator uses incognito, the
  test verifies a different scenario (cross-session URL reopen)
  which is OUT of scope here.
- **Order vs. dependency contradiction.** `project-url.md`
  trailing metadata is `order: 4` but the scenario depends on
  `projects-copy-clone.md` which has `order: 5`. The chain rev 3
  `unresolved_ambiguities` records this; dependency graph
  follows the named-variant evidence, treating `order` as
  advisory only. No additional action needed at migration time.
- **Source-text-correction recurrence.** The "with layout" → "save
  personal view customizations" correction was applied 2026-04-29
  per `mig-2026-04-29-source-text-correction`. If the source is
  edited again to revert to "with layout" wording (which lacks a
  producer), this scenario's migration would need to be re-run
  with the corrected wording. Flag for retro: should source-text
  correction history be auto-detected by the Migrator and surfaced
  as a Direct-answer if the source diverges from the corrected
  state?
