# Migration Report — project-url.md

## Step mapping

The original is 4 numbered steps (Markdown source uses `1.` for all
items, auto-numbered as 1-4 on render) + an "Expected" line + trailing
JSON metadata `{ "order": 4 }`. All 4 steps are preserved; the
"Expected" line becomes an explicit verification step.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Go to Browse > Dashboards" | Scenarios > "URL deep-link ..." step 1 | preserved |
| 2. "Click the project (from the previous step): original / copied with the link / copied with clone / saved with personal view customizations" | Scenarios > step 2 (with 4 variant sub-bullets) + parameterisation note | preserved as a parameterised over 4 variants. Source-text-correction context: "with layout" → "saved with personal view customizations" was already applied to source on disk per `mig-2026-04-29-source-text-correction`; migrated body preserves the corrected 4-variant list verbatim. |
| 3. "Go to the Context Panel > Links and copy URL" | Scenarios > step 3 | preserved |
| 4. "Open a new tab in browser and insert the URL" | Scenarios > step 4 | preserved (split: open tab + paste/navigate) |
| "Expected: the corresponding project should be opened" | Scenarios > step 5 (verification) | preserved as explicit verification step. D-STEP-02 axis-clarity: original's bare "Expected" line becomes a numbered verification step listing what specifically is verified (tables, viewers, view layout render; no console errors). |
| Trailing JSON `{ "order": 4 }` | (dropped from body) | metadata-not-step (per orchestrator chain analysis convention; original `order: 4` captured in `scenario-chains/projects.yaml` rev 2 `order_from_files`) |

No original numbered step or expected-result is silently dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/projects.yaml` rev 2
  `output_plan.project-url.md.target_layer = playwright`. The
  scenario exercises new-tab navigation, browser session continuity,
  URL deep-link parsing, and full project re-render — all
  Playwright-driven. `api-contract` cannot exercise the new-tab UI
  flow.
- **Why this `priority`:** chose `regression` per A-STRUCT-MECH-06
  enum (`smoke | regression | edge | perf`). The URL-deep-link reopen
  flow is regression-prone (URL serialisation/deserialisation, query
  param round-trip, view-layout reconstruction across sessions) and
  the 4-variant matrix verifies uniform behavior — NOT a single
  golden-path smoke. Per-cycle Invariant 1 honored: Migrator
  SKILL.md:222 (`p0|p1|p2|p3`) NOT consulted.
- **Why this `strategy`:** `end_to_end_fixtures` per
  `scenario-chains/projects.yaml` rev 2
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
  - No `project-url-spec.ts` exists at any path under `public/`.
    The Automator stage will create the spec fresh in a future
    cycle. This scenario's migration is therefore unconstrained by
    an existing spec's choices.
  - Adjacent specs (`browser-spec.ts`, `deleting-spec.ts`,
    `opening-spec.ts`, `uploading-spec.ts`) all use the
    `loginToDatagrok`, `softStep`, `evalJs`, `closeAll` convention
    — that convention applies here too. Read-only inspection
    confirms `Date.now()` suffix pattern is the section-wide
    naming convention.
- **Helpers reused / candidate:** the migration body identifies a
  candidate Playwright-layer helper:
  `helpers.playwright.projects.buildVariantsComposite(page, originalSource, copyOptions)`
  to coordinate fixture setup across two producers (the original
  source from upload/uploading, and the link/clone/personal-view
  variants from `projects-copy-clone.md`). Flagged for B14 propose-
  only flow; not applied autonomously.
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
    reproduction path differs.
  - Other bugs (GROK-19212, GROK-19103, GROK-19403, GROK-18345,
    GROK-19728) do not intersect this URL-deep-link flow.
  `related_bugs: []` is the correct call.
- **Decision log queried:** yes — `decision-log.yaml` rev 8 read.
  The `mig-2026-04-29-source-text-correction` entry directly
  applies to this scenario: source step 2 was edited to use the
  corrected "personal view customizations" terminology before
  this migration ran. Migrated body honors the correction
  verbatim. The `mig-2026-04-29-fixture-synthesized-inline` entry
  also applies: it confirms that "save personal view
  customizations" is a save-flow flag produced by
  `projects-copy-clone.md`, not an external fixture — so the
  Setup section's fixture prerequisite list correctly attributes
  it to `projects-copy-clone.md` as the producer.
- **Per-cycle override invariants (all 3) status:**
  - Invariant 1 (priority enum source-of-truth): honored —
    `priority: regression` is canonical per A-STRUCT-MECH-06.
  - Invariant 2 (existing -spec.ts / -api.ts READ-ONLY):
    SATISFIED — no project-url spec file exists; trivially
    satisfied. Sibling specs read for convention only and not
    modified.
  - Invariant 3 (atlas-aware sub_features_covered for share):
    correctly does NOT apply — step 3's "Context Panel > Links"
    is a clipboard-copy of a deep-link URL, NOT the right-click
    Share dialog. The URL-deep-link mechanism is covered by
    `projects.url-params.build-share-link` (already in the
    sub_features_covered list); `projects.shell.share-via-
    context-menu` is correctly OMITTED.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — every numbered step + the Expected verification of the
original is preserved; the parameterisation over 4 variants is
preserved.

## Deferred items (NOT opt-outs)

- **Composite fixture build.** The 4 variants come from 2 different
  producers (`upload-project.md` / `uploading.md` for "original";
  `projects-copy-clone.md` for the 3 copy-mode variants). The
  Automator's `beforeAll` block needs to either (a) run the
  upstream scenarios sequentially in the chain, OR (b) build the
  fixture state via direct API replay. Real prerequisite, not
  effort.
- **Browser session continuity across new tabs.** Step 4 opens the
  URL in a "new tab in browser". In Playwright, this can be
  implemented as a new `BrowserContext.newPage()` (preserves
  authentication cookies) OR a new `Browser.newContext()` (fresh
  session, would re-authenticate). The original scenario implies
  the SAME authenticated session continues — the new tab inherits
  cookies. Automator decides between the two implementations at
  spec time.
- **URL clipboard-copy mechanism.** Step 3 says "copy URL" from
  Context Panel > Links. In Playwright, clipboard access requires
  permissions setup; an alternative is to read the URL value
  directly from the DOM (the displayed link's `href` attribute or
  inner-text) without going through the OS clipboard. Automator
  decides at spec time.
- **"Original project" identity disambiguation.** The migrated
  body's Setup section says the "original" can be either the
  `demog` project from `upload-project.md` OR any
  `Test_Case<N>_Sync` from `uploading.md`. Automator decides which
  specific project to use as "original" based on chain-run state
  and fixture availability.

## Edge cases

The original lists no explicit edge cases. Implicit edge cases
derivable from the steps:

- **URL with stale view-layout reference.** If the original
  project's view layout is later modified, the copied URL still
  encodes the OLD layout state. New-tab open should either render
  the encoded state OR fetch the latest server state — behavior
  depends on URL parameter design (snapshot vs. live).
- **URL with session-bound permissions.** If the URL-copy session
  has different permissions than the new-tab session (unlikely
  in same-browser scenario but possible if cookies are invalidated
  mid-test), the new-tab open may fail with a permission error.
- **Variant-specific render differences.** The 4 variants may
  render differently on URL reopen — e.g. "copied-with-link"
  retains the link to the original tables, while "copied-with-
  clone" has its own table copies. The verification step says
  "the corresponding project should be opened" which is satisfied
  if the project loads with its variant-specific content. The
  exact render-correctness assertion per variant is left to
  Automator.
- **URL with embedded URL params (cross-variant).** If a variant
  was saved with URL parameters (e.g. a filter or row-selection),
  those params should be encoded in the deep-link URL and applied
  on reopen. Implicit; verified by the project rendering with the
  same state on URL load.

(none additional)

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
  `grok-browser/references/` for Context Panel structure and the
  Links section selector at spec time.
- **URL format & query parameters.** The original is silent on
  the URL format. The atlas sub_feature `projects.url-params.build-
  share-link` describes URL construction from project nqName +
  param mapping + FuncCall values. Automator should NOT assume a
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
- **Source-text-correction recurrence.** The "with layout" → "save
  personal view customizations" correction was applied 2026-04-29
  per `mig-2026-04-29-source-text-correction`. If the source is
  edited again to revert to "with layout" wording (which lacks a
  producer), this scenario's migration would need to be re-run
  with the corrected wording. Flag for retro: should source-text
  correction history be auto-detected by the Migrator and surfaced
  as a Direct-answer if the source diverges from the corrected
  state?
