---
feature: powerpack
sub_features_covered:
  - powerpack.search.power
  - powerpack.search.power.dispatch
  - powerpack.search.power.init
  - powerpack.welcome.view
  - powerpack.welcome.suggestion-nav
  - powerpack.search
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Powerpack/power-search-enter.md
migration_date: 2026-05-23
related_bugs:
  - GROK-18656
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - power-search-enter-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T14:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses cleanly as YAML and contains all four
          required fields with the expected types. feature is the
          scalar string powerpack; sub_features_covered is a 6-item
          list of dotted atlas IDs (powerpack.search.power,
          powerpack.search.power.dispatch, powerpack.search.power.init,
          powerpack.welcome.view, powerpack.welcome.suggestion-nav,
          powerpack.search); target_layer is playwright; coverage_type
          is regression.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains the parent H2 ## Scenarios framing plus three
          H3 scenario headings (### Scenario 1: short-query QA repro,
          ### Scenario 2: 8-path exercise, ### Scenario 3:
          suggestion-nav regression guard). The mechanical "starts with
          ## " match is satisfied by both the parent H2 and the H3
          children; the framing H2 sections (## Setup, ## Scenarios,
          ## Notes) bracket the body cleanly.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Each of the three scenario headings is followed by numbered
          steps. Scenario 1 carries 5 numbered steps; Scenario 2 lays
          out a 4-step probe template (focus, type, press Enter,
          verify) which is the unit-of-iteration applied to the 5
          listed query rows; Scenario 3 carries 5 numbered steps.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          No empty scenarios. Every scenario carries at least one
          numbered step (counts 5, 4, 5 respectively).
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer value playwright is a member of the canonical
          enum {playwright, apitest, manual-only}.
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type value regression is a member of the canonical
          test-kind enum {smoke, regression, edge, perf}. No
          severity-axis values (p0..p3) appear anywhere in the
          frontmatter.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type label regression is declared at file-level
          frontmatter and therefore binds all three scenarios. The
          label is drawn from the unified test-kind enum shared with
          atlas interactions[].coverage_type and
          edge_cases[].coverage_type (the 2026-04-30 unified rename) —
          not from the severity axis used by atlas
          critical_paths[].priority or bug-library priority.
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          Repeated preconditions (open Datagrok with PowerPack
          installed, navigate to the Welcome view, open the browser
          DevTools console) are factored into the ## Setup section as
          three numbered preconditions. Each scenario starts at a
          focus-the-input step rather than re-doing the page
          navigation, so the shared preconditions are not duplicated
          across scenario bodies.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer is absent from the YAML frontmatter (the Notes
          section discusses pyramid_layer bug-focused in prose as a
          derivation comment under Rule 3, but does not set it as a
          frontmatter key). The mode rule
          "PASS-by-vacuity when pyramid_layer is absent from
          frontmatter" applies.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Scenarios use concrete query strings (literal QA, a, new,
          user, a Project regex pattern, 1+1, dem) and real atlas
          source paths. Verified against
          public/packages/PowerPack/src/search/power-search.ts (L21
          initSearch, L26 powerSearch) and
          public/packages/PowerPack/src/welcome-view.ts L118
          suggestionMenuKeyNavigation — all functions exist at the
          declared lines. All six entries in sub_features_covered
          resolve to real atlas IDs in
          references/feature-atlas/powerpack.yaml (powerpack.welcome.view
          L87, powerpack.welcome.suggestion-nav L94, powerpack.search
          L310, powerpack.search.power L466, powerpack.search.power.init
          L473, powerpack.search.power.dispatch L480). No
          angle-bracket placeholders, no square-bracket placeholders,
          no generic stand-ins like column1, the column, or
          some-file.csv.
      - check: A-BUG-01
        status: PASS
        evidence: |
          Single-scenario scope satisfied. Atlas known_issues block at
          references/feature-atlas/powerpack.yaml lines 1349-1445 carries
          9 curated bug entries with the schema test_coverage.exists
          false plus empty test_coverage.paths (the equivalent of the
          test_coverage needed predicate). GROK-18656 appears at atlas
          line 1425 with affects_sub_features exactly matching this
          scenario's six sub_features, and is addressed via the
          related_bugs frontmatter entry plus explicit body references
          (the H1 names the bug, the body opening paragraph cites it,
          and the Notes section describes the reproduction). The
          remaining 8 known_issues entries target sub_features outside
          this scenario file (widgets, dashboards, dialogs,
          formula-lines, add-new-column, io); coverage of those bugs
          by other scenarios in the section is a chain-level concern
          owned by Gate F under F-BUG-COVERAGE-01, not Gate A.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          No scenario opted out for effort or complexity reasons. The
          Notes Deferrals line reads None explicitly, and the three
          scenarios fully cover the bug's reproduction surface
          (Scenario 1 canonical QA repro, Scenario 2 five additional
          dispatch shapes, Scenario 3 suggestion-nav regression
          guard).
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No TODO add later or to be done in next phase deferrals
          appear anywhere in the body. The Notes Deferrals bullet is
          the explicit string None. The Coverage-map note records a
          STEP B fallback (coverage-map for powerpack not present at
          authoring time) but does not defer scenario work — it
          documents the absent gap-cross-check artifact only.
  d:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T13:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T15:42:00Z
    spec_runs:
      - spec: power-search-enter-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 30
        failure_keys: []
---

# PowerPack — Power Search ENTER-key dispatch (GROK-18656 regression)

Bug-focused regression scenario for `GROK-18656`: pressing **Enter**
in the Welcome View search input must safely dispatch the Power Search
engine across all 8 dispatch paths
(`powerSearch` in
`public/packages/PowerPack/src/search/power-search.ts#L26`:
view-search → regex-entity → project → JS-eval → function-call →
table-query → specific-widget → function-evaluation). Before the fix,
short or ambiguous queries (e.g., `QA`) caused the dispatch to throw
(`reported error: null`) because one of the 8 path branches did not
null-check its intermediate state. Each path must either produce a
result, produce a graceful empty-result, or surface a user-friendly
error — none may throw.

Atlas surface exercised: `powerpack.search.power` (the parallel free-
text power-search engine), `powerpack.search.power.dispatch`
(the 8-path top-level dispatch — primary regression surface),
`powerpack.search.power.init` (wires the engine to the platform input;
covered implicitly by every Enter-key dispatch on the welcome view),
`powerpack.welcome.view` (`welcomeView` autostart-immediate; the host
of the search input that this scenario types into),
`powerpack.welcome.suggestion-nav` (`suggestionMenuKeyNavigation` —
the Enter-on-highlighted-suggestion path that must remain functional
after the fix), and `powerpack.search` (the umbrella search
infrastructure).

## Setup

1. Open Datagrok with PowerPack installed (default platform load —
   `powerPackInit` runs at startup and `initSearch()` wires the
   power-search engine to the welcome-view search input).
2. Navigate to the Welcome view (Datagrok logo / `Home` icon in the
   sidebar). The Welcome view's search input must be visible and
   focused (per `welcomeView` autostart-immediate behavior).
3. Open the browser DevTools Console (F12) so console errors are
   visible while pressing Enter — the GROK-18656 reproduction
   surfaces as a thrown error visible in the console / a Datagrok
   error balloon.

## Scenarios

### Scenario 1: Pressing Enter on short query `QA` (canonical GROK-18656 reproduction)

1. **Focus the Welcome view search input.** Click into the search
   field at the top of the Welcome view.
2. **Type the query `QA`** (two uppercase letters — the canonical
   GROK-18656 reproduction query).
3. **Press Enter.** Power Search dispatches across the 8 paths.
4. **Verify no error fires.** No exception in the browser DevTools
   Console; no Datagrok red error balloon; no `reported error: null`
   message. The Welcome view does not crash or freeze.
5. **Verify a sensible dispatch result is rendered.** One of:
   - Matching entities surface in the suggestion menu / results
     pane (the view-search, regex-entity, function-call, table-query,
     or widget path produced a hit).
   - A graceful empty-results message is shown (none of the 8 paths
     produced a hit, but the dispatcher completed cleanly).
   - A user-friendly error balloon is shown (e.g. "no results
     found" or similar) — never an uncaught throw.

Expected:
- Pressing Enter on `QA` does NOT throw. This is the direct
  GROK-18656 invariant.
- The dispatch result is one of the three sensible outcomes
  enumerated in Step 5.

### Scenario 2: Pressing Enter across additional dispatch shapes (8-path exercise)

This scenario exercises multiple of the 8 dispatch paths within
`powerSearch` to guard against null-safety regressions in any single
branch (per atlas `powerpack.search.power.dispatch`: 8-path complex
code; any path may harbor a similar regression).

For each of the following queries, repeat the same 4-step probe:

1. **Focus the Welcome view search input** (clear any prior query).
2. **Type the query** (exact string as listed).
3. **Press Enter.**
4. **Verify no error fires and a sensible dispatch result is
   rendered** (per Scenario 1 Step 5 acceptance criteria).

Queries to exercise (one per row):

- **`a`** — single character; exercises the short-query path (the
  most aggressive null-safety surface; the same shape that triggered
  `QA` in GROK-18656).
- **`new`** — common word that overlaps the New Users sub-provider
  suggestion path; exercises the function-call / widget dispatch
  branches.
- **`user`** — username-style; exercises the regex-entity branch
  (entity name match against users) and the function-call branch.
- **`Project[0-9]+`** — regex-pattern-shaped entity name; exercises
  the regex-entity dispatch branch explicitly.
- **`1+1`** — function-evaluation expression; exercises the JS-eval
  / function-evaluation dispatch branch.

Expected (each query, independently):
- No exception in the DevTools Console; no `reported error: null`.
- The Welcome view remains responsive; the search input retains
  focus or a sensible follow-up view is loaded.
- A result, empty-results message, or user-friendly error is
  rendered — never a thrown exception.

### Scenario 3: Enter-on-highlighted-suggestion path is preserved (suggestion-nav regression guard)

The GROK-18656 fix must not break the `suggestionMenuKeyNavigation`
Enter-on-highlighted-suggestion path
(`public/packages/PowerPack/src/welcome-view.ts#L118` —
atlas `powerpack.welcome.suggestion-nav` declares
`press Enter on highlighted suggestion → click suggestion`).

1. **Focus the Welcome view search input.**
2. **Type a query that surfaces a suggestion menu** (e.g., type
   `dem` so the Demog dataset / sample data surfaces in the
   suggestion list).
3. **Press ArrowDown** to highlight the first suggestion in the
   menu (verify the highlight moves to the first suggestion —
   atlas declared interaction).
4. **Press ArrowDown again** to move highlight to the second
   suggestion (or ArrowUp to move back to the first — verify
   highlight navigation works in both directions).
5. **Press Enter on the highlighted suggestion.** The suggestion
   is activated (equivalent to clicking it): the corresponding
   entity is opened or the corresponding view is navigated to.

Expected:
- ArrowDown / ArrowUp navigation in the suggestion menu works
  as declared in atlas `powerpack.welcome.suggestion-nav`
  interactions.
- Enter on a highlighted suggestion activates that suggestion
  (clicks it) — the pre-existing behavior preserved by the
  GROK-18656 fix.
- No collision between the suggestion-nav Enter handler and the
  power-search Enter dispatch — the suggestion-nav path takes
  precedence when a suggestion is highlighted (per the atlas
  interaction declaration); the power-search dispatch fires only
  when no suggestion is highlighted.

## Notes

- **target_layer rationale.** `playwright`. The scenario exercises
  the Welcome view search input's Enter-key handler and the resulting
  Power Search dispatch — a UI surface backed by DOM keyboard events,
  rendered suggestion menus, and rendered result panes. Headless
  JS-API (`apitest`) cannot drive a real Enter-key dispatch through
  the welcome-view input element nor verify suggestion-menu
  rendering. No pixel-precision drag or native file picker is
  involved, so `manual-only` is not required.
- **coverage_type rationale.** `regression`. Bug-focused
  (`pyramid_layer: bug-focused` per Rule 3 — canonical GROK-18656
  reproduction surface); guards against re-regression on the Power
  Search Enter-key dispatch path. Not `smoke` (this is not a section
  golden path — the section's smoke is the top-level
  `add-new-column.md` per chain `ui_coverage_plan`). Not `edge`
  (the bug surfaces on a common entry path — typing a short query
  and pressing Enter — not on a boundary value). Not `perf`.
- **Pyramid layer.** `bug-focused` per Rule 3 — discriminator
  test: GROK-18656 fails Scenario 1 (and the `QA` row of Scenario 2)
  before the fix (`reported error: null` thrown on the dispatcher).
  After the fix, Scenarios 1-3 all pass. Atlas `coverage_type`
  on edge-case entries is `[]` (Phase 0 bootstrap), so this
  scenario's `coverage_type: regression` is derived from STEP E
  heuristics (general bug-focused coverage of a common feature
  shape) rather than from an atlas `edge_cases[]` entry — no
  cross-check mismatch.
- **Atlas sub_features traceability.**
  - `powerpack.search.power` — the parallel free-text power-search
    engine (`public/packages/PowerPack/src/search/power-search.ts#L1`);
    the umbrella sub_feature exercised by every Enter-key dispatch
    in this scenario.
  - `powerpack.search.power.dispatch` — `powerSearch(s, host,
    inputElement)`
    (`public/packages/PowerPack/src/search/power-search.ts#L26`);
    the 8-path top-level dispatch that GROK-18656 broke. Primary
    regression surface.
  - `powerpack.search.power.init` — `initSearch()`
    (`public/packages/PowerPack/src/search/power-search.ts#L21`);
    wires the engine to the platform input. Covered implicitly by
    every Enter-key dispatch on the welcome view (if `initSearch`
    had not run, the input would not dispatch through `powerSearch`).
  - `powerpack.welcome.view` — `welcomeView`
    (`public/packages/PowerPack/src/package.ts#L134`); the
    autostart-immediate home view that hosts the search input.
  - `powerpack.welcome.suggestion-nav` —
    `suggestionMenuKeyNavigation`
    (`public/packages/PowerPack/src/welcome-view.ts#L118`); the
    ArrowDown / ArrowUp / Enter-on-highlighted-suggestion handler.
    Scenario 3 guards against accidental regression of this path
    by the GROK-18656 fix.
  - `powerpack.search` — the umbrella search infrastructure
    (`public/packages/PowerPack/src/package.ts#L219`).
- **Related bug.** `GROK-18656` (p2, regression-risk).
  Reproduction: open Welcome view → type `QA` → press Enter →
  error fires (`reported error: null`). Expected: Power Search
  engine must safely handle Enter-key dispatch from the Welcome
  View search input across all 8 dispatch paths. Short queries
  like `QA` must dispatch gracefully (results, empty results, or
  user-friendly error — never throw).
- **Chain context.** This scenario is the section's bug-focused
  witness for GROK-18656 — Critic F's
  `bug_focused_candidates[]` entry for the bug (added in
  `cycle-2026-05-20-powerpack-coverage`) had empty `spans[]`
  pending this scenario's authoring. The chain's
  `order_from_files[]` is updated to include this scenario at
  the end of the section's order.
- **Deferrals.** None. The three scenarios fully cover the bug's
  reproduction surface (Scenario 1: canonical `QA` repro;
  Scenario 2: 5 additional dispatch shapes exercising multiple
  of the 8 paths; Scenario 3: preserved suggestion-nav Enter
  path). No deferral required.
- **Coverage map.** Coverage map for PowerPack
  (`references/coverage-map/powerpack.yaml`) is not present at
  authoring time — gap-vs-coverage cross-check skipped per
  STEP B fallback. The Critic F coverage-gap dispatch
  (`gap: bug-uncovered :: GROK-18656`) drove this scenario's
  authoring directly.
