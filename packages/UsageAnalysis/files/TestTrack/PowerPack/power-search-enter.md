---
feature: powerpack
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [powerpack.cp.power-search-dispatch, welcome-view-power-search]
realizes: []
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

Regression test for GROK-18656: pressing **Enter** in the Welcome View
search input must safely dispatch the Power Search engine across all
of its lookup paths (view search, entity name match, project match, JS
expression evaluation, function call, table query, widget match,
function evaluation). Before the fix, short or ambiguous queries (e.g.
`QA`) caused the dispatch to throw an error, because one of the paths
did not check for a missing intermediate value. Each path must either
produce a result, a graceful empty result, or a user-friendly error —
never an uncaught exception.

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

This scenario exercises several more of the search dispatch paths to
guard against null-safety regressions in any of them.

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

The GROK-18656 fix must not break the existing behavior where pressing
Enter on a highlighted suggestion activates it.

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

- **Related bug.** GROK-18656: opening the Welcome view, typing `QA`,
  and pressing Enter threw an error (`reported error: null`). The
  Power Search engine must handle Enter-key dispatch gracefully for
  short queries like this — producing results, an empty-results
  message, or a user-friendly error, but never throwing.
