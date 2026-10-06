---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [projects.cp.url-parameterized-share]
realizes: [views.projects]
realized_as:
  - project-url-spec.ts
pyramid_layer: integration
ui_coverage_responsibility:
  - context-panel-links-url-copy
  - new-tab-open-url
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/projects/project-url.md
migration_date: 2026-05-20
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.projects.buildVariantsComposite
unresolved_ambiguities:
  - click-the-project-single-click-vs-double-click
  - context-panel-links-location
  - url-format-query-parameters
  - new-tab-vs-incognito-tab
  - order-vs-dependency-contradiction
  - source-text-correction-recurrence
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      The existing `project-url-spec.ts` exercises only the `demog`
      representative source (1 of the 4 project variants) via direct URL
      navigation; the copied-with-link / copied-with-clone /
      personal-view-customizations variants are not deep-linked in the spec.
      The URL build/apply/shell-open contract is source-agnostic, so the
      single-variant walk preserves the invariant. NOTE: gate_verdicts.b is
      FAIL (spec unstable, [B-RUN-PASS, B-STAB-01]) — a separate open issue,
      not addressed by this scope reduction.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
  d:
    verdict: PASS
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: batch-6.6-pilot-2026-05-20-projects-project-url
    timestamp: 2026-05-20T12:35:00Z
    spec_runs:
      - spec: project-url-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 70
        failure_keys: [B-RUN-PASS, B-STAB-01]
---

# Project URL — open a project from its link

Copies the link of a saved project from the **Context Panel**, opens
it in a new browser tab, and checks that the right project opens. This
is done for the original project and for its copies.

## Setup

1. Log in as the test user.
2. The projects from `projects-copy-clone.md` must exist: `copyClone`,
   `copyCloneLink` and `copyCloneClone`. If they do not, run that test
   first.

## Scenario

1. **The original project.**
   - Go to **Browse > Dashboards**.
   - Type `copyClone` into the search box.
   - Click the `copyClone` tile once.
   - In the **Context Panel**, under **Details**, click **Links...**.
   - In the **Links to copyClone** dialog, click the copy icon next to
     **URL**.
   - Open a new browser tab.
   - Paste the URL into the address bar.
   - **Verify:** the URL is the server address, then `/p/`, your login
     and `.CopyClone`.
   - Press Enter.
   - **Verify:** the new tab shows the grid, the bar chart and the
     scatter plot.
   - **Verify:** there is no line chart.
   - **Verify:** the grid is sorted by `AGE` ascending.
   - **Verify:** the `DIS_POP` column is hidden.
   - **Verify:** no error balloon appears.
   - Close the tab.

2. **The copy with link.**
   - Go to **Browse > Dashboards**.
   - Type `copyCloneLink` into the search box.
   - Click the `copyCloneLink` tile once.
   - In the **Context Panel**, under **Details**, click **Links...**.
   - In the **Links to copyCloneLink** dialog, click the copy icon next to
     **URL**.
   - Open a new browser tab.
   - Paste the URL into the address bar.
   - **Verify:** the URL is the server address, then `/p/`, your login
     and `.CopyCloneLink`.
   - Press Enter.
   - **Verify:** the new tab shows the grid, the line chart, the bar
     chart and the scatter plot.
   - **Verify:** no error balloon appears.
   - Close the tab.

3. **The copy with clone.**
   - Go to **Browse > Dashboards**.
   - Type `copyCloneClone` into the search box.
   - Click the `copyCloneClone` tile once.
   - In the **Context Panel**, under **Details**, click **Links...**.
   - In the **Links to copyCloneClone** dialog, click the copy icon next to
     **URL**.
   - Open a new browser tab.
   - Paste the URL into the address bar.
   - **Verify:** the URL is the server address, then `/p/`, your login
     and `.CopyCloneClone`.
   - Press Enter.
   - **Verify:** the new tab shows the grid, the histogram, the bar
     chart and the scatter plot.
   - **Verify:** no error balloon appears.
   - Close the tab.

4. **Cleanup.**
   - In **Browse > Dashboards**, right-click the `copyClone` tile and
     choose **Delete Project**.
   - Click **DELETE**.
   - Wait until the dialog closes.
   - Repeat for `copyCloneLink` and `copyCloneClone`.

## Expected results

- The link from **Links...** opens exactly that project in a new tab.
- Copies open as themselves, not as the original.
