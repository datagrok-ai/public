---
feature: notebooks
original_path: public/packages/UsageAnalysis/files/TestTrack/notebooks/browser.md
migration_date: 2026-06-17
source_text_fixes:
  - removed-order-json-footer
  - demog-apply-to-condition-promoted-to-setup-note
candidate_helpers: []
unresolved_ambiguities:
  - filter-templates-icon-selector-unclear
  - apply-to-demog-requires-demog-table-open-in-fixture
scope_reductions: []
target_layer: playwright
realized_as:
  - browser-spec.ts
coverage_type: regression
priority: p0
realizes: [browse-and-open-html]
produced_from: migrated
related_bugs:
  - GROK-11693
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T18:00:00Z
    failure_keys: []
    review_round: 1
  d:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T17:30:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T16:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T21:30:00Z
    failure_keys: []
  b:
    verdict: FAIL
    cycle_id: 2026-06-17-notebooks-migrate-01
    timestamp: 2026-06-17T22:30:00Z
    spec_runs:
      - spec: browser-spec.ts
        result: failed
        attempts: 3
        duration_seconds: 105
        failure_keys: [ B-RUN-PASS, B-STAB-01 ]
        run_mode: headless-cold
---

# Notebook Browser

Exercises the Notebooks browser view: navigation, filter templates, context panel
accordion, Apply-to (requires matching open table), back-navigation, and New Notebook
creation from the browser ribbon.

**Chain dependency:** Requires the `demog` notebook to already exist on the server
(produced by `create.md`). Step 4 (Apply to) additionally requires `demog.csv` to be
open in the platform before navigating to the Notebooks browser.

See: `public/help/compute/jupyter-notebook.md#Notebooks browser`

## Setup

- Server has a notebook named `demog` persisted (created by the `create.md` scenario).
- `System:DemoFiles/demog.csv` is open as an active table so that the notebook is
  applicable and the **Apply to > demog** context-menu item appears.

## Scenarios

### S1 — Navigate to the Notebooks browser

**Steps:**

1. Go to **Browse > Platform > Notebooks** (or via **ML | Notebooks | Browse Notebooks**).

**Expected result:**

- The Notebooks browser view opens.
- At least the `demog` notebook card is visible in the gallery.

---

### S2 — Filter templates

**Steps:**

2. Near the search field, locate the **Filter templates** icon and click it to open the
   filter-template picker; apply at least one template to filter the list of notebooks.

**Expected result:**

- The notebook list is filtered according to the selected template.
- Deselecting or clearing the filter restores the full list.

> **Ambiguity:** The exact selector for the "Filter templates" icon is not documented
> in the refdoc; it may be the same as the `[name="icon-filter"]` Toggle filters
> control or a distinct control. Needs recon to confirm the selector before
> automation.

---

### S3 — Context panel accordion tabs

**Steps:**

3. Click the `demog` notebook card to select it and inspect all tabs on the **Context
   Panel** (Details, Actions, Activity, Sharing, Chats).

**Expected result:**

- The Context Panel opens and displays the `demog` notebook's Details accordion (name,
  description, Created / Modified dates, tags).
- All accordion sections (Details, Actions, Activity, Sharing, Chats) are accessible
  and render without console errors.
- In particular, the **Sharing** tab renders the notebook's sharing info without a
  NullError. (Regression guard for GROK-11693: error in "Sharing" tab on Property Panel
  for a notebook created via Open in Notebook.)

---

### S4 — Apply to open table

**Steps:**

4. Right-click the `demog` notebook and select **Apply to > demog** from the context
   menu (the sub-menu entry appears because the `demog.csv` table is open and matches
   the notebook's linked tables).

**Expected result:**

- The `demog` notebook is executed against the `demog` table.
- The rendered HTML output view opens showing the executed notebook results.
- No console errors during apply.

> **Ambiguity:** The test fixture must open `demog.csv` before navigating to the
> Notebooks browser; the original scenario states `(available if the demog dataset is
> open)` as an inline condition. This must be a fixture step run before the scenario.
>
> Note: Executing the notebook against the table requires a live Docker container for
> Jupyter notebook execution. If the container is not running this step may fail with
> an infra-level error; treat as environmental flakiness in that case.

---

### S5 — Back-navigate to the Notebook Browser

**Steps:**

5. Navigate back to the **Notebook Browser** (e.g., via the browser back-button or
   **ML | Notebooks | Browse Notebooks**).

**Expected result:**

- The Notebooks browser view re-opens.
- The notebook list is restored to the pre-filter state (no active filter from S2).

---

### S6 — New Notebook from ribbon

**Steps:**

6. In the browser ribbon, click **New Notebook** and wait for the new notebook view to open.

**Expected result:**

- A new empty notebook is created on the server and opens in the editor view.
- The new notebook appears in the Notebooks browser when navigating back.
- No console errors during creation.

---

## Notes

- **GROK-11693** (Sharing tab NullError for a notebook created via Open in Notebook) is
  directly exercised in Step 3. This scenario is the primary regression guard for that bug.
- See `public/help/compute/jupyter-notebook.md` for the help-page navigation reference.

---
{
  "order": 3
}
