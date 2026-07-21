---
feature: histogram
realizes_atlas:
  - histogram.cp.range-filter-bounds
realizes:
  - viewers.histogram
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-18948
    status: fixed
  - id: GROK-19581
    status: fixed
  - id: GROK-19760
    status: fixed
  - id: github-2329
    status: fixed
  - id: github-2301
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: "df.filter.trueCount reflects only the rows whose AGE falls within
      [30, 60] — less than the full row count."
  - anchor: "Scenario 1 Step 5"
    expectation: "df.filter.trueCount is unchanged after enabling Split Stack (the
      filter survives stacking, GROK-18948)."
  - anchor: "Scenario 1 Step 7"
    expectation: "Enabling Normalise to Filter with the sub-range active does not
      raise a new console error (github-2329)."
  - anchor: "Scenario 1 Step 8"
    expectation: "Narrowing Min to 40 via the range input field further reduces
      df.filter.trueCount and the viewer re-renders without error
      (github-2329)."
  - anchor: "Scenario 2 Step 2"
    expectation: "Setting Max below Min (e.g. Max = 20 when Min = 40) does not
      crash: df.filter.trueCount drops below the prior valid-range count and no
      new console error is raised (GROK-19581, GROK-19760). OBSERVED BEHAVIOR
      (recon 2026-07-21, demog.csv), NEEDS PRODUCT CONFIRMATION: the build
      applies the inverted range verbatim — it neither refuses nor corrects it —
      collapsing the filter to almost no rows with no user-facing signal. The
      tickets guarantee only no-crash; whether silently accepting an inverted
      range is intended is a product question. The spec asserts the observed
      collapse + no-error floor, not a rejection the product does not perform."
  - anchor: "Scenario 2 Step 3"
    expectation: "Setting Min to a value below the column minimum (e.g. -999) is
      applied verbatim, not clamped, and widens the lower bound to the data
      extent: df.filter.trueCount rises above the prior count, the effective range
      stays non-inverted (min <= max), and no new console error is raised
      (GROK-19581)."
  - anchor: "Scenario 2 Step 4"
    expectation: "Setting Max to a value above the column maximum (e.g. 999) is
      applied verbatim, not clamped, and widens the upper bound to the data
      extent: df.filter.trueCount rises above the prior count, the effective range
      stays non-inverted (min <= max), and no new console error is raised
      (GROK-19760)."
  - anchor: "Scenario 3 Step 2"
    expectation: "df.filter.trueCount returns to the full row count after restoring
      Min and Max to the full column range."
  - anchor: "Scenario 3 Step 3"
    expectation: "Disabling filtering leaves the viewer in a no-error state (canvas
      repaint teardown floor)."
realized_as:
  - histogram-range-filter-bounds-spec.ts
---

# Histogram — Range Filtering and Bound Validation

## Setup

1. Close all open tables and viewers.
2. Open `System:DemoFiles/demog.csv`.
3. Add a Histogram viewer to the table view.
4. In the properties panel (Context Panel > Value), set **Value** to `AGE`.
5. In the properties panel (Context Panel > Filter), enable **Show Range Inputs** (exposes the Min / Max text fields).
6. Enable filtering by turning on **Filter** in the properties panel.

## Scenarios

### Scenario 1: Sub-range filter survives Split Stack; Normalise to Filter path

Steps:
1. In the Min field, type `30`; in the Max field, type `60`. Press Enter after each.
2. Read `df.filter.trueCount` via the JS API helper (e.g. `grok.shell.t.filter.trueCount`).
3. Verify `df.filter.trueCount` reflects only the rows whose AGE value falls within [30, 60] — the count must be less than the total row count of demog.csv.
4. In the properties panel under Category, set **Split Column Name** to `SEX` and enable **Split Stack**.
5. Verify `df.filter.trueCount` is unchanged after enabling Split Stack — the range filter must survive stacking (GROK-18948). Disable Split Stack afterwards (no-error teardown).
6. Clear **Split Column Name** (set to None).
7. Enable **Normalise to Filter** — verify no new console error is raised and the viewer re-renders (github-2329).
8. Set **Min** to `40` via the range input field (keep Max at 60). Verify `df.filter.trueCount` decreases further (more rows excluded) and no new console error is raised (github-2329 — normalise-on-range-change path).

Expected:
- df.filter.trueCount reflects only the rows within the [30, 60] AGE sub-range.
- df.filter.trueCount is unchanged after enabling Split Stack (the range filter survives stacking, GROK-18948).
- Enabling Normalise to Filter and then narrowing the range via the Min input does not raise a new console error and the viewer re-renders (github-2329).

### Scenario 2: Invalid range inputs are rejected without crashing

Steps:
1. Set **Min** to `40` and **Max** to `60` via the range input fields (establish a known-valid sub-range).
2. Set **Max** to `20` (below the current Min of 40) — a Min >= Max violation. Verify: the inverted range is applied (not refused) and collapses the sub-range so `df.filter.trueCount` drops below the prior valid-range count, with no new console error (GROK-19581, GROK-19760).
3. Restore a valid Max (e.g. `60`), then set **Min** to `-999` (below the column minimum). Verify the value is applied verbatim (not clamped) and widens the lower bound to the data extent — `df.filter.trueCount` rises, the effective range stays non-inverted (min <= max), and no new console error is raised (GROK-19581).
4. Set **Min** back to `40`, then set **Max** to `999` (above the column maximum). Verify the value is applied verbatim (not clamped) and widens the upper bound to the data extent — `df.filter.trueCount` rises, the effective range stays non-inverted (min <= max), and no new console error is raised (GROK-19760).

Expected:
- Setting Max below Min applies the inverted range and collapses the filter (trueCount drops); the build does not refuse or correct the inversion, but raises no console error (GROK-19581, GROK-19760).
- Setting Min below the column minimum is applied verbatim (not clamped) and widens the lower bound to the data extent (trueCount rises), effective range non-inverted, no new console error (GROK-19581).
- Setting Max above the column maximum is applied verbatim (not clamped) and widens the upper bound to the data extent (trueCount rises), effective range non-inverted, no new console error (GROK-19760).

### Scenario 3: Round-trip — restore full range and disable filtering

Steps:
1. Disable **Normalise to Filter** (revert to default).
2. Set **Min** and **Max** back to the full AGE column range (the column's actual min and max as reported by `df.getCol('AGE').min` / `.max`). Verify `df.filter.trueCount` returns to the full row count of demog.csv.
3. Disable **Filter** in the properties panel. Verify the viewer remains in a no-error state and the histogram canvas repaints (no-error teardown floor).

Expected:
- df.filter.trueCount returns to the full demog.csv row count after restoring Min and Max to the column range (round-trip).
- Disabling filtering leaves the viewer in a no-error, cleanly repainted state.
