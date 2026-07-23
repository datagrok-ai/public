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
  - anchor: "S1: Min=30/Max=60 filters to AGE sub-range (< full row count)"
    expectation: "df.filter.trueCount reflects only the rows whose AGE falls within
      [30, 60] — less than the full row count."
  - anchor: "S1: filter.trueCount unchanged after Split Stack (GROK-18948)"
    expectation: "df.filter.trueCount is unchanged after enabling Split Stack (the
      filter survives stacking, GROK-18948)."
  - anchor: "S1: Normalise to Filter raises no error (github-2329)"
    expectation: "Enabling Normalise to Filter with the sub-range active does not
      raise a new console error (github-2329)."
  - anchor: "S1: Min=40 narrows filter further, no error (github-2329)"
    expectation: "Narrowing Min to 40 via the range input field further reduces
      df.filter.trueCount and the viewer re-renders without error
      (github-2329)."
  - anchor: "S2: establish valid Min=40/Max=60 sub-range"
    expectation: "After typing Min = 40 and Max = 60 into the range input fields,
      the filtered row count (df.filter.trueCount) lands strictly between zero and
      the full row count — a known-valid sub-range is in effect before the
      invalid-bound steps."
  - anchor: "S2: Max below Min collapses filter without crash (GROK-19581, GROK-19760)"
    expectation: "Setting Max below Min (e.g. Max = 20 when Min = 40) does not
      crash: the filter card shows the validation warning \"max should be
      greater than min\", df.filter.trueCount drops below the prior valid-range
      count, and no new console error is raised (GROK-19581, GROK-19760). This
      is intended behavior (confirmed by the operator): the UI flags the
      inverted range to the user and still applies the collapsed sub-range. The
      spec asserts the collapse + no-error floor; the warning text is a
      DOM-assertable signal available if a stronger check is wanted later."
  - anchor: "S2: Min below column min widens range without crash (GROK-19581)"
    expectation: "Setting Min to a value below the column minimum (e.g. -999) is
      applied verbatim, not clamped, and widens the lower bound to the data
      extent: df.filter.trueCount rises above the prior count, the effective range
      stays non-inverted (min <= max), and no new console error is raised
      (GROK-19581)."
  - anchor: "S2: Max above column max widens range without crash (GROK-19760)"
    expectation: "Setting Max to a value above the column maximum (e.g. 999) is
      applied verbatim, not clamped, and widens the upper bound to the data
      extent: df.filter.trueCount rises above the prior count, the effective range
      stays non-inverted (min <= max), and no new console error is raised
      (GROK-19760)."
  - anchor: "S3: restore full range returns filter to full row count (round-trip)"
    expectation: "df.filter.trueCount returns to the full row count after restoring
      Min and Max to the full column range."
  - anchor: "S3: disable Filter leaves viewer no-error (teardown floor)"
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
2. Read the filtered row count (the number of rows currently passing the filter, shown in the table's status bar).
3. Verify the filtered row count reflects only the rows whose AGE value falls within [30, 60] — the count must be less than the total row count of demog.csv.
4. In the properties panel under Category, set **Split Column Name** to `SEX` and enable **Split Stack**.
5. Verify the filtered row count is unchanged after enabling Split Stack — the range filter must survive stacking (GROK-18948). Disable Split Stack afterwards (no-error teardown).
6. Clear **Split Column Name** (set to None).
7. Enable **Normalise to Filter** — verify no new console error is raised and the viewer re-renders (github-2329).
8. Set **Min** to `40` via the range input field (keep Max at 60). Verify the filtered row count decreases further (more rows excluded) and no new console error is raised (github-2329 — normalise-on-range-change path).

Expected:
- The filtered row count reflects only the rows within the [30, 60] AGE sub-range.
- The filtered row count is unchanged after enabling Split Stack (the range filter survives stacking, GROK-18948).
- Enabling Normalise to Filter and then narrowing the range via the Min input does not raise a new console error and the viewer re-renders (github-2329).

### Scenario 2: Out-of-range and inverted inputs are handled without crashing

Steps:
1. Set **Min** to `40` and **Max** to `60` via the range input fields (establish a known-valid sub-range).
2. Set **Max** to `20` (below the current Min of 40) — a Min >= Max violation. Verify: the filter card shows the validation warning "max should be greater than min", the inverted range still applies and collapses the sub-range so the filtered row count drops below the prior valid-range count, with no new console error (GROK-19581, GROK-19760). Intended behavior — the UI warns the user and still applies the collapsed range.
3. Restore a valid Max (e.g. `60`), then set **Min** to `-999` (below the column minimum). Verify the value is applied verbatim (not clamped) and widens the lower bound to the data extent — the filtered row count rises, the effective range stays non-inverted (min <= max), and no new console error is raised (GROK-19581).
4. Set **Min** back to `40`, then set **Max** to `999` (above the column maximum). Verify the value is applied verbatim (not clamped) and widens the upper bound to the data extent — the filtered row count rises, the effective range stays non-inverted (min <= max), and no new console error is raised (GROK-19760).

Expected:
- Setting Max below Min shows the "max should be greater than min" warning and still applies the inverted range, collapsing the filter (the filtered row count drops), with no console error (GROK-19581, GROK-19760) — intended behavior.
- Setting Min below the column minimum is applied verbatim (not clamped) and widens the lower bound to the data extent (the filtered row count rises), effective range non-inverted, no new console error (GROK-19581).
- Setting Max above the column maximum is applied verbatim (not clamped) and widens the upper bound to the data extent (the filtered row count rises), effective range non-inverted, no new console error (GROK-19760).

### Scenario 3: Round-trip — restore full range and disable filtering

Steps:
1. Disable **Normalise to Filter** (revert to default).
2. Set **Min** and **Max** back to the full AGE column range (the column's actual minimum and maximum values). Verify the filtered row count returns to the full row count of demog.csv.
3. Disable **Filter** in the properties panel. Verify the viewer remains in a no-error state and the histogram canvas repaints (no-error teardown floor).

Expected:
- The filtered row count returns to the full demog.csv row count after restoring Min and Max to the column range (round-trip).
- Disabling filtering leaves the viewer in a no-error, cleanly repainted state.

## Automation notes

- The Min / Max range text inputs are `.d4-filter-input-min` / `.d4-filter-input-max` inside the viewer root `[name="viewer-Histogram"]`; they appear once **Show Range Inputs** is on, and a typed value is committed with Enter.
- "Filtered row count" is read as `grok.shell.tv.dataFrame.filter.trueCount`; the full-range restore in Scenario 3 uses the column's actual extent from `df.getCol('AGE').min` / `.max`.
- Scenario 1: the [30, 60] sub-range check tolerates a one-row bin-edge difference between the filtered count and the number of AGE values inside the range.
