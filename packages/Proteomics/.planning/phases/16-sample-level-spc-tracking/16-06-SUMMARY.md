---
phase: 16
plan: 06
wave: 3
status: complete
completed: 2026-06-08
---

# Plan 16-06 — SPC Dashboard viewer + sidebar + Define-baseline modal + drill-down

## What's Built

### `src/viewers/spc-dashboard.ts` (~530 LOC, new)

| Public export | Role |
|---------------|------|
| `openSpcDashboard()` | Orchestrator. Loads runs, builds the instrument picker (suffix `(baseline locked)` / `(no baseline)`), caches the last-picked instrument under `userDataStorage` key `spc.last_picked_instrument`, then opens a sibling TableView via `grok.shell.addTableView` (NEVER docks onto the protein df's view). For baseline-less instruments the right pane carries the empty-state banner + `Define baseline...` CTA. With a baseline locked, four I-charts + four MR-charts are docked. |
| `createSpcChartPanel(runsDf, baseline, metric)` | Pure factory returning `{iChart, mrChart}`. I-chart gets UCL/CL/LCL formula lines (Plan 16-01 spike PASS path: `style: 'dashed'` for the bounds; if the spike returns FAIL at runtime, the line additions are wrapped in try/catch so the chart still renders — fallback to the constant-column overlay is captured as a follow-up). MR-chart gets UCL = D4 × mean(MR) + CL = mean(MR); LCL = 0 is the X-axis per UI-SPEC. |
| `computeControlLines(baseline, metric)` | Pure math used by `SPC:formula_lines` test AND by the I-chart factory above. Returns `{ucl, cl, lcl}` or NaN-bearing lines when the baseline metric is unavailable (so the factory can skip the addLine call). |
| `resolveDrillDown(row, deps)` | Pure async helper used by `SPC:drilldown_resolved` / `SPC:drilldown_missing` tests AND by the live click handler. `deps.projectsFind` is injectable so tests can mock without a server. Returns `{kind: 'opened', id}` on success and `{kind: 'toast', message: ...}` for both the null-id and find-returns-null cases. |
| `showDefineBaselineDialog(instrumentId, runs, onLocked, existing?)` | Modal owning the Define vs. Rebuild flow. First-7-checked default for Define; existing.included_run_ids checked for Rebuild. Iterate-3σ toggle + 4×8 Nelson rule grid + live false-alarm-rate disclosure. Per-metric outlier removal (UNION exclusion across the 4 metrics — D-04 conservative policy) followed by `computeBaselineStats` over retained runs. `onLocked(baseline)` callback persists via `saveBaseline`. |

### `src/package.ts`

Plan 16-05's stub body replaced with `await openSpcDashboard()`. The `Proteomics | Visualize | SPC Dashboard...` menu item now opens the real dashboard.

### Behind the scenes

- **MR columns** (`MR_<metric>`) added via `df.columns.addNewFloat(...)` per memory `feedback_dg_column_bulk_init` (cheap at SPC's ~104 rows/year but consistent with the bulk-init guidance); row 0 explicitly set to null per memory `feedback_dg_column_init_null_sentinel`.
- **Drill-down subscription** uses `tv.dataFrame.onCurrentRowChanged` (single-producer per CONTEXT.md anti-pattern note — no `CampaignSelectionBus`).
- **Sidebar grid** wires each `ui.input.bool` to an `onChanged` subscription that updates `baseline.rules_enabled` in memory, persists via `saveBaseline`, and re-evaluates Nelson rules over the whole slice via a local `recomputeStatusesInPlace` writeback — point colors update without re-opening the view.
- **`userDataStorage`** soft-deprecation comment links to `subcellular-location.ts:246` for the precedent and flags v1.5 migration to `grok.userSettings`.

## RED → GREEN

| Test | Asserts |
|------|---------|
| `SPC:dashboard_renders` | `createSpcChartPanel` returns non-null |
| `SPC:formula_lines` | `computeControlLines` returns `{ucl: 22.3+1.5, cl: 22.3, lcl: 22.3-1.5}` |
| `SPC:drilldown_resolved` | `resolveDrillDown` returns `{kind: 'opened', id: 'proj-abc'}` when projectsFind succeeds |
| `SPC:drilldown_missing` | `resolveDrillDown` returns `{kind: 'toast', message: ...'Source for ...r-zzz...}` when projectsFind returns null |

## Acceptance gate verification

| Gate | Result |
|------|--------|
| `npm run build` exit 0 | ✓ (3 perf warnings only) |
| `openSpcDashboard` export | 1 |
| `createSpcChartPanel` export | 1 |
| `showDefineBaselineDialog` export | 1 |
| `grok.shell.addTableView` | 1 |
| `grok.dapi.userDataStorage` | 2 (get + put) |
| `grok.dapi.projects.find` | 1 |
| `CampaignSelectionBus` | 0 (anti-pattern guard) |
| `addLine` formula-line calls | 7 (3 I-chart × +CL/UCL guards + 2 MR-chart) |
| `WHEELER_D4` | 2 (import + usage) |
| `openSpcDashboard` in package.ts | 2 (import + call) |

## Threat mitigations

- T-16-V4-ACL — `grok.dapi.projects.find(id)` returns null when ACL denies access; the resolver surfaces the same biologist-readable toast as the missing-source case.
- T-16-02 — Already mitigated by Plan 03's `loadBaseline` (try/catch on parse). A separate "baseline file corrupt" toast is a v1.5 hardening — for v1.4 a corrupt file is treated as no baseline (empty-state banner appears).

## Deferred (Plan 16-07)

`aggregateParetoCounts` for the P2 Pareto panel. The rules-tripped tag is set unconditionally by Plan 16-02's `setSpcStatus` so the Pareto can land as a v1.4.1 hotfix without unwinding any of the work here.

## Key files

```yaml
key-files:
  created:
    - src/viewers/spc-dashboard.ts
  modified:
    - src/package.ts
```

## Self-Check: PASSED
