---
created: 2026-05-04T13:42:23.818Z
title: Add progress indication to QC dashboard
area: ui
files:
  - packages/Proteomics/src/viewers/qc-dashboard.ts:21-55
---

## Problem

Running the QC dashboard takes a noticeable amount of time (MA computation, loess
trend, per-group CV, correlation matrix, then building and tiling 4–5 viewers),
and from the user's perspective there is no visual indication that anything is
happening. Surfaced during the 2026-05-04 Cytokinetics demo prep — the long pause
between menu click and dashboard appearance reads as "broken or hung" to a viewer
who hasn't seen the package before.

`openQcDashboard` already creates a `DG.TaskBarProgressIndicator` at line 22:

```ts
const pi = DG.TaskBarProgressIndicator.create('Creating QC dashboard...');
```

…so something is being signalled. The complaint is that this signal isn't
prominent enough to register: the task bar indicator is a thin bar at the bottom
of the screen, easy to miss, and the message never updates as the work proceeds.

## Solution

Investigate and address in roughly this order:

1. **Verify the existing PI is firing.** Reproduce the slow QC run and confirm
   the task bar shows "Creating QC dashboard..." for the duration. If it doesn't
   appear at all, fix the existing call before adding more.

2. **Make the existing indicator phase-aware.** Call `pi.update(...)` between the
   compute steps so the user sees the stages move forward:
   - "Computing MA values..."
   - "Computing loess trend..."
   - "Computing per-group CV..."
   - "Building dashboard..."
   This usually flips perception from "hung" to "working" without any new UI.

3. **Consider an inline indicator** on the dashboard panel itself (a centered
   spinner + status text inside the docked area where the viewers will land), so
   the feedback is co-located with the place the user is looking. Reference the
   diagnosis in `packages/Proteomics/.planning/debug/showheatmap-hangs.md` —
   `showHeatmap` has the same class of problem and the proposed fix there
   (TaskBarProgressIndicator + inline feedback) likely applies in the same shape.

Acceptance: a viewer with no prior context should be able to tell the QC
dashboard is working within ~1 second of clicking the menu item, and should
never see > 3 seconds of UI silence during the build.
