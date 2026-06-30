---
created: 2026-05-04T13:50:00.000Z
title: Reorganize QC dashboard layout for navigability
area: ui
files:
  - packages/Proteomics/src/viewers/qc-dashboard.ts:135-143
---

## Problem

The QC dashboard tiles 7 viewers around the central grid for a total of 8
visible panels:

1. Central data grid (the source DataFrame)
2. MA Plot (scatter)
3. MA Trend (scatter)
4. CV Plot (scatter)
5. Sample Correlation (CORR_PLOT)
6. Missing Values (grid)
7. Missing % per Sample (bar chart)
8. Intensity Distributions (box plot)

They're docked into a single workspace with no grouping by purpose, no visual
hierarchy, and labels that don't tell a story — "MA Plot / MA Trend / CV Plot"
sit adjacent without explaining how they relate. To a viewer encountering it
for the first time (e.g. the Cytokinetics audience the 2026-05-04 demo is being
prepared for), it reads as a wall of plots rather than a guided diagnostic
flow. On laptop-projector resolutions the panels are also small enough that
several plots become unreadable.

## Solution

Investigate and propose. Options to weigh:

1. **Group by analytical purpose using `DG.TabControl`.** Four natural tabs:
   - *Distributions* — Intensity Distributions box plot, Missing % per Sample
   - *Variance* — CV Plot, Sample Correlation
   - *MA Diagnostics* — MA Plot, MA Trend
   - *Missingness* — Missing Values grid, Missing % per Sample bar
   This collapses 7 panels into 4 named tabs and gives each plot enough screen
   space to read at projection resolution.

2. **Hero + drawer.** Pick 2–3 highest-value viewers (likely Intensity
   Distributions + Sample Correlation + MA Plot) for the primary docked
   region; dock the rest in a secondary panel that can be collapsed.

3. **Add inline captions / section headers** so each group communicates
   *what question it answers* (e.g. "Are the samples comparable?" above the
   distribution plots, "Where are values missing?" above the missingness pair).

4. **Title the dashboard window itself** with the source DataFrame name so the
   user knows what experiment these QC panels describe — relevant given
   Phase 11 already wired DataFrame naming for other viewers.

Whichever path is taken, the acceptance bar is: a viewer who has never seen
the package can, within ~10 seconds, identify which panel answers
"are my samples comparable" without being told.

## Related

- `packages/Proteomics/.planning/todos/pending/2026-05-04-add-progress-indication-to-qc-dashboard.md`
  — the perceived-hang complaint that surfaced this dashboard's UX issues
  in the same demo prep
- `packages/Proteomics/.planning/todos/pending/2026-03-10-improve-chart-titles-and-axis-labels-across-all-viewers.md`
  — adjacent UX item; consider tackling together since both touch the same viewers
