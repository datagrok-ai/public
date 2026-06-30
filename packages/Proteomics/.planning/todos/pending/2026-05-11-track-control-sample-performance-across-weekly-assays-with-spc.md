---
created: 2026-05-11T14:10:00.000Z
title: Track control sample performance across weekly assays with SPC
area: analysis
resolves_phase: 16
resolves_requirements: SPC-01, SPC-02, SPC-03, SPC-04, SPC-05, SPC-06, SPC-07, SPC-08
files:
  - packages/Proteomics/src/viewers/qc-dashboard.ts
  - packages/Proteomics/src/analysis/experiment-setup.ts
---

## Problem — new use case

The screening campaign runs the **same assay design** week over week — same
cell line, same vehicle/DMSO control, same workflow — but with **different
compounds** swapped into the treatment arm each week. This means the control
arm is, by intent, the *same biological experiment* repeated weekly.

If the control arm drifts (instrument re-calibration, reagent batch change,
sample handling shift, cell-state drift), every cross-week compound comparison
silently inherits that drift. Today the package has no way to look at the
control samples across weeks and ask "are we still measuring the same thing?"

## Short answer to the question — yes, SPC

**Statistical Process Control** is the canonical framework for exactly this:
detect when a stable process leaves its expected operating range. Proteomics
groups in pharma already apply SPC to sample-level QC metrics; we should bring
it into this package as a first-class capability rather than a side analysis.

The relevant SPC primitives:

- **Shewhart / X-bar charts** — per-protein control intensity tracked across
  weeks with control limits at ±3σ from the process mean
- **EWMA / CUSUM** — more sensitive to small sustained drifts than Shewhart;
  useful for catching slow reagent or instrument decay
- **Western Electric / Nelson rules** — pattern-based detection (e.g. "8
  consecutive points on the same side of the mean" → drift, even within
  ±3σ)
- **Common-cause vs special-cause** — separates the noise budget you accept
  from the events you must act on

## What "control performance" means here

Two layers, both worth instrumenting:

1. **Sample-level metrics** across weeks (one chart per metric):
   - Median intensity of control samples
   - % missing values in control samples
   - Sample-to-sample correlation within the control replicates
   - Total identified proteins in control

2. **Protein-level metrics** across weeks:
   - Per-protein log-intensity in control, tracked over time
   - Useful for spotting drift in specific proteins (e.g. housekeeping
     proteins drifting suggests sample handling; specific pathway proteins
     drifting suggests cell-state shift)
   - Optionally flag proteins with the largest excursions for review

## Likely shape of the solution

1. **Data model addition.** A "control" sample needs to be identifiable as
   such across weeks. Either extend `Annotate Experiment` with a "control
   reference type" tag, or piggyback on the existing group-1 / group-2
   annotation if the convention is "group-1 is always control."

2. **Cross-run aggregation.** Build a longitudinal table indexed by
   `(week, control replicate, protein)` from a set of selected weekly
   analyses. Reuses the screening-campaign concept being captured in the
   adjacent compound-comparison todo.

3. **SPC viewer.** A new viewer that, given the longitudinal table and a
   metric:
   - Plots the metric across weeks with control limits drawn
   - Color-codes points by Western Electric rule violations
   - Lets the user click a flagged point to drill into the offending week's
     analysis
   - Has a phase-aware mode (e.g. "we changed reagent lot at week 5" — split
     the process mean before/after)

4. **Top-N drift surface.** A summary table listing the proteins or sample
   metrics with the largest unexplained drift, ranked by how many SPC rules
   they currently violate. This is what a QA reviewer wants to open first.

5. **Optional automation.** Tag a weekly run as `control-passed` /
   `control-flagged` so the cross-week compound comparison (adjacent todo)
   can refuse to combine runs whose controls were drifting.

## Open questions for investigation

- Process-mean baseline: rolling window vs fixed reference period? How many
  initial runs are needed before the chart is trustworthy?
- Which SPC rule set to default to — Shewhart only (simple) vs full Nelson
  (sensitive)? Should the user pick per metric?
- Per-protein SPC at proteome scale (thousands of proteins × many weeks) gets
  expensive and multiplies false positives via multiple testing. Likely need
  to either (a) restrict to a curated panel of QC proteins (housekeeping +
  pathway-specific anchors) or (b) apply an FDR correction to the SPC rule
  violations.
- Are control replicates within a single week treated as one point (averaged)
  or as multiple points (giving within-week dispersion as part of the chart)?
  Latter is more informative but harder to plot legibly.
- How does this interact with existing within-run QC dashboard? Same viewer,
  different mode? Separate menu entry? Likely separate — within-run is per
  sample, cross-run is per week.

## Scope warning

This is a **new analytical capability** with mature statistical underpinnings
(SPC is decades old) but no current foothold in the package. Scope-wise it's
smaller than the cross-team-sharing or compound-comparison use cases —
roughly one phase of work if scoped tightly to (1) longitudinal table builder
+ (2) one SPC viewer + (3) one default rule set. Could ship in the same
milestone as the screening-campaign data-model work since both depend on
introducing "screening run" as a first-class concept.

## Related

- `packages/Proteomics/.planning/todos/pending/2026-05-11-compare-most-impactful-compounds-across-weekly-screening-runs.md`
  — depends on the same screening-campaign data-model concept (week, run,
  control vs treatment annotation). This SPC use case can be thought of as
  the *quality gate* before the compound-comparison view is meaningful — if
  the controls are drifting, the compound delta is suspect.
- `packages/Proteomics/src/viewers/qc-dashboard.ts` — existing within-run QC
  dashboard. The cross-run SPC view is its longitudinal sibling.
