---
created: 2026-05-04T13:55:00.000Z
title: Add explanatory help to all analytics dialogs and viewers
area: ui
files:
  - packages/Proteomics/src/analysis/normalization.ts
  - packages/Proteomics/src/analysis/imputation.ts
  - packages/Proteomics/src/analysis/differential-expression.ts
  - packages/Proteomics/src/analysis/enrichment.ts
  - packages/Proteomics/src/analysis/experiment-setup.ts
  - packages/Proteomics/src/analysis/pca.ts
  - packages/Proteomics/src/viewers/qc-dashboard.ts
  - packages/Proteomics/src/viewers/volcano.ts
  - packages/Proteomics/src/viewers/heatmap.ts
---

## Problem

When a user opens an analytics dialog (Normalize, Impute Missing Values,
Differential Expression, Enrichment Analysis, Annotate Experiment) or sees a
result viewer (volcano, PCA, heatmap, QC panels), nothing in the UI explains
what the function is trying to do or what to look for in the output.

Some field-level tooltips already exist —
`differential-expression.ts` (method picker, FC threshold, p-value), `imputation.ts`
(downshift, width), `enrichment.ts` (FC, p-value) — so the per-input layer is
partially done. Missing layer is the **dialog-level / viewer-level** "what
question does this answer" explanation, which is what a new user (e.g. the
Cytokinetics audience the 2026-05-04 demo is being prepared for) actually needs
before they touch a single input.

A scientist seeing the QC dashboard for the first time should be able to learn,
without leaving the platform, why CV plots and MA trends matter for label-free
proteomics. A reviewer seeing a volcano plot should be told what the axes mean
and what "significance" means in this package's hands.

## Solution

Investigate, but the likely shape:

1. **Standardize a help affordance** — small (?) icon next to each dialog
   title, and one in each viewer's header bar. Clicking opens a popover (or a
   `grok.shell.info(...)` panel) with:
   - 2–3 sentences on what the function does
   - What inputs the user is choosing between
   - What to look for in the result (e.g. "spike-ins should appear in the upper
     corners of the volcano")
   - A link to deeper docs when they exist

2. **Source of truth for the text** — the function's `//description:` metadata
   is already a one-liner; extend with a longer-form `meta.help` or a co-located
   `*.help.md` file the dialog reads at open time. Co-located markdown is
   easier to edit and review than escaped TS strings.

3. **Apply consistently** — at minimum every Proteomics menu item that
   performs analytics or shows a result. Build it as a small helper
   (`buildHelpButton(text)`) so adding one to a new dialog/viewer is a single
   line.

4. **Acceptance:** every Proteomics menu item has a one-click route to a
   plain-language explanation of what it computes and how to read the result,
   without leaving the platform.

## Scope warning

This touches every analytics dialog and viewer in the package — likely a phase
of its own rather than a single ad-hoc commit. Consider promoting to the
backlog when the next milestone is being scoped.

## Related

- `packages/Proteomics/.planning/todos/pending/2026-03-10-improve-chart-titles-and-axis-labels-across-all-viewers.md`
  — the labels-and-titles cousin of this todo; both improve viewer
  intelligibility and would naturally batch into one milestone
- `packages/Proteomics/.planning/todos/pending/2026-05-04-reorganize-qc-dashboard-layout-for-navigability.md`
  — the layout reorganization shares the same "help users navigate the
  dashboard" goal
