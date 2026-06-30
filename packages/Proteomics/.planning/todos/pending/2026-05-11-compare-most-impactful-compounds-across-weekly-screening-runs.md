---
created: 2026-05-11T14:05:00.000Z
title: Compare most-impactful compounds across weekly screening runs
area: analysis
resolves_phase: 17, 18
resolves_requirements: CAMP-01, CAMP-02, CAMP-03, CAMP-04, CAMP-05, CAMP-06, CAMP-07, CAMP-08, CAMP-09, CAMP-10, CAMP-11, CAMP-12, CAMP-13, CAMP-14, CAMP-15, CAMP-16
files:
  - packages/Proteomics/src/analysis/differential-expression.ts
  - packages/Proteomics/src/analysis/experiment-setup.ts
  - packages/Proteomics/src/viewers/volcano.ts
---

## Problem — new use case

A pharma proteomics group runs a **weekly screening campaign**: each week N
compounds are tested in a single batch, each producing its own DE analysis
(compound vs vehicle/DMSO). At the end of each week, the team identifies the
**most impactful compound** — the one whose treatment produced the largest
proteome-level response. They then want chemists to **compare this week's
top compound vs last week's top compound** side by side, to understand whether
the chemistry is converging on a useful effect or drifting.

Today this comparison is done by exporting DE results from two separate
analyses into Excel, hand-aligning on protein/gene, and eyeballing differences.
Chemists never see the compounds in their native chemistry context (no
structure side-by-side, no SAR view). The package currently has no concept of
"compound" or "screening run" — every analysis is one-off and standalone.

Cytokinetics 2026-05-04 meeting surfaced this; expect similar asks from any
pharma running phenotypic or proximity-based proteomic screens.

## What "most impactful" needs to mean

Open question worth nailing down before any coding. Candidates:

- **Count** — number of proteins with adj.p < α and |log2FC| > τ (simple,
  threshold-dependent)
- **Magnitude** — sum or mean of |log2FC| across significant hits
- **Top-N strength** — average |log2FC| of the top-N most significant proteins
- **Pathway-level** — count of enriched GO/KEGG terms at FDR < α
- **User-selectable** — let the screening lead pick the rule per campaign

Whichever rule is chosen, it has to be **reproducible** — i.e. baked into the
DataFrame or project as a tag, not recomputed ad hoc.

## Likely shape of the solution

1. **Lift the data model from "one analysis" to "screening campaign."**
   Introduce first-class concepts for:
   - **Screening run** — a weekly batch (date, list of compounds tested, the
     vehicle/DMSO reference)
   - **Compound** — name + structure (SMILES or molfile); ties into the
     existing Datagrok Chem package's molecule semantic type and cell renderer
   - **Per-compound DE** — the DE columns the package already produces, tagged
     with the originating compound + week

2. **Add a "Top compound this week" scoring step.** After all per-compound DEs
   in a screening run complete, compute the impact score (rule above) and
   write it back as a tag on the analysis. Surface a ranked table:
   *"Week of 2026-05-11: top-impact compounds"*.

3. **Build a cross-run compound comparison view.** Pick week1 + week2 (or any
   two analyses tagged with compound identity). Render:
   - Two volcano plots side by side, linked selection, shared thresholds
   - A diff table: proteins significant in A only, B only, both (with
     log2FC_A, log2FC_B, delta)
   - The two compound structures rendered at the top via the Chem cell
     renderer (chemist-oriented framing)
   - Optional: enrichment overlap (terms common to both, unique to each)

4. **Chemist-oriented affordances.** A chemist scrolling the diff table
   should be able to click a protein and see *which compound* moved it more.
   A chemist looking at the two structures should be able to see *what
   substructure differs* (Chem package's MCS / substructure highlight).

## Open questions for investigation

- Does the existing `Annotate Experiment` flow accommodate a *compound*
  axis, or do we need a new annotation step ("Annotate Screening Campaign")?
- How do compounds get into the platform — uploaded as a separate SDF/CSV per
  campaign, or queried from an internal compound DB connector?
- Is the comparison view a new dialog, a new top-level menu entry, or a
  context-action on a "screening campaign" entity?
- How does the chemist navigate from week1's top compound to week2's top
  compound — manual selection, a "compare to previous week" shortcut, or
  always-on side-by-side once two campaigns exist?
- Storage / sharing: is the cross-week comparison itself a saveable artifact
  the chemist can come back to? Ties into the
  `2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target`
  todo — chemists are a parallel consumer team.

## Scope warning

This is a **new use case** that introduces multi-experiment concepts the
package doesn't have today (screening campaign, compound entity, cross-run
comparison). Likely a phase of its own — possibly part of a "Screening &
SAR" milestone that pairs with the cross-team-sharing use case to form a
discovery-workflow story. Worth surfacing to `/gsd-review-backlog` when the
next milestone is being scoped.

## Related

- `packages/Proteomics/.planning/todos/pending/2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md`
  — companion downstream-consumer use case (biologics review vs chemist
  review). Both want curated read-only views of DE results filed under a
  workflow-relevant identifier (target / compound), so they may share
  primitives (trimmed DataFrame contract, "Share for Review" affordance)
- Datagrok **Chem** package — required for compound structure rendering,
  SMILES handling, MCS / substructure highlights
