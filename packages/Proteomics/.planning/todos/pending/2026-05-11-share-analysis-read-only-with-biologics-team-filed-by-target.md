---
created: 2026-05-11T14:00:00.000Z
title: Share analysis read-only with biologics team, filed by target
area: sharing
promoted_to: Phase 999.3 (ROADMAP.md backlog)
promoted_on: 2026-06-03
resolves_phase: 15
resolves_requirements: PUB-01, PUB-02, PUB-03, PUB-04, PUB-05, PUB-06, PUB-07, PUB-08, PUB-09, PUB-10, PUB-11, PUB-12, PUB-13
files:
  - packages/Proteomics/src/analysis/differential-expression.ts
  - packages/Proteomics/src/viewers/volcano.ts
---

> **Promoted to backlog phase 999.3 on 2026-06-03** with paired user-stories at the ROADMAP level:
> - *As a proteomics expert, I should be able to publish a read-only analysis that biologists will interrogate.*
> - *As a biologist, I should be able to interrogate a proteomics experiment by examining the volcano plot and enrichments without concern that I will change its content.*
>
> This todo retains the detailed solution shape, Datagrok-platform integration notes, and open questions. Phase 999.3 is the ROADMAP-visible reference; this file is its scoping appendix.

## Problem — new use case

Once the proteomics team has produced a DE analysis, they want to **hand it
off to the biologics team in read-only mode**, filed under a specific drug
target. The biologics reviewer should be able to open the volcano plot, browse
significant hits, and (eventually) leave comments — but should *not* be able to
re-run normalization, change DE parameters, edit columns, or download the raw
intensity matrix.

Two constraints make this non-trivial:

1. **Minimum data shared** — biologics needs only what the volcano review
   requires (protein IDs, gene symbols, log2FC, p-value, adjusted p-value,
   significance flag — possibly accompanying enrichment results). The raw
   per-sample intensities, peptide counts, and intermediate normalized
   columns should *not* travel with the share.

2. **Target-based filing** — the artifact lives under "target X" in a
   biologics-facing organization, not under the proteomics author's personal
   workspace. A biologist navigating to "target X" should see all proteomics
   analyses associated with that target without knowing which proteomics
   scientist produced each one.

Cytokinetics meeting on 2026-05-04 surfaced this; expect similar asks from
any pharma customer with a target-driven discovery org.

## Likely shape of the solution

Datagrok already has the primitives for this — **Spaces** for grouping,
permissions for read-only access, and projects for serializing a view. The
work is in composing them and exposing a one-click "share for review" flow
that does the right thing.

1. **Choose the container.** Most likely: a shared **Space** per target (or a
   top-level "Targets" Space with **subspaces** per target). Biologics group
   gets View permission on the Space; only proteomics group can publish into
   it.

2. **Define the minimum shareable artifact.** Probably a Datagrok project
   containing:
   - A trimmed DataFrame (protein ID, gene symbol, log2FC, p, adj.p, sig flag,
     optionally enrichment terms)
   - A pre-built volcano viewer layout
   - The dashboard title, group names from `Annotate Experiment`, and the DE
     method used (audit context: "limma 2-tailed, FC≥1.5, adj.p<0.05")
   - No raw intensities, no peptide counts, no original file
   Decide whether to embed the trimmed DataFrame in the project or store it
   as a separate read-only table the project references.

3. **Add a "Share for Review" entry to the package menu.** Probably
   **Proteomics → Share Analysis for Review…** that prompts for:
   - Target (autocomplete against existing target Spaces, or create new)
   - Reviewer group (default: biologics)
   - Free-form note / context for the reviewers
   Then builds the trimmed project, files it under the target subspace,
   and applies View-only permissions to the reviewer group.

4. **Reviewer-side affordances.** When biologics opens the shared project:
   - Volcano renders immediately with significance thresholds applied
   - A "details" panel shows the audit context (method, thresholds, groups,
     source proteomics author, share date)
   - Comments thread / annotations enabled, edits disabled
   - "Request re-run with different parameters" button that pings the
     proteomics author rather than letting biologics edit the analysis itself

## Open questions for investigation

- Does Datagrok's project serialization let us strip columns *before*
  publishing, or do we need to clone the DataFrame first? (Probably the
  latter — cleaner semantics.)
- How are Spaces/subspaces created programmatically? `grok.dapi.spaces`
  vs `grok s spaces ...` CLI. Verify both paths.
- Permissions model: per-Space ACL vs per-Entity ACL. Which composes better
  with a target taxonomy that may grow to hundreds of entries?
- "Filed by target" — does target need to be a first-class entity in this
  package (semantic type, lookup) or is it sufficient to use the Space name?
- How does the biologics reviewer surface "I'm done reviewing, here are my
  notes"? Comments on the project? A separate workflow object?

## Scope warning

This is a **new use case** spanning UI (sharing dialog), data engineering
(trimmed DataFrame contract), and platform integration (Spaces, permissions,
projects). Likely a phase of its own — possibly the spine of a new v1.5
milestone "Cross-team review workflow." Worth surfacing to `/gsd-review-backlog`
when the next milestone is being scoped.

## Related

- `packages/Proteomics/.planning/todos/pending/2026-03-03-expand-import-to-handle-local-drive-datagrok-folders-and-existing-dataframes.md`
  — adjacent to this in spirit (file/share-system integration), but inbound
  rather than outbound
- Datagrok platform docs: Spaces, Permissions, Projects (no in-repo
  reference; see datagrok.ai help)
