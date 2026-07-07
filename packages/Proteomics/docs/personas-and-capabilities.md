# Personas & capabilities

Who does what with a Proteomics dataset, what each persona is *allowed* to do, and
**how that boundary is enforced** — today and where we are taking it.

This doc is the reference behind the README's ["What you can achieve"](../README.md#what-you-can-achieve)
section. The README leads with outcomes; this leads with the capability boundary and its
enforcement, because a client requirement (see [Target state](#target-state)) needs that
boundary to be real, not conventional.

## The two personas

| Persona | Role | Outcome they want |
|---|---|---|
| **Proteomics analyst** (data generator) | Builds the analysis: import → annotate → normalize → impute → differential expression | A defensible, reproducible, FDR-controlled result they can stand behind and hand off |
| **Biology scientist** (data consumer / reviewer) | Interprets a finished analysis: volcano, enrichment, per-protein annotation | To understand *what changed and what it means* — without touching (or being able to change) the analysis itself |

The dividing line is **"the interpretation"** — the set of decisions that determine what the
result *says*: group assignment, normalization, imputation, the DE method, and the
significance thresholds. The analyst owns these. The biologist reads the consequences and
explores them, but must not alter them.

## Capability matrix

Legend: ✅ allowed · 🚫 must not · 👁 read-only interaction (view/filter/select, no state change)

| Capability | Analyst | Biologist | Changes the interpretation? |
|---|---|---|---|
| Import a search-engine file | ✅ | — (receives a shared snapshot instead) | Produces it |
| Annotate Experiment (assign groups) | ✅ | 🚫 | **Yes** |
| Set Log2 Scale | ✅ | 🚫 | **Yes** |
| Normalize | ✅ | 🚫 | **Yes** |
| Impute Missing Values | ✅ | 🚫 | **Yes** |
| Differential Expression | ✅ | 🚫 | **Yes** |
| Compute SPC Status | ✅ | 🚫 | **Yes** |
| Volcano Plot | ✅ | 👁 | No (view) |
| Heatmap / PCA / Group-Mean Correlation / QC Dashboard | ✅ | 👁 | No (view) |
| Enrichment Analysis & Charts | ✅ | ✅ (derives biology; non-destructive) | No |
| UniProt protein panel | ✅ | ✅ | No |
| Filter / select / cross-link in any viewer | ✅ | ✅ | No |
| Share Analysis for Review | ✅ | 🚫 | N/A (publishes) |
| Save changes to a shared snapshot | ✅ (own copy) | 🚫 | **Yes** |

The 🚫 rows are the ones CK cares about: a biologist must not be able to re-run the pipeline
or re-thresholds and thereby produce a *different* interpretation of the same data.

## How the boundary is enforced today

The package differentiates on **three axes** — none of which is the user's identity.

### 1. By data shape — the only code-enforced line (`sampleOnly` grey-out)

`src/menu.ts` gates every analysis-mutating action behind
`sampleLevelDisabledReason(df)`. On a Spectronaut **Candidates** file (pre-computed DE, no
per-sample intensities) these items grey out with an explanatory tooltip:

> *Not applicable to a Spectronaut Candidates file — it already carries computed
> differential expression and has no per-sample intensities. Use Volcano Plot, Enrichment
> Analysis or the UniProt panel, or import the matching Spectronaut Report for sample-level
> analyses.*

Gated (sample-level): Annotate Experiment, Set Log2 Scale, Normalize, Impute, Differential
Expression, Compute SPC Status, Heatmap, PCA, Group-Mean Correlation, QC Dashboard, Show All
Visualizations.

Always available: Import, Volcano Plot, Enrichment Analysis + Charts, Export Enrichment
Inputs, SPC Dashboard, Share for Review, UniProt panel.

**This gates on the *data*, not the *person*.** It happens to keep a biologist who only ever
receives Candidates files in the interpretation lane — but it does nothing to stop a
biologist who opens a **Report** file from running the whole pipeline.

### 2. By workflow stage — the menu grouping mirrors the roles

`Import → Annotate → Analyze` is the analyst building the result; `Visualize → Share` is the
consumer interpreting it. The volcano auto-opening on both DE-completion and Candidates
import lands the interpreter on the deliverable without making them run the pipeline. This
is a *convention*, not a restriction.

### 3. By access — view-only sharing (the one real identity boundary)

`Share Analysis for Review` publishes a frozen, versioned snapshot and grants the reviewer
group **View only**. The publish flow's verify-and-rollback gate aborts if the reviewer
group has any Edit / Share / Delete right at the project, child-space, or umbrella-space
level. So through the **shared snapshot**, a biologist genuinely cannot change the
interpretation — they get a read-only copy.

## The gap

Enforcement today is **data-shape + access**, not **identity**. Concretely:

- A biology user who opens a **Report** file (not a Candidates file, not a shared snapshot)
  sees the full analyst menu enabled and can run/redo DE, re-normalize, re-threshold, and
  re-share. Nothing checks who they are.
- The view-only guarantee exists **only** for analyses that arrived via Share for Review. A
  biologist working on a live table is unconstrained.

CK does not want their biologists to be able to change the interpretation *at all* — which
data-shape gating cannot guarantee, because it depends on which file they happened to open.

## Target state

Add an **identity/role** axis on top of the existing data-shape and access gating, so the
🚫 capabilities are unavailable to a biologist regardless of file type.

Direction (to be specified in an implementation phase):

1. **A "Proteomics Analysts" Datagrok group.** Analyst-only actions require membership.
2. **Menu enforcement.** Extend the `sampleOnly` predicate in `src/menu.ts` to an
   `analystOnly` predicate that also greys out (with a clear reason) when the current user
   is not in the analysts group — so the interpretation-changing items disappear for
   biologists on *any* table, not just Candidates.
3. **Handler enforcement.** Guard the mutating handlers in `src/package.ts` (the same place
   `requireSampleLevelData` already lives) so the capability is enforced even if a function
   is invoked directly, not only via the menu — the menu grey-out is UX, the handler guard
   is the security boundary.
4. **Reviewer mode is the default for biologists.** A biologist's normal entry point stays
   the view-only shared snapshot; the role gate is the backstop for the live-table case.

Open questions for the implementation phase: whether Enrichment re-runs count as "changing
the interpretation" (currently treated as non-destructive exploration), how the analysts
group is provisioned per client/deployment, and whether to expose an explicit read-only
"Reviewer mode" toggle vs. inferring it from group membership.

## See also

- [`README.md` → What you can achieve](../README.md#what-you-can-achieve) — the outcome-led persona intro
- [`CLAUDE.md` → Menu architecture](../CLAUDE.md) — how the ribbon menu and `isEnabled` grey-out are built
- `src/menu.ts` — `sampleLevelDisabledReason`, the `sampleOnly` gate
- `src/package.ts` — `requireSampleLevelData`, the handler-level guards
- `src/publishing/` — the Share for Review view-only publish flow
