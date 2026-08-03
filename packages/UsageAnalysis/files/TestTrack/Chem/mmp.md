---
feature: chem
target_layer: playwright
coverage_type: edge
priority: p0
realizes_atlas: [chem.cp.mmp-analysis]
realizes: [chem.analyze.matched-molecular-pairs]
realized_as:
  - mmp-spec.ts
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/mmp.md
migration_date: 2026-05-11
source_text_fixes:
  - three-activities-both-activities-two-activities
candidate_helpers:
  - helpers.playwright.chem.openMmpEditor
  - helpers.playwright.chem.selectMmpActivityColumns
  - helpers.playwright.chem.switchMmpViewerTab
unresolved_ambiguities:
  - two-vs-three-activities-inconsistency-in-source-body
  - first-7-entries-in-transformation-tab-what-is-an-entry
  - diff-types-scaling-methods-fragment-cutoff-defaults
scope_reductions: []
related_bugs: [GROK-18517]
---

# Chem | Matched Molecular Pairs on mmp_demo (GROK-18517 regression surface)

End-to-end MMP analysis walk on the bundled `mmp_demo.csv` demo dataset: open the dataset, run
`Chem | Analyze | Matched Molecular Pairs...`, select both activity columns (`CYP3A4` and
`hERG_pIC50`) in the `MMPEditor` dialog, then exercise all four MMP viewer tabs
(Transformation, Fragments, Cliffs, Generation) and verify no minified JS runtime error
fires. This scenario IS the GROK-18517 reproduction path verbatim — the bug was
"`J.aS(...).b7 is not a function`" (minified runtime error tree-shook out of the MMP code path
on the bundled demo dataset). The discriminator invariant: **MMP generation succeeds without
a minified-runtime error and all four viewer tabs render**.

`mmp_demo.csv` is considered the canonical dataset for exercising MMP analysis end-to-end.

## Setup

1. **Provision linked dataset.** `System:DemoFiles/chem/mmp_demo.csv` (per source JSON footer
   `order: 14`, `datasets: ["System:DemoFiles/chem/mmp_demo.csv"]`). The CSV is bundled in
   the Chem package's demo file share — no external provisioning required. The dataset has
   four columns: `smiles` (molecules), `CMPD_CHEMBLID` (ChEMBL identifiers), `CYP3A4`
   (numeric activity), and `hERG_pIC50` (numeric activity). Two activity columns total — the
   `MMPEditor` dialog will surface both for selection.
2. **Confirm Chem package is loaded** so that the **Chem > Analyze > Matched Molecular
   Pairs...** top-menu entry is registered (`mmpAnalysis`; package source
   `Chem/src/package.ts#L2273`) and the `MMPEditor` dialog + `Matched Molecular Pairs Analysis`
   viewer (role `viewer`, hidden from gallery) surface on demand.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. The dataset is opened fresh
   inside the scenario and the MMP analysis is in-session only (`produces` field of chain
   entry: "MMP analysis on mmp_demo.csv with Transformation/Fragments/Cliffs/Generation
   tabs exercised across both activities (in-session only)"); the scenario does not persist
   the MMP viewer or analysis result.

## Scenarios

### Run MMP on mmp_demo with both activities and walk all four viewer tabs

Single linear walk (6 steps mirroring the original 6 numbered bullets, with expected-result
verifications woven inline per D-STEP-02; the source-text defect "three activities" in
original steps 4 and 5 silently fixed to "both activities" — see migration report § Source-
text fixes):

1. Open `System:DemoFiles/chem/mmp_demo.csv` (close any previously open views first to
   start from a clean state). Wait for the table view to render with the `smiles` column
   populated and the RDKit cell renderer applied. Verify the table view has 4 columns —
   `smiles`, `CMPD_CHEMBLID`, `CYP3A4`, `hERG_pIC50` — and that the two numeric activity
   columns (`CYP3A4`, `hERG_pIC50`) are detected as numeric.
2. From the top menu, run **Chem > Analyze > Matched Molecular Pairs...**
   (`mmpAnalysis`). The `MMPEditor` custom dialog opens (custom dialog that gathers
   molecules column, activity columns, diff types, scaling methods, fragment cutoff).
   Verify the dialog renders without console errors.
3. In the `MMPEditor` dialog, select **both** activity columns — `CYP3A4` and `hERG_pIC50`
   — in the activity-columns selector. Confirm the molecules column auto-binds to `smiles`
   and the diff-types / scaling-methods / fragment-cutoff fields show sensible defaults.
   Click **OK** to run analysis.

   **GROK-18517 discriminator (regression guard):** the MMP analysis must complete without
   firing a minified JS runtime error (`J.aS(...).b7 is not a function` or equivalent
   tree-shook-function signature). On success, the `Matched Molecular Pairs Analysis`
   viewer renders; on failure (bug reproduction condition), an error balloon surfaces with
   the minified-function-name message and the viewer never renders. **The scenario fails
   the moment a minified runtime error fires during MMP generation** — this is the primary
   correctness invariant of this bug-focused migration.
4. With the MMP viewer rendered, click the first 7 entries (rows) in the Transformation
   tab and observe per-row changes in the tab's molecule rendering / activity-difference
   columns. Verify each click updates the visible transformation pair (parent → product
   with the substitution highlighted) and the activity-difference values across both
   `CYP3A4` and `hERG_pIC50` are surfaced. No console errors per click.
5. Switch to the **Fragments** tab and verify both activity columns (`CYP3A4`,
   `hERG_pIC50`) are shown — the fragments tab summarizes per-fragment activity
   differences across both numerical activity columns selected at editor time. Verify the
   fragments grid renders without console errors.
6. Switch to the **Cliffs** tab and verify the filter controls for both activity columns
   (`CYP3A4`, `hERG_pIC50`) are present and functional — adjusting a filter (e.g. activity
   range slider) should affect the cliffs scatter visualization. Verify the cliffs view
   renders the activity-difference scatter without console errors.
7. Switch to the **Generation** tab and verify the generation results render — the
   generation view shows candidate molecule variants produced by applying observed MMP
   transformations. Verify the generation grid populates without console errors.

Implicit cross-step invariants:

- **GROK-18517 invariant:** no minified-runtime error fires at any step during MMP
  generation or tab navigation. The bug's reproduction surface is steps 2 + 3 (open
  editor + click OK with both activities selected); the four-tab walk in steps 4-7
  exercises the post-generation surfaces and would also fail if MMP generation half-
  completed.
- **Both activities propagate.** Both `CYP3A4` and `hERG_pIC50` selected at editor time
  must appear in the Fragments tab summary and the Cliffs tab filter controls — one
  activity column dropping out of either tab is a coverage gap regression.
- **No console errors throughout.** Each tab switch and each transformation-entry click
  must complete without console errors. The minified runtime error of GROK-18517
  surfaces in the JS console — console-error introspection is the most sensitive
  detector.

## Notes

- **GROK-18517 reference.** *Chem: MMP: failed to generate analysis for mmp_demo dataset* (priority
  p2, fixed). The underlying expectation is broader than this one dataset: every bundled Chem demo
  dataset (mmp_demo, scaffold_demo, similarity_demo, etc.) should successfully run its
  corresponding analysis workflow end-to-end. This scenario covers the `mmp_demo` case; the
  `scaffold_demo` / `similarity_demo` cases are covered by other scenarios.
- **Why 7 transformation entries.** Clicking the first 7 rows in the Transformation tab is a
  deliberate breadth choice — enough to surface a per-row regression while keeping the scenario's
  runtime bounded.
- **Why this matters.** `mmp_demo.csv` is the primary onboarding path for the MMP feature — a
  regression here is high-visibility, since demo workflows are usually the first thing a new user
  or evaluator tries.
- **Sibling spec.** A Playwright spec already exists at `mmp-spec.ts`; it is extended (or rewritten,
  if existing coverage is too partial) to cover the full 7-step walk described here.
