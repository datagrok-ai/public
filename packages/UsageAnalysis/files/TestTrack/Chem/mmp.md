---
feature: chem
sub_features_covered: [chem.analyze.mmp, chem.analyze.mmp.top-menu, chem.analyze.mmp.editor, chem.analyze.mmp.viewer, chem.demos.mmpa]
target_layer: playwright
coverage_type: edge
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/mmp.md
migration_date: 2026-05-11
migration_report: mmp-migration-report.md
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

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: simple`, `pyramid_layer: bug-focused`, `target_layer: playwright`, strategy
`simple`. Atlas critical path `chem.cp.mmp-analysis` (p1) and the GROK-18517 entry both name
`mmp_demo` as the "golden run" for MMP. UI coverage owned
(`ui_coverage_delegated_to: null`) over the MMP editor select-activities flow and all four
viewer tabs.

## Setup

1. **Provision linked dataset.** `System:DemoFiles/chem/mmp_demo.csv` (per source JSON footer
   `order: 14`, `datasets: ["System:DemoFiles/chem/mmp_demo.csv"]`). The CSV is bundled in
   the Chem package's demo file share — no external provisioning required. The dataset has
   four columns: `smiles` (molecules), `CMPD_CHEMBLID` (ChEMBL identifiers), `CYP3A4`
   (numeric activity), and `hERG_pIC50` (numeric activity). Two activity columns total — the
   `MMPEditor` dialog will surface both for selection (atlas `chem.analyze.mmp.editor` line
   421).
2. **Confirm Chem package is loaded** so that the **Chem > Analyze > Matched Molecular
   Pairs...** top-menu entry is registered (`mmpAnalysis`, atlas `chem.analyze.mmp.top-menu`
   line 406, package source `Chem/src/package.ts#L2273`) and the `MMPEditor` dialog +
   `Matched Molecular Pairs Analysis` viewer (atlas `chem.analyze.mmp.viewer` line 414,
   role `viewer`, hidden from gallery) surface on demand.
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
   (`mmpAnalysis`, atlas `chem.analyze.mmp.top-menu` line 406). The `MMPEditor` custom
   dialog opens (atlas `chem.analyze.mmp.editor` line 421 — "custom dialog that gathers
   molecules column, activity columns, diff types, scaling methods, fragment cutoff").
   Verify the dialog renders without console errors.
3. In the `MMPEditor` dialog, select **both** activity columns — `CYP3A4` and `hERG_pIC50`
   — in the activity-columns selector. Confirm the molecules column auto-binds to `smiles`
   and the diff-types / scaling-methods / fragment-cutoff fields show sensible defaults.
   Click **OK** to run analysis.

   **GROK-18517 discriminator (regression guard):** the MMP analysis must complete without
   firing a minified JS runtime error (`J.aS(...).b7 is not a function` or equivalent
   tree-shook-function signature). On success, the `Matched Molecular Pairs Analysis`
   viewer renders (atlas `chem.analyze.mmp.viewer` line 414); on failure (bug reproduction
   condition), an error balloon surfaces with the minified-function-name message and the
   viewer never renders. **The scenario fails the moment a minified runtime error fires
   during MMP generation** — this is the primary correctness invariant of this bug-focused
   migration.
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

- **`coverage_type: edge`** — bug-focused migration anchored on GROK-18517 ("Chem: MMP:
  failed to generate analysis for mmp_demo dataset"). The scenario IS the bug's
  reproduction path verbatim: open `DemoFiles/chem/mmp_demo` → Chem > Analyze > Matched
  Molecular Pairs → select both activities → click OK → exercise viewer tabs. Per chain
  YAML `output_plan.mmp.md.reason`: "Bug-focused scenario (GROK-18517 mmp_demo
  reproduction) layered on canonical MMP viewer walk (4 tabs)." Per chain YAML
  `mmp.md.pyramid_layer = bug-focused` rationale: "scenario body IS the GROK-18517
  reproduction path verbatim ... Discriminator test: GROK-18517 WOULD fail this scenario
  before fix". `coverage_type: edge` is the natural fit for a bug-focused regression
  guard against a specific failure mode (minified runtime error tree-shook out of the
  MMP code path on the bundled demo dataset).
- **GROK-18517 cross-references.** Bug entry per
  `bug-library/chem.yaml#GROK-18517` — priority `p2`, status `fixed`, `affects`:
  `[chem.analyze.mmp, chem.analyze.mmp.top-menu, chem.analyze.mmp.editor,
  chem.demos.mmpa]` — all four sub-features are covered by this scenario.
  `edge_case_for_atlas` invariant: "every Chem demo dataset (mmp_demo, scaffold_demo,
  similarity_demo, etc.) successfully runs its corresponding analysis workflow end-to-
  end — regression coverage for demo-vs-current-code compatibility." This scenario
  satisfies the `mmp_demo` arm of that invariant; sister `scaffold_demo` / `similarity_demo`
  coverage is owned by other Chem demo scenarios.
- **Atlas critical path.** `chem.cp.mmp-analysis` (p1) sub-features used:
  `chem.analyze.mmp`, `chem.analyze.mmp.top-menu`, `chem.analyze.mmp.viewer`,
  `chem.analyze.mmp.editor`, `chem.demos.mmpa`. Atlas description: "OK runs analysis →
  MMP viewer renders fragments / generations / cliffs tabs without minified-runtime
  errors (GROK-18517 regression surface). Demo dataset `mmp_demo` is the golden run."
  This scenario IS the critical path walk.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`chem-add-mmp`, `chem-mmp-editor-select-activities`,
  `chem-mmp-viewer-transformation-tab`, `chem-mmp-viewer-fragments-tab`,
  `chem-mmp-viewer-cliffs-tab`, `chem-mmp-viewer-generation-tab`) is exercised via UI
  driving — top-menu walk, MMPEditor dialog interaction (activity column selector + OK
  button), and four tab switches with content verification. Direct invocation of
  `mmpAnalysis` via JS API does not exercise the `MMPEditor` dialog gate or the
  minified-runtime-error reproduction surface (the bug surfaced during MMP generation
  with the editor-driven flow); UI driving is mandatory for this scenario.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/mmp-spec.ts` (per
  `existing-test-index.yaml`). The Automator will extend it to cover the full 7-step
  walk per this migrated scenario (or write fresh if existing coverage is partial).
- **First 7 entries pattern.** Original step "click first 7 entries" preserved verbatim
  as a specific Transformation-tab walk depth. The number is deliberate per the original
  scenario authoring — covering enough breadth in the transformation grid to surface any
  per-row regression while staying bounded for spec runtime. Automator at spec time
  drives 7 row clicks via the MMP viewer's transformation grid surface.
- **Bundled demo dataset criticality.** Per GROK-18517 expected: "Demo workflows must
  always work because they are the primary onboarding path." The bundled `mmp_demo.csv`
  is the primary onboarding path for MMP — a regression here is high-visibility.
- **Order in chain.** `order: 14` per source JSON footer — last scenario in the Chem
  section.
- **Source-text fixes applied.** Original steps 4 and 5 reference "all three activities"
  but the bundled `mmp_demo.csv` has only TWO activity columns (`CYP3A4` and
  `hERG_pIC50` — verified by reading the dataset header). The "three" is an authoring
  defect; original step 3 correctly references "two activities". Migrated body resolves
  by silently fixing "three" → "both" / "two" per the actual dataset shape (per chain
  rev 2 directive footer note (c) Olena 2026-05-11 — silent canonical-phrasing fixes
  permitted with full disclosure in migration report § Source-text fixes).
