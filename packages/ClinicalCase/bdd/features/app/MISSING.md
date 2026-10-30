# Clinical Case — what the features do not claim yet

## ClinicalCase

### Cannot be put back (ClinicalCase)

- **Opening a study** (a study link on the start page, or a study node in the tree) and every view under it:
  Summary, Timelines, Patient Profile, Adverse Events, Laboratory, AE Risk Assessment, Survival Analysis,
  Distributions, Correlations, Time Profile, Tree map, Medical History, Visits, Validation, and the Domains
  group. `initClinicalStudyData` (`src/utils/app-utils.ts`) writes `<study>.d42` into
  `System:AppData/ClinicalCase/SDTM/<study>/` (`readClinicalData`, ~line 618), and `ClinicalStudy.validate`
  (`src/clinical-study.ts:211`) runs `ClinicalCase:RunCoreValidate` — the CDISC rules engine in the
  package's container (`dockerfiles/`) — and writes `validation_results.json` there. A container run never
  becomes a feature (bdd CLAUDE.md), and the library has no step that deletes a server file without writing
  it first: the only one, `a file {string} with text {string} is on the server`, writes text, and the app
  would then read that text as a `.d42` or as validation results. Survival Analysis also runs
  `scripts/survivalPlot.r`. What would make it a feature: a library step that remembers a server folder and
  deletes what was added to it at feature end (`the files added to the folder {string} are deleted at
  feature end`), plus a stand that runs the ClinicalCase container behind `the stand runs the {string}
  service`.
- **Import** in the Import study view: writes the chosen files and `studyConfig.json` into
  `System:AppData/ClinicalCase/<standard>/<study>/` and then opens the study (above).
- **Run ClinicalCase** in the SDTM folder viewer: opens the domain tables and the app, whose listing writes
  the study configs (the app feature puts those back only because it knows when the listing happens).

### Needs a step (ClinicalCase)

- The Import study view with a valid study: the files input takes several files (`study.json` + `dm.csv`)
  at once, and `user uploads {string} through {element}` sets one. A step that sets a list of project
  files on one chooser would let this run without Import:

```gherkin
    When user uploads "../files/SDTM/CDISCPILOT01/study.json, ../files/SDTM/CDISCPILOT01/dm.csv" through "Open files" icon
    Then the current view should show the text "Summary"
    And the current view should show the text "CDISCPILOT01"
    And "Import" button should be enabled
```

  Signal: the Summary table `ui.tableFromMap` appends under the files (`createImportStudyView`).

### Findings 30 Oct 2026 (ClinicalCase)

- **Import study... on the start page closes the app and opens nothing when the Clinical Case tree node has
  not been expanded.** Hand steps (local stand, fresh page): Browse > Apps, click **Clinical Case** (the
  label, not the arrow) — the start page opens; click **IMPORT STUDY...**. Measured: the Clinical Case view
  closes and Home becomes current; no Import view (twice, runs p5c/p5d). Expected: the "Import Study -
  Clinical Case" view, as it opens when the node was expanded first (the app feature's third scenario
  passes that way). Source: `openApp`, `src/utils/app-utils.ts:191-195` — the button closes
  `currentOpenedView` (the start page itself) and selects the tree item `Import study`, which only
  `cdiscAppTB` (the app tree browser, run when the node expands) creates; `openApp` adds the study groups but
  not that item, so `items.find(...)` is undefined. README: "Click **Import study** tree item".
