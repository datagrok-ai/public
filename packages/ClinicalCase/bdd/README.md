# ClinicalCase — BDD features

The start of the Clinical Case application and the package's SDTM folder viewer, as features on
`@datagrok-libraries/bdd`. Every step is the library's: this project has no bindings of its own.

| Feature | What it walks |
|---|---|
| `app/clinical-case-app.feature` | Browse > Apps > Clinical Case: the start page with the demo studies CDISCPILOT01 (306 subjects) and nida-ctn-0001 (138), the studies and Import study under the tree node, the Import study view with Import disabled and its error for a file without a study definition |
| `files/sdtm-folder.feature` | Browse > Files > App Data > ClinicalCase > SDTM > CDISCPILOT01: the folder viewer's "Folder contains SDTM data", the study id from dm.csv, the Run ClinicalCase button |

Listing the studies writes `studyConfig.json` into each study folder of `System:AppData/ClinicalCase/SDTM/`;
the app feature puts the same text there right after the listing through
`a file {string} with text {string} is on the server`, which deletes both files at feature end.

Still out, with the reason, in `features/app/MISSING.md`: opening a study and every view under it (the
CDISC rules engine runs in the package's container, and `<study>.d42` and `validation_results.json` are
written into the study folder with no library step to take them away), Import, Run ClinicalCase, and a
finding (Import study... on the start page closes the app when the tree node was never expanded).

## Running

```bash
cd public/packages/ClinicalCase/bdd
node ../../../libraries/bdd/bin/grok-bdd.js compile
DATAGROK_URL=http://localhost:8888 node ../../../libraries/bdd/bin/grok-bdd.js run
```

The package must be published on the stand (`grok publish localhost --release`).
