# Apps — what the features leave out

## Won't do

- Edge inputs in the demos and guides: they show a happy path; the viewers' edge cases are their own features'.
- What happens inside an app (Hit Design's rows, an app's own dialogs): the app's own features.
- **Plates**: the plate list and search results (the seed differs per stand: it adds the Excel plates where DemoFiles
  holds `hts/xlsx_plates`); CREATE with a plate and template SAVE / Clone / Delete (they write into the package
  database); a well hovered or clicked (no well areas); Excel plates and the Assay Plates demo. Readers and the search
  syntax are package tests (`plate-tests.ts`).
- **Clinical Case**: the app (start page, studies tree, Import study, a study's views). Listing writes
  `studyConfig.json` into `System:AppData/ClinicalCase`, and a study runs the CDISC rules engine in the package's
  container. The app's opening is the package test `clinicalCase.appOpensWithountExceptions`; its SDTM folder launcher:
  `browse/MISSING.md`.
- **Model Hub (Compute2)**: the gallery search (the platform gallery's, EDA `models/models-gallery.feature`),
  Tools > Browse Runs... and comparing runs (Compute2 tests `Apps: Browse runs`, `RunComparison: *`). Compare Runs...
  opening is in `browse/browse-apps-and-dashboards.feature`.
- **Proteomics Import** (MaxQuant, Generic Matrix, Spectronaut, FragPipe): parsing, filtering, column suggestions and
  semantic types are the package's tests (`src/tests/parsers.ts`, `generic-parser.ts`). The menu and file chooser are
  claimed by Jmp `open.feature` and PowerPack `io/xlsx-open.feature`.
