# Proteomics import features — what is still missing

## Proteomics import

Feature: `import-files.feature` (30 Oct 2026, local stand).

### Questions (Proteomics import)

- **No help page.** The Proteomics package has no page under `public/help/`; the feature follows the package
  README ("Supported inputs", "Import"). Should the import (and the pipeline after it) get a help page, so a
  feature has an "expected" that is not the package's own README?
- **The table view's ribbon Proteomics menu.** `menu.ts` builds a Proteomics group on a Proteomics table
  view's ribbon (Import, Annotate, Analyze, Visualize, Share, items greyed out until their step applies). After
  the import the claims stop at the table; the ribbon menu and its grey-out are the next feature's subject
  (each step is a dialog: Annotate Experiment, Normalize, Impute, Differential Expression).

### Not written (Proteomics import)

- **Spectronaut Report / Candidates, FragPipe**: the same file chooser path as MaxQuant...; what each parser keeps
  and computes is the package's parser tests (`src/tests`). The Candidates import opens its own analysis view
  (`openCandidatesAnalysisView`), which would be the one worth a scenario.
