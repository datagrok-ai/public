# Plates app — what the features could not say yet

## Plates

### Needs a step (Plates)

- **CREATE with an imported plate, SAVE of a new template, Clone / Delete of a template.** Each writes
  plates or templates into the package database (`plts` schema), and no library step removes them. A
  step pair is needed, e.g. `the plates with barcode starting with {string} are removed at feature end`
  and `the plate template {string} is removed now and at feature end` (deleting through
  `Plates:` CRUD functions, read back gone), plus a CSV fixture for the Create view's import. Then:

  ```gherkin
  Scenario: CREATE saves an imported plate, and Apps > Plates lists it
    Given the plate template "bdd-template-{time}" is removed now and at feature end
    When user clicks on Apps---Plates---Create tree node inside browse tree
    And user uploads "fixtures/bdd-plate.csv" through "Import" input
    And user clicks on "CREATE" button
    Then an info balloon containing "Plate Created Successfully" should have been shown
    And no errors should have been logged
  ```

- **`user opens the {string} app` does not find Plates.** The step looks the app up by
  `DG.Func.find({tags: ['app']})`; Plates' app is declared by the `@grok.decorators.app` role
  (`options.role = "app"`, no `#app` tag), so the step answers "no Plates app is registered". The
  features open the app from the Browse tree instead; the step should also take functions whose role is
  `app`.

### Questions (Plates)

- **Search Plates' Template picker does not filter the plates found.** With "Cell counting" or
  "Dose-response" picked the result grid keeps all 40 plates; only the property inputs change
  (`plates-search-view.ts` `refreshResults` builds matchers from the inputs only). Is the template meant
  to restrict the search to its plates? If it is, the scenario claims 20 after the pick.
- **The Assay Plates demo needs `System:DemoFiles/hts/xlsx_plates`**, which the local stand lacks (the
  `data` submodule holds no `hts` folder). The demo (`assayPlatesDemo`) and the Excel plate viewer are
  untested until the stand has those files.
