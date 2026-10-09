# Models — what the features could not say yet

Translated on 2026-10-09 with the bindings as they are (no library, binding or core change). Sources:
`files/TestTrack/Models/*.md` and `playwright-public/Models/*.test.ts`. Probed and run on a local stand
built from core master `1847d054e8`.

## 1. Steps wanted

### 1.1 A click with a key held, on any element — Ctrl+click two gallery cards, then Compare

Wanted: `When user clicks on {element} holding {key}` in `bindings/common/steps.ts` (a Playwright
`click({modifiers})` on the located element; the area version already exists for viewers). With it,
models-testdemog-lifecycle-smoke.md block 4 (GROK-19550, multi-select + Compare) is:

```gherkin
  Scenario: Two cards picked with Ctrl+click are compared side by side
    When user clicks on "BDD-Kestrel-{run}" label in gallery
    And user clicks on "BDD-Osprey-{run}" label in gallery holding Control
    Then the context panel should show "2 models"
    And "Actions" pane in context panel should contain text "Compare"
    When user clicks on "Actions" pane in context panel
    And user clicks on "Compare" label in "Actions" pane in context panel
    Then the "Compare models" view should be current
    And the table should have 2 rows
    And the table should have the columns "Name, Description, Method, Source"
    And no errors should have been logged
```

Checked by hand in a probe (Playwright `click({modifiers: ['Control']})`): the panel shows "2 models"
with Items 2 and Actions (Compare, Share..., Tag..., Edit..., Delete, Save as Zip, Create URL
alias...), and Compare opens "Compare models" with exactly Name, Description, Method, Source — not the
per-metric and image columns the old md expected.

### 1.2 Removing a table the UI saved to the server

The Is applicable to... quick filter of the gallery saves the chosen open table to the server when it is
not there yet (`predictive_models_view.dart` `_isApplicableFilter`, `dapi.tables.save`). No step deletes
such a table at feature end, so the filter is not translated. Wanted: `Given the tables named {string}
are deleted now and at feature end` (the library already lists tables as a cleanup source in
`serverEntities`). With it, and once 2.2 is answered:

```gherkin
  Scenario: Is applicable to... keeps the models that fit the chosen table
    Given the tables named "readings, levels" are deleted now and at feature end
    When user clicks on "Toggle filters" icon
    And user clicks on "Is applicable to..." tag
    Then the open menu should list "levels"
    When user picks "levels" from the open menu
    Then "BDD-Osprey-{run}" label in gallery should be visible
    And "BDD-Kestrel-{run}" label in gallery should be absent
```

### 1.3 The training table of a model deleted through the gallery

The model sweep (`namedCleanup`, `bindings/platform/steps.ts`) finds a model's uploaded training table
and its wrapper project through the live model. A feature that deletes its model through the gallery
(apply-and-delete does; models-gallery leaves Delete to it for that reason) leaves both on the server: the server's model delete
(`ml_service.dart` `deleteModel`) keeps the `trainedOn` table. Wanted: `N predictive model(s) named …
should be on the server` (or the save itself) records the training table and the wrapper project for
the feature-end cleanup, so they go whatever deletes the model.

### 1.4 A signal for the Positive class cutoff re-render

Moving Predict probability's cutoff re-renders the confusion matrix and the metrics after a 200 ms
debounce (`predictive_modeling_previews.dart`) without setting `aria-busy` on the preview, so no claim
can wait for it. With the busy state set around it:

```gherkin
    When user enters "0.8" into "Positive class cutoff" input
    Then model preview should be ready
    And "Accuracy" table row should not contain text "0.900"
```

## 2. To check by hand (candidate findings — nothing filed)

### 2.1 Skip unique categories: no model trains after it is ticked

Steps:
1. Open a table with an identifier-like text column and a numeric one, for example paste
   `id_like,feature1,target` / `row_1,1,3.1` / `row_2,2,1.2` / `row_3,3,4.4` / `row_4,4,2.5` /
   `row_5,5,5.3` / `row_6,6,1.7` as a new table.
2. ML > Models > Train Model..., Predict = target, Features = id_like and feature1.
3. Insights & Tips says id_like has too many unique categories; tick **Skip unique categories**.

Expected: id_like is dropped and the model trains on feature1 (preview with metrics, SAVE enabled).
Actual: Model Engine and its parameters appear and SAVE turns active, but the preview stays empty
and the console logs `Invalid argument (namesOrColumns): Not supported type: null` from
`ColumnList.toColumnList` in `PredictiveModelingEngine.apply` (`predictive_modeling_engines.dart:58`).
Picking another engine changes nothing. Walked by hand 2026-10-09: the same error and stack
(Screenshot_40); SAVE is enabled over a model that never trained.

Why (from the code): the preview calls `apply` with a column map built from the reduced features
(feature1 only), whose length differs from `model.input` (id_like, feature1); `apply` then rebuilds the
map from `model.input`, finds no id_like in the reduced frame and maps it to null. One-hot encoding goes
through the same path but is caught by the `<name>=` lookup.

Ready for the feature once a ticket exists (train-data-checks, after the tip):

```gherkin
    # GROK-NNNNN
    @known-failure
    Scenario: Skipping the identifier-like column lets the model train
      When user checks "Skip unique categories" input
      Then model preview should be ready
      And "R squared" table row should be visible
      And no errors should have been logged
```

### 2.2 Is applicable to... does not narrow the gallery, and saves the open table

Steps: train and save one model on table A (columns f1, f2, target) and one on table B (g1, g2,
score), keep both tables open, Browse > Platform > Predictive models, the filter icon, **Is applicable
to...**, pick B.

Expected: only the model trained on B. Actual (local stand): the search reads `isApplicableTo = "<id>"`
and the counter stays 2. The filter only writes that text into the search, which the server reads as
a smart filter; nothing turns it into the `tableId` parameter of `getModels` (`ml_service.dart`), and
even that branch only sorts the list, with a one-argument comparator
(`models.sort((pmi, _) => pmi.isApplicable(tableInfo) ? 0 : 1)`). Picking a table that was never saved also saves it to the server (1.2). Whether "Is applicable
to" was meant to filter or to sort is the question.

### 2.3 Text of the tips

- "Column 'id_like' contains **contain** too many unique categories." (`tooManyUniqueCategories`,
  `predictive_modeling_validators.dart`).
- "No models **registred** that match this type of data" under the Features when no engine fits.

### 2.4 Details pane: "Trained on" shows an id

A model trained in the view and saved shows `Trained on 12a571e0-bde0-…` in its Details pane, not a
table name (the training table is built without a name in `fillModelParameters`). Intended?

### 2.5 A second model saved from the same Train Model view takes the first one's name

Train, save as `A-LR`; change Model Engine, save as `A-PLS`. On the server the second model is
`name = ALR_1`, `friendlyName = A-PLS` (its card is `div-ALR-1`). The friendly name shows everywhere a
person looks, so this matters only where `name` is used (URL, scripting by name). Intended?

### 2.6 Apply to from a card leaves the previous model current (low confidence)

Steps: in Predictive models click model A (the context panel shows A), right-click model B's card,
Apply to > a table B fits. Expected (`runPredictiveModellingApply` sets `AppEvents.currentObject` to the
dialog's model): the context panel shows B while the dialog is open. Actual in the feature run: the
current object stayed `Model "A"`. Changing the Model choice inside the dialog does make the new model
current (apply-and-delete), so only the dialog's first model is affected. Not instrumented yet — the
context panel's own guards may drop the change; walk it by hand before anything else.

### 2.7 Questions on the old cases' expectations

- models-validators-edge.md says a categorical feature "does not block training". With no engine
  taking a text feature, the view shows no Model Engine and SAVE stays disabled until One-hot
  encoding is ticked. The feature claims that; the md needs updating if it is intended.
- models-one-hot-suffix-collision.md expects the saved model's inputs as `featureA=Yes, featureA=No,
  …`. The model keeps the table's columns (featureA, featureB) as its inputs and expands them inside;
  Apply maps (2/2) and predicts. The feature claims the latter.
- models-lifecycle-csv-table.md expects a "save-success notification" after Save; none is shown.
- One-hot on a fresh table that lacks a category: `oneHotEncoded` builds the 0/1 columns from the
  apply table's own categories, so a table whose featureA holds only Yes may hand the model one input
  short. Not probed yet — apply the One-hot model of train-options to a table with featureA = Yes only.

## 3. API tests wanted (none exists — rule "Anything an API test can check")

| Behaviour | Where | Old source |
|---|---|---|
| Engine discovery: functions with `meta.mlrole` train/apply pair up by `mlname` | ApiTests `src/dapi/models.ts` (new) | models-engines-discovery-api |
| `/ml/zip/{id}` is a zip; blobs, images and the image list round-trip | ApiTests `src/dapi/models.ts` | models-api-misc-api |
| `dapi.models` save / suggested / delete; a deleted model's blobs go | ApiTests `src/dapi/models.ts` | the dapi count checks of every Models spec |
| The data checks (class imbalance ±20 %, categorical features, >0.8 unique, Pearson > 0.9, missing values) | xamgle Dart test beside `predictive_modeling_validators.dart` | models-validators-edge, GROK-3525 |
| One-hot on apply rebuilds `<name>=<category>` columns from the table's columns | EDA `src/tests/model-serialization-tests.ts` | models-one-hot-suffix-collision |
| Ignore missing trains on the rows left (10 rows with 2 missing targets → 8), and the per-column counts the missing-values tip prints for two or more columns | xamgle Dart test beside `predictive_modeling_validators.dart` | models-bug-grok-3525 |

## 4. Not translated, with the reason

- Chemprop (`chemprop.md`, `chemprop-ui.md`): trains and predicts in the chem-chemprop Docker container.
- `predictive-models.md` scenario 3 (apply to Tools > Dev > Open test dataset, random walk): the
  random-walk columns are none of the model's inputs, so Apply offers no model; the old spec fell back
  to the first option and asserted only a column count.
- `browse/modelhub.test.ts`: the Compute Model Hub, not predictive models — the Browse and DiffStudio
  features cover it.
