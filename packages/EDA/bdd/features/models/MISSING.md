# Models — what the features could not say yet

Translated on 2026-10-09 with the bindings as they are; the review of 2026-10-11 added a library step
(a click with a key held) and two core fixes (§2). Sources:
`files/TestTrack/Models/*.md` and `playwright-public/Models/*.test.ts`. Probed and run on a local stand
built from core master `1847d054e8`.

## Resolved

- **Ctrl+click two cards, then Compare** (models-testdemog-lifecycle-smoke block 4, GROK-19550): the library has
  `user clicks on {element} holding {key}` now, and `models-gallery` ends on the comparison. Under the test a
  card clicked within 2 s of expanding a context pane (Performance, before Run Evaluation) did not become
  current: the accordion's guard (`AppEvents.propertyEdited`) dropped it, as it drops a change a property edit
  causes. A click on a card is deliberate; `ui.setUpSelection` (d4 `ui.dart`) releases the guard, as a click in
  the grid already does.

- **Skip unique categories trained nothing.** Ticked, it offered the engines and enabled SAVE over an empty
  preview (`Invalid argument (namesOrColumns): Not supported type: null` from `PredictiveModelingEngine.apply`),
  and a model saved that way failed on apply with `Incorrect parameters count` (EDA `orderedFeatures`): the
  model listed the skipped column among its inputs while its coefficients were for the rest. Skipping is a
  choice of features made at training now — `skipUniqueCategories` records the columns it drops and
  `fillModelParameters` leaves them out of the model's inputs, so nothing is skipped again on apply.
  `train-data-checks` ticks it and trains; the model's inputs leaving the column out is the xamgle test
  `Predictive modeling | Train View | Skip unique categories` (saved and applied to the same table by hand:
  Inputs 1/1, it predicts).
- **Tip typos**: "contains contain too many unique categories" and "No models registred" are fixed
  (`predictive_modeling_validators.dart`, `predictive_modeling_view.dart`); `train-data-checks` claims the
  whole tip.
- **The training table of a model deleted through the gallery** (2026-10-12): `N predictive model(s) named … should be
  on the server` records each model's training table and wrapper project for the feature-end cleanup
  (`recordModelArtifacts`, library `platform/steps.ts`), so they go whatever deletes the model; the EDA model
  features leave no table, project or model behind. The server's model delete keeps its `trainedOn` table — by
  design, the lead's call (2026-10-11).
- **A signal for the Positive class cutoff re-render** (2026-10-12): the preview is `aria-busy` from the move until
  the debounced re-render ends (`predictive_modeling_previews.dart`); `train-options` claims Accuracy 0.800 at a
  cutoff of 0.8. Earlier previews no longer listen to the cutoff once replaced (`PredictiveModelPreview.dispose`).
- **A second model saved from the same Train Model view took the first one's name** (2026-10-12): the name was kept
  from the server's answer to the first save; it follows the friendly name again (`fillModelParameters`). xamgle
  test `Train View | Second model saved takes its own name`.
- **One-hot on a table that lacks a category** (2026-10-12): confirmed — a table with featureA = Yes only handed
  the model `featureA=Yes` alone. Training records each input's categories (`options['oneHotCategories']`) and
  apply rebuilds the training columns from them. xamgle test `Apply | One-hot rebuilds the training categories`.
- **§3, the API tests** (2026-10-12): the data checks and Ignore missing are xamgle tests ("Core: xamgle:
  Predictive modeling", 14/14); engine discovery, `/ml/zip`, blobs, images and save / suggested / delete are
  ApiTests `Dapi: models`. One-hot on apply is the xamgle test above (the one-hot logic lives in core, and the JS
  API cannot make a model with preprocessing options).

## Still open

Nothing planned.

## Won't do

- **Is applicable to...** (the old §1.1 and §2.1): the filter only writes `isApplicableTo = "<id>"` into the search,
  which neither narrows the gallery nor sorts it by applicability, and saves the chosen table again (a foreign-key
  error the second time). Filed on GitHub by Olesia; the feature waits for the answer to "filter or sort".
- **"Trained on" shows an id** (the old §2.2): an enhancement, filed by Olesia.
- **The old cases' expectations** (the old §2.5): models-validators-edge ("a categorical feature does not block
  training"), models-one-hot-suffix-collision (inputs listed as `featureA=Yes, …`) and models-lifecycle-csv-table
  (a save-success notification) describe what the platform does not do; the features claim what it does, and the
  md cases are their owner's to update.
- **Chemprop** (`chemprop.md`, `chemprop-ui.md`): trains and predicts in the chem-chemprop Docker container.
- **`predictive-models.md` scenario 3** (apply to Tools > Dev > Open test dataset, random walk): the random-walk
  columns are none of the model's inputs, so Apply offers no model.
- **`browse/modelhub.test.ts`**: the Compute Model Hub, not predictive models — the Browse and DiffStudio features
  cover it.

### Notes

- **Apply to from a card** (the old §2.4) — not a finding: walked by hand on 2026-10-09, the panel follows the
  dialog's model.
- **A strong negative correlation is never reported** by the Highly correlated check (`highlyCorrelated` flags
  only r > 0.9) — a candidate finding for the models' owner.
