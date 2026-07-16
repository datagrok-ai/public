---
feature: models
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [train-classification-end-to-end, train_and_save_model, apply_model_on_table, compare-multiple-models, delete_model]
realizes: [views.models]
realized_as:
  - models-testdemog-lifecycle-smoke-spec.ts
related_bugs: [GROK-2381, GROK-19177, GROK-19550, GROK-846]
---

# TestDemog predictive model lifecycle (smoke)

A smoke run of the full predictive-model lifecycle — train, apply,
browse/search/compare, and delete — on a single model (**TestDemog**)
trained on the public `demog.csv` demo dataset.

## Setup

A clean Datagrok session. All blocks operate on `demog.csv` from
`System:DemoFiles`. Block 1 trains and saves **TestDemog**; Blocks 3–5 consume
it. Downstream specs are expected to provision the TestDemog fixture in
`beforeAll` via `js-api-replay` rather than serializing UI dependence on Block 1.

For Block 4, ≥1 additional model must be present in the catalog alongside
TestDemog. The spec MUST provision a second model in `beforeAll`.

## Scenarios

### Block 1: Train classification model (TestDemog)

1. Open `demog.csv` from **System:DemoFiles** (or by clicking the star icon if
   pinned).
2. Go to **ML > Models > Train Model...** — the `PredictiveModelingView` opens
   (left dock = parameters form; right dock = Interactive modeling preview).
3. Configure: **Predict** `SEX` (categorical / 2-class), **Features** `WEIGHT`
   and `HEIGHT`. **Immediately tick Ignore missing.** **Verify:** the Model
   Engine field appears (e.g. `Eda: XGBoost`) and the **SAVE** button becomes
   enabled — the view auto-trains once a missing-value strategy is chosen.

   > **UI behaviour note:** the Model Engine field and SAVE button remain
   > disabled until either **Ignore missing** or **Impute missing** is ticked.
   > Ticking **Ignore missing** hides **Impute missing** and triggers automatic
   > training; ticking **Impute missing** first opens the KNN config dialog
   > instead.

4. **Verify:** `Impute missing` is no longer visible (hidden once `Ignore missing`
   is checked). `Predict probability` is visible for the 2-class categorical
   target `SEX`.
5. Tick **Predict probability**. **Verify:** re-train completes; preview shows
   probability / threshold output for the SEX target.
6. Save the model as **TestDemog**. **Verify:** the model is discoverable in
   **Browse > Platform > Predictive models**.

> **Bug anchor — GROK-2381:** Step 6 is the post-train notification surface
> (success path). No bug-focused spec asserts the failure notification path
> yet.

### Block 2: Train regression model (numerical target)

1. Reuse the same `demog.csv` table.
2. Go to **ML > Models > Train Model...** again.
3. Configure: **Predict** `WEIGHT` (numerical), **Features** `HEIGHT`. **Tick
   Ignore missing.** **Verify:** engine loads and SAVE button becomes enabled.
4. **Verify:** `Predict probability` is NOT visible (gated on categorical
   targets only).
5. Save under a distinct name (suggested: `TestDemog_Regression`).

> **Predict probability** is gated on 2-category categorical targets;
> do not toggle it on this block.

### Block 3: Apply TestDemog model

*Prerequisite: TestDemog exists (Block 1 or `beforeAll` `js-api-replay` fixture).*

1. Open `demog.csv` from **System:DemoFiles**. Capture the initial column count.
2. Go to **ML > Models > Apply Model...** — the **Apply predictive model** dialog
   opens. **Verify:** dialog title reads "Apply predictive model"; Model select is
   populated via `dapi.ml.suggested(tableInfo)`; the **OK** button is present.
3. Select **TestDemog** from the Model select. **Verify:** the `ColumnsMapInput`
   reflects inputs `WEIGHT`, `HEIGHT` mapped against `demog.csv` columns.
4. Click **OK**. **Verify:** a new prediction column appears; column count is
   strictly greater than the initial count; the appended column carries the
   `Tags.PredictiveModel` markup.

> **Bug anchor — GROK-19177:** Step 2 is the dialog-open anchor (happy path).
> Bug-focused spec asserts the empty-models-list guard (0 suggested entries →
> OK is blocked, not an unhandled exception).

### Block 4: Browse, search, context panel, filter templates, multi-select Compare

*Prerequisite: TestDemog exists + ≥1 additional model in the catalog.*

1. Navigate to **Browse > Platform > Predictive models**.
2. Type `TestDemog` in the search field. **Verify:** the TestDemog card surfaces.
3. On the **Context Panel**, walk every tab and verify it renders model metadata.
4. Click the **Filter templates** icon. **Verify:** the filter-templates panel
   opens with content.
5. Clear the search field so all catalog models are visible.
6. Verify ≥2 models are present (precondition for multi-select).
7. **CTRL+click** at least two model cards. **Verify:** ≥2 cards are selected.
   Assert via the JS-side selection model — DOM-class assertion is unreliable
   post-fix per GROK-19550 refdoc gotcha (no `.selected` class applied).
8. On the **Context Pane**, open the **Commands** tab and click **Compare**.
   **Verify:** a new view opens with a "Compare models" DataFrame (Name /
   Description / Method / Source + per-metric columns + image cells at
   width 200).

> **Bug anchor — GROK-19550:** Step 7 is the multi-select regression target.
> Assert via JS selection model or the downstream Compare result table.

### Block 5: Delete TestDemog

*Prerequisite: TestDemog exists. Section-terminal cleanup — `delete-spec.ts`
MUST run last among the four realized specs.*

1. Navigate to **Browse > Platform > Predictive models** (`type: models`,
   path `/models`).
2. Locate **TestDemog**. **Verify:** at least one card with name "TestDemog" is
   present.
3. Right-click the **TestDemog** card and select **Delete**. **Verify:** a
   confirm-delete modal opens.
4. Click **Delete** in the confirmation dialog. **Verify:** dialog closes; no
   error balloon or console exception.
5. **Verify:** the TestDemog card no longer appears in the Predictive models
   browser.

> **Bug anchor — GROK-846:** Step 4 is the server-side delete trigger
> (UI-visible arc only). The FK-constraint state-cleanup invariant has no
> dedicated spec yet.