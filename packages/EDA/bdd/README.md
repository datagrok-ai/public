# EDA — BDD features

The Test Track cases of EDA (`packages/UsageAnalysis/files/TestTrack/EDA/**/*.md`) and the package's
Playwright specs (`playwright/`), as features on `@datagrok-libraries/bdd`. Each feature says in its
description what it translates and what it left out. Every step is the library's: this project has
no bindings of its own.

| Feature | Case | What it walks |
|---|---|---|
| `analyze/pca.feature` | pca | the Features picker's All, three components, then Center and Scale |
| `analyze/pls.feature` | pls | the dialog, All refused while price is both target and predictor, the run without it |
| `analyze/multivariate-analysis.feature` | multivariate-analysis | the run, its columns and tables, the selection shared by the grid and the scatter plots, a coefficient bar selecting its loading |
| `analyze/anova.feature` | anova | the dialog, the conclusion on the box plot, the table of the test |
| `analyze/control-comparisons.feature` | package spec only | the dialog, the box plot, the per-group sizes of the result table |
| `analyze/filtered-group-comparison.feature` | GROK-20795 | `@known-failure`: a filtered run should count the filtered rows |
| `models/train-on-cars.feature` | linear-regression, pls-regression, xgboost2 | the target, fifteen features picked in the column picker, three engines, XGBoost's clickers and sliders |
| `models/train-on-iris.feature` | softmax, xgboost1 | the same for classification, Softmax's sliders, XGBoost's clicker and slider |
| `pareto-front.feature` | pareto-front-viewer, steps 5–7 | the label picked on cars and on demog, an axis set by hand |
| `pareto-front-objectives.feature` | pareto-front-viewer, steps 2–4 and 7 | the property categories, the objectives offered, the min/max conflict warning; `@known-failure` for an empty column |

Not translated: `playwright-public/Sharing/share-model-permissions.md` and the package's
`share-model-permissions.test.ts`, which need a second signed-in user the library does not have.

## Running

```bash
cd public/libraries/bdd && npm ci && npm run build && npx playwright install chromium
cd ../../packages/EDA && npm ci && npx grok-bdd link
DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev npx grok-bdd run --reporter=list
npx grok-bdd run generated/models/train-on-cars.test.ts   # one feature
```

Ten features, 26 scenarios: under a minute on four workers against dev (2026-09-11), and 28 of 28
with `--repeat-each=2` three times in a row. The stand needs EDA published and `cars.csv`,
`demog.csv` and `iris.csv` in `System:DemoFiles`. Nothing is saved to the server.

## What the platform gave these features

- The column picker ("Select columns...") is a grid viewer: `text of cell N of __name` names the
  column of a row, `cell N of x` is its checkbox. Its Search input filters the rows; a row keeps its
  number.
- The Train Model view trains on every change. The model card reports the parameters the model was
  trained with and its scores, which is what the claims read; the view remembers the last session's
  hyperparameters, so each scenario sets what it reads.
- A Dart property grid category tells its state only by its icon; the library's expand reads it.
- For two seconds after a property is edited in the context panel the platform ignores a change of
  the current object, which is why the objectives are one journey over one viewer.
