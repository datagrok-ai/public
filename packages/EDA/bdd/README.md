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
| `models/train-on-cars.feature` | linear-regression, pls-regression, xgboost2 | the target, fifteen features picked in the column picker, three engines, a selection in the result chart, XGBoost's clickers and sliders |
| `models/train-on-iris.feature` | softmax, xgboost1 | the same for classification, Softmax's sliders, all five XGBoost controls |
| `models/share-model.feature` | share-model-permissions (owner's side) | a model trained and saved, its Sharing pane, the Share dialog, a share to the second account and its revoke |
| `pareto-front.feature` | pareto-front-viewer, steps 5–6 | the label picked on cars and on demog, and none on iris |
| `pareto-front-objectives.feature` | pareto-front-viewer, steps 1–4 and 7 | the property categories, the objectives offered, the min/max conflict warning, an axis and the labels chosen in the panel; `@known-failure` for an empty column |

Not translated: the recipient's side of `share-model-permissions` (seeing, applying and being refused
the shared model), which needs a second signed-in session; the Spaces features share the same way.

## Running

```bash
cd public/libraries/bdd && npm ci && npm run build && npx playwright install chromium
cd ../../packages/EDA && npm ci && npx grok-bdd link
DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev npx grok-bdd run --reporter=list
npx grok-bdd run generated/models/train-on-cars.test.ts   # one feature
```

Eleven features compile to thirteen Playwright tests. The stand needs EDA published and `cars.csv`,
`demog.csv` and `iris.csv` in `System:DemoFiles`. Sharing needs a dev key for setup of `bddsecond`.
The feature uses a unique `{run}` model name, disables notifications, and removes its model,
sharing wrapper and newly created training table at teardown. Cleanup failure fails the run.

The model tests require the core Train Model preview readiness attributes (`aria-busy` and
`aria-invalid`). Rebuild the Dart client after the companion core changes before running these tests. Model
artifact cleanup also needs the core project-deletion fix, which checks permission before removing
the wrapper relation. Missing-help HTML rejection and IPv4 API routing in the host-dev nginx
configuration fix two local-stand failures uncovered by these tests.
The filtered group-count defect (GROK-20795) and the empty Pareto objective offer remain narrowly
marked `@known-failure`: setup and positive result checks must pass normally.

## What the platform gave these features

- The Train Model view trains on every change. The preview reports when the latest training,
  predictions, charts and history finish; a card's parameter text alone can precede completion.
  Each scenario establishes its hyperparameters and exercises the resulting chart selection.
- The column picker, the property grid categories, the two-second current-object drop after a
  property edit and the Share dialog of a model are platform facts: `libraries/bdd/CLAUDE.md`,
  "Facts that cost a run each".
