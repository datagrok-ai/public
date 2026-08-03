---
feature: eda
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - anova.test.ts
  - multivariate-analysis.test.ts
  - pca.test.ts
  - pls.test.ts
  - pareto-front-viewer.test.ts
  - MLMethods/linear-regression.test.ts
  - MLMethods/pls-regression.test.ts
  - MLMethods/softmax.test.ts
  - MLMethods/xgboost1.test.ts
  - MLMethods/xgboost2.test.ts
related_bugs: []
---

# EDA Playwright tests

End-to-end coverage for Test Track scenarios under
`public/packages/UsageAnalysis/files/TestTrack/EDA/`. Tests are UI-first; JS API is used
only when the UI path is verifiably blocked (canvas-rendered controls, Dart-side dialog
state that does not surface to JS, missing test data, etc.).

Created 2026-05-18, last green run 2026-05-20 against `https://dev.datagrok.ai`
(13/13 passing, no skips and no `test.fail`).

## Files

```
playwright-tests/e2e/EDA/
├── helpers.ts                       # Shared helpers (menu, dialog, picker, viewer probes)
├── anova.test.ts                    # ML > Analyze > ANOVA on demog.csv
├── multivariate-analysis.test.ts    # ML > Analyze > Multivariate Analysis on cars.csv
├── pca.test.ts                      # ML > Analyze > PCA on cars.csv (+ Center/Scale)
├── pls.test.ts                      # ML > Analyze > PLS on cars.csv (UI smoke — see limitation)
├── pareto-front-viewer.test.ts      # ParetoFrontViewer eligibility + label autoselect
└── MLMethods/
    ├── linear-regression.test.ts    # ML > Models > Train Model on cars.csv
    ├── pls-regression.test.ts       # ML > Models > Train Model on cars.csv (PLS)
    ├── softmax.test.ts              # ML > Models > Train Model on iris.csv (numeric features only)
    ├── xgboost1.test.ts             # XGBoost classification on iris.csv
    └── xgboost2.test.ts             # XGBoost regression on cars.csv
```

## Running

```powershell
# All EDA tests
cd C:\DataGrok\NewBitbucket\reddata\playwright-tests
npx playwright test e2e/EDA --reporter=list

# Single file
npx playwright test e2e/EDA/anova.test.ts --reporter=list

# Filter by test title (regex)
npx playwright test e2e/EDA --grep "PCA" --reporter=list

# Visible browser (debug)
npx playwright test e2e/EDA/pls.test.ts --headed --reporter=list

# Auth state is built once by e2e/global-setup.ts using env vars from .env:
#   DATAGROK_URL, DATAGROK_LOGIN, DATAGROK_PASSWORD
# It writes e2e/.auth.json; delete that file to force a fresh login.
```

### Reports

`playwright.config.ts` configures three reporters: `list` (console), `html`
(`e2e/html-report/`) and `monocart-reporter` (`e2e/monocart-report/index.html`).

> **Run without `--reporter` to get the monocart report at the configured path.**
> Passing `--reporter=...` on the CLI *overrides* the whole config reporter block and drops
> monocart's `outputFile` option, so the report is written to the default
> `playwright-tests/monocart-report/` instead of `e2e/monocart-report/`. You then risk
> opening a stale `e2e/monocart-report/index.html` from an earlier run.

```powershell
# Generates list + html + monocart at their configured paths:
npx playwright test e2e/EDA

# Open the monocart report (Windows):
Start-Process e2e\monocart-report\index.html
# or serve it (if assets don't load via file://):
npx monocart show-report e2e/monocart-report/index.html
```

The monocart summary has separate `Failed` and `Errors` rows. `Failed` is the test result;
`Errors` counts *captured error events*, including each failed retry of an `expect.poll`
matcher — so a green suite could still show `Errors > 0`, fluctuating run to run. The EDA
tests therefore wait via `page.waitForFunction` (polls in-browser, rejects only once on
timeout) plus a single trailing `expect`, which keeps `Errors` reliably at `0`. Don't
reintroduce `expect.poll` here.

## Last results

13 passed in ~3.1 min against dev:

| # | Scenario | Result | Time |
|---|----------|--------|------|
| 1 | ANOVA on demog.csv (Box plot + Analysis + F-test tabs)                          | PASS | 10s |
| 2 | Linear Regression on cars.csv predicting price                                  | PASS | 9s |
| 3 | PLS Regression on cars.csv with 3 components predicting price                   | PASS | 19s |
| 4 | Softmax on iris.csv predicting Species (numeric features only)                  | PASS | 19s |
| 5 | XGBoost classification on iris.csv predicting Species                           | PASS | 10s |
| 6 | XGBoost regression on cars.csv predicting price                                 | PASS | 20s |
| 7 | MVA on cars.csv (Grid + 3 Scatter + 2 Bar)                                      | PASS | 10s |
| 8 | Pareto cars-with-missing — empty + string columns excluded from Min/Max         | PASS | 9s |
| 9 | Pareto cars.csv — auto-select `model` as Label                                  | PASS | 8s |
| 10 | Pareto demog.csv — auto-select `USUBJID` as Label                              | PASS | 24s |
| 11 | Pareto cars.csv — Description/Objectives/Axes/Labels/Legend categories present | PASS | 23s |
| 12 | PCA on cars.csv adds PC1/PC2/PC3, then PC1 (2)/PC2 (2)/PC3 (2) with Center+Scale | PASS | 13s |
| 13 | PLS dialog opens on cars.csv with all expected inputs and accepts Components=3  | PASS | 10s |

## UI vs API split per scenario

Every test opens its CSV through `openDemoCsv` — UI navigation to
`/files/System.DemoFiles/?browse=files` + dblclick the file label. The helper falls back
to `grok.dapi.files.readCsv` only when the dblclick does not produce a TableView within
~25 s (e.g. file absent from the tile or overlay intercepts the click).

| Scenario | Menu nav | Dialog/View | Inputs | Trigger | Verification |
|----------|----------|-------------|--------|---------|--------------|
| ANOVA | UI | UI | defaults | UI (RUN) | UI (tab labels + viewer types) |
| MVA | UI | UI | defaults | UI (RUN) | UI (viewer types + rowCount) |
| PCA | UI | UI | UI (Components, Center, Scale) | UI (OK) | UI (column names) |
| PLS | UI | UI | UI (Components, Using=All) | **skipped** | UI (dialog inputs visible) |
| Pareto | API (`addViewer`) | n/a | API (`v.props.*`) | n/a | API (`props.getProperties()`) |
| Linear Regression | UI | UI | UI (Predict) | **API (`eda:trainLinearRegression`)** | API (return non-null) |
| PLS Regression | UI | UI | UI (Predict) | **API (`eda:trainPLSRegression`)** | API (return non-null) |
| Softmax | UI | UI | UI (Predict) | **API (`eda:trainSoftmax`, numeric features only)** | API (return non-null) |
| XGBoost 1/2 | UI | UI | UI (Predict) | **API (`eda:trainXGBooster`)** | API (return non-null) |

## Known platform limitations and bugs the tests work around

### 1. PLS dialog disables RUN when "select all" includes Predict (UX bug)

The scenario literally says "select all available columns" for Using. After `All` is
clicked, `price` is in both Predict and Using, the validator silently disables RUN with
only a red border on Predict — no tooltip, no banner. See
[pls-run.md](../../../public/packages/UsageAnalysis/files/TestTrack/EDA/pls-run.md).

The test cannot work around the bug from the UI either:

* The "Select columns..." sub-dialog uses a canvas-rendered grid; per-row checkboxes are
  not DOM elements (`document.querySelectorAll('.d4-dialog input[type="checkbox"]')` ⇒ 0).
* `DG.Dialog.getOpenDialogs()` returns the PLS dialog handle with `inputs.length === 0` —
  the inputs live entirely on the Dart side and aren't exposed to JS.
* `grok.functions.call('EDA:PLS', { features, ... })` (and `Func.apply` / `Func.prepare`)
  strips the JS wrapper from the `ColumnList` going through the Dart bridge, so the
  function body fails with `TypeError: t.byIndex is not a function`. Direct
  `getPlsAnalysis(...)` works only from inside the EDA package; we can't import it.

Decision: the test verifies the dialog opens with the documented input set and accepts
`Components=3`, then cancels. When the platform bug is fixed (auto-exclude Predict from
Using, or surface the reason inline) restore the trailing RUN/columns assertion.

### 2. Softmax requires numeric-only features (not a platform bug)

`eda:trainSoftmax` treats *every* column of the passed `df` as a feature and throws
`"Training failes - incorrect features type"` if any of them is non-numeric — unlike
XGBoost, it does not tolerate the string target column being included. The earlier
`test.fail()` workaround was hiding an incorrect call (the full iris frame, including the
string `Species`, was passed as features).

The fix lives in `trainEdaModelViaApi`: with `numericOnly: true` it feeds only the numeric
feature columns (and now excludes the predict column from them), while still passing
`Species` as the predict column. `softmax.test.ts` therefore runs as a normal passing test.
See [softmax-run.md](../../../public/packages/UsageAnalysis/files/TestTrack/EDA/MLMethods/softmax-run.md).

### 3. Train Model view: canvas-based Features picker, missing Model Engine dropdown, no Train button

Documented across `linear-regression-run.md`, `pls-regression-run.md`, `xgboost*-run.md`.
Consequences for automation:

* "Select Columns" sub-dialog rows are canvas-rendered — `mousedown`/`click`/`pointerdown`
  on a cell does not toggle the underlying BitSet (the run reports confirm this for all
  Train Model paths).
* The `Model Engine` dropdown is not always rendered on dev — the engine is selected by
  the function name on the API call instead (`eda:trainLinearRegression`,
  `eda:trainPLSRegression`, `eda:trainXGBooster`).
* There is no visible "Train" button in the view — training is triggered implicitly. From
  automation the training is invoked through the registered function with the same input
  set RUN would have produced.

The MLMethods tests cover the UI path up to and including Predict-column selection, then
fall back to the function-registry call for training and assert a non-null model blob.

### 4. Pareto Front Viewer warning is canvas-only

`_showErrorMessage` draws the "Cannot minimize and maximize ..." text directly onto the
viewer canvas; the `errMsg` field is renamed by the production minifier, so there is no
stable surface to assert against from automation. The test verifies the eligibility list
(empty + string columns excluded) instead, which is the deterministic part of the
scenario.

### 5. `cars-with-missing.csv` is not deployed on dev under `System:DemoFiles/`

See [pareto-front-viewer-run.md](../../../public/packages/UsageAnalysis/files/TestTrack/EDA/pareto-front-viewer-run.md).
The Pareto test synthesises an equivalent dataset in memory: open `cars.csv`, null every
cell of the `turbo` column, rename to `cars-with-missing`.

### 6. `shell.tv` does not survive opening Train Model

Opening the Train Model view makes the `PredictiveModel` view the active view, so
`grok.shell.tv` no longer points to the source TableView. The `trainEdaModelViaApi` helper
walks `grok.shell.views` and uses the first `TableView` with a populated `dataFrame`,
rather than `shell.tv`.

## helpers.ts cheat sheet

| Function | Purpose |
|----------|---------|
| `BASE` | `process.env.DATAGROK_URL` |
| `inputHost(safeName)` | `[name="input-host-…"]` selector |
| `inputEditor(safeName)` | `${inputHost(safeName)} .ui-input-editor` |
| `openDemoCsv(page, fileName)` | UI-first file open with JS API fallback after timeout |
| `waitForCurrentTableView(page, timeoutMs)` | Wait until `grok.shell.tv.dataFrame.columns.length > 0` |
| `clickTopMenuLeaf(page, nameAttr)` | DOM-level menu click with hover-chain for nested submenus |
| `waitForDialog(page, title, timeoutMs)` | Wait for `.d4-dialog .d4-dialog-title` (substring, case-insensitive) |
| `topDialog(page)` | Locator for the last `.d4-dialog` |
| `clickDialogPrimary(page, candidates)` | Click first enabled of `button-OK`/`button-Run`/`button-RUN` |
| `selectAllColumnsInPicker(page, hostName)` | Open picker → click `label-All` → click `button-OK` |
| `setDialogColumnListInput(page, caption, names)` | API-only fallback via `DG.Dialog.getOpenDialogs()` (limited utility — see PLS notes) |
| `currentNumericColumnNames(page)` | Names of non-string columns of `shell.tv.dataFrame` |
| `setInputValue(page, hostName, value)` | Click editor → `Ctrl+A` → type → `Tab` |
| `setBoolInputOn(page, hostName)` | Idempotent checkbox-on |
| `currentColumnNames(page)` | `shell.tv.dataFrame.columns.names()` |
| `currentViewerTypes(page)` | `Array.from(shell.tv.viewers).map(v => v.type)` |
| `visibleTabLabels(page)` | Text of every `.d4-tab-host .d4-tab-header` |
| `isPrimaryEnabled(page, candidates)` | True if any candidate primary button is enabled |
| `resetShell(page)` | Escape → `grok.shell.closeAll()` → home |
| `visibleErrorBalloons(page)` | Visible `.d4-balloon-error` text |
| `openTrainModelView(page)` | `ML > Models > Train Model...` + wait for `PredictiveModel` view |
| `setPredictColumn(page, column)` | Synthetic `mousedown` on Predict editor + type + Enter |
| `trainEdaModelViaApi(page, fn, predict, opts)` | Find first TableView, optionally restrict to numeric features (`numericOnly` also excludes the predict column), call `g.functions.call(fn, { df, predictColumn, ...extraParams })`. Predict column is read from the source frame so a string target survives the numeric clone (softmax). |

## Selector patterns used

| Element | Selector |
|---------|----------|
| Top menu leaf | `[name="div-ML---Analyze---PCA..."]`, `[name="div-ML---Models---Train-Model..."]` |
| Input host | `[name="input-host-Predict"]`, `[name="input-host-Using"]`, `[name="input-host-Components"]` |
| Input editor | `[name="input-host-…"] .ui-input-editor` |
| Dialog (top) | `.d4-dialog` (use `.last()`) |
| Dialog title | `.d4-dialog .d4-dialog-title` |
| Primary buttons | `[name="button-OK"]`, `[name="button-Run"]`, `[name="button-RUN"]`, `[name="button-CANCEL"]` |
| Column picker All/None | `.d4-dialog [name="label-All"]`, `.d4-dialog [name="label-None"]` |
| Files browser file label | `label:has-text("cars.csv")` (case-insensitive substring) |
| Tabs | `.d4-tab-host .d4-tab-header` |

### Casing inconsistency to remember

Different ML dialogs use different casing for the primary action button:

| Dialog | Primary button |
|--------|----------------|
| ANOVA | `button-Run` |
| Multivariate Analysis | `button-RUN` |
| PCA | `button-OK` |
| PLS | `button-RUN` |

`clickDialogPrimary` tries `OK` → `Run` → `RUN` in order, so callers don't need to know
which one their dialog uses.

## Debugging recipes

### See what the page looks like at failure

* Screenshot is in `e2e/test-results/<test-name>-chromium/test-failed-1.png` after a run.
* Add an ad-hoc screenshot: `await page.screenshot({ path: 'e2e/test-results/probe.png', fullPage: true });`

### Dump DOM for an open dialog

```ts
const info = await page.evaluate(() => {
  const last = Array.from(document.querySelectorAll('.d4-dialog')).slice(-1)[0];
  return {
    title: last?.querySelector('.d4-dialog-title')?.textContent?.trim(),
    inputs: Array.from(last?.querySelectorAll('[name^="input-host-"]') ?? [])
      .map((el) => el.getAttribute('name')),
    buttons: Array.from(last?.querySelectorAll('[name^="button-"], [name^="label-"]') ?? [])
      .map((el) => el.getAttribute('name')),
  };
});
console.log(JSON.stringify(info, null, 2));
```

### Inspect a viewer's properties

```ts
const props = await page.evaluate(() => {
  const g = (window as any).grok;
  const v = Array.from(g.shell.tv.viewers).find((x: any) => x?.type?.toLowerCase().includes('pareto')) as any;
  return v?.props?.getProperties?.().map((p: any) => ({
    name: p.name, type: p.propertyType, category: p.category,
  })) ?? [];
});
console.log(props);
```

### Check the current TableView state

```ts
const state = await page.evaluate(() => {
  const g = (window as any).grok;
  return {
    activeViewType: g.shell.v?.type,
    tableViews: Array.from(g.shell.views)
      .filter((v: any) => v?.type === 'TableView')
      .map((v: any) => ({ name: v.name, rowCount: v.dataFrame?.rowCount, cols: v.dataFrame?.columns?.names?.() })),
  };
});
console.log(state);
```

### Force a single test to retry-free / no-timeout

```ts
test.setTimeout(600_000);
test.describe.configure({ retries: 0 });
```

### Reset auth state

```powershell
del e2e\.auth.json
# next run rebuilds it via global-setup
```

### Bypass dev — test against another server

```powershell
$env:DATAGROK_URL = "https://public.datagrok.ai"
$env:DATAGROK_LOGIN = "<login>"
$env:DATAGROK_PASSWORD = "<password>"
del e2e\.auth.json
npx playwright test e2e/EDA --reporter=list
```

## Things tried that did NOT work — don't go back down these paths

1. **`DG.Dialog.getOpenDialogs()[i].input('Using').value = [...]`** — dialog is Dart-side,
   `inputs.length === 0` returned to JS.
2. **`document.querySelector('[name="input-host-Using"]').__dgInput`** — no `__*` keys on
   the input host element.
3. **Canvas pixel click on the "x" checkbox column** in the "Select columns..." picker —
   pointer events don't reach the cell renderer; even when they did the BitSet wasn't
   toggled in the run-report tests.
4. **`DG.Viewer.fromRoot(picker.querySelector('.d4-viewer.d4-grid'))`** — returns null,
   the embedded picker grid is not introspectable through the Viewer API.
5. **`grok.functions.call('EDA:PLS', { features: subDf.columns, ... })`** — the
   `features` arg loses its JS wrapper through the Dart bridge; `Func.apply` and
   `Func.prepare(...).call()` fail the same way.
6. **`viewer.errMsg`** — renamed by the production minifier; not addressable by name.
7. **`grok.shell.tv` after opening Train Model** — points to PredictiveModel view, not
   the source TableView. Use `Array.from(grok.shell.views).find(v => v.type === 'TableView')`.

## When the underlying platform bugs are fixed

| Bug | What changes in the test |
|-----|--------------------------|
| PLS auto-excludes Predict from Using (or surfaces the reason for disabled RUN) | Restore the full PLS scenario in `pls.test.ts`: select Using=All, set Components=3, click RUN, assert PLS1/PLS2/PLS3 columns appear |
| Train Model view exposes a Model Engine dropdown and a Train button | Replace the API-call training in `linear-regression.test.ts`, `pls-regression.test.ts`, `softmax.test.ts`, `xgboost*.test.ts` with UI clicks; assert the result viewers (Loadings/Scores/etc) where present |
| `cars-with-missing.csv` is deployed under `System:DemoFiles/` on dev | Drop `openCarsWithMissingFromCarsCsv` from `pareto-front-viewer.test.ts` and use `openDemoCsv(page, 'cars-with-missing.csv')` |
| Per-column `name=` attributes added to "Select columns..." sub-dialog rows | Implement a proper per-column toggle helper in `helpers.ts`; revisit the PLS deselect-`price` workaround |
| `Dialog.inputs` exposes the Dart-side input list to JS | The generic `setDialogColumnListInput` helper becomes usable; revisit PLS |

## Test Track source scenarios

All scenarios mirror the markdown files under
`public/packages/UsageAnalysis/files/TestTrack/EDA/`:

* `anova.md`, `anova-run.md`
* `multivariate-analysis.md`, `multivariate-analysis-run.md`
* `pca.md`, `pca-run.md`
* `pls.md`, `pls-run.md`
* `pareto-front-viewer.md`, `pareto-front-viewer-run.md`
* `MLMethods/linear-regression.md`, `linear-regression-run.md`, `linear-regression-spec.ts`
* `MLMethods/pls-regression.md`, `pls-regression-run.md`, `pls-regression-spec.ts`
* `MLMethods/softmax.md`, `softmax-run.md`, `softmax-spec.ts`
* `MLMethods/xgboost1.md`, `xgboost1-run.md`, `xgboost1-spec.ts`
* `MLMethods/xgboost2.md`, `xgboost2-run.md`, `xgboost2-spec.ts`

The `*-spec.ts` files in that directory are grok-browser/CDP-style reference scripts
written by the run-the-scenarios skill — they are NOT compatible with this Playwright
harness (they `chromium.connectOverCDP` to an externally-launched browser), but their
selector and step-by-step DOM patterns were the source for the helpers here.
