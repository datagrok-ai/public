# EDA Playwright tests

End-to-end coverage for Test Track scenarios under
`public/packages/UsageAnalysis/files/TestTrack/EDA/`. The ML *interaction* under test is
driven through the UI (Top-Menu navigation, analysis dialogs, RUN); **data loading and
between-test cleanup go through the JS API** (see "Shared page" below). The API is also used
where a UI path is verifiably blocked (canvas-rendered controls, Dart-side dialog state not
exposed to JS, missing test data).

Created 2026-05-18. Green 13/13 both locally against `https://dev.datagrok.ai` and on the
Jenkins **Test-Playwright** CI (the `ui_tests` stack) — no skips, no `test.fail`.

## Files

```
public/playwright-public/EDA/
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

## Shared page + API data load (why the suite is fast)

The suite reuses **one logged-in page for the whole worker** instead of a fresh browser
context per test. `helpers.ts` exports a custom `test` whose worker-scoped fixture boots the
Datagrok SPA exactly once (a single navigation to `BASE`), then overrides Playwright's
built-in `page`, so every `async ({ page }) => …` body transparently receives that shared
page — no test bodies had to change.

Within each test there is **no full-page navigation**:

* `openDemoCsv(page, name)` loads data via `grok.dapi.files.readCsv` + `grok.shell.addTableView`
  (~0.4 s). The ML menu attaches to the new TableView exactly as for a UI-opened file.
* `resetShell(page)` (in `afterEach`) calls `grok.shell.closeAll()` to drop views/tables in
  place — no reload.

This removed the two ~3 s reloads that previously dominated each test (a Files-browser
navigation to open the CSV + a home navigation to reset). Result: **~2.5 min → ~47 s**
locally (dev) and **~3.0 min → ~1.5 min** on CI — ~13 tests at ~2.6 s each (PCA ~5.7 s: it
runs two analysis passes).

Because the page is shared, the files use `test.describe.serial` (a failure skips the rest
of that file instead of running against dirtied state). Keep `closeAll()` for reset — don't
reintroduce per-test `page.goto()`.

## Running

These tests live in the standalone `public/playwright-public` Playwright project and share
its `playwright.config.ts`, `e2e/global-setup.ts` auth, and CSV reporting with the other
suites there (connections, queries, scripts, …).

### Via the grok runner (canonical)

```powershell
cd C:\DataGrok\NewBitbucket\reddata\public\playwright-public
npm install                          # first time
npx playwright install chromium      # first time
grok test --skip-puppeteer --host dev --category EDA              # all EDA tests
grok test --skip-puppeteer --host dev --category EDA --test PCA   # filter by title (regex)
```

`--category EDA` scopes the run to this subfolder; `--test` becomes Playwright's `--grep`
(test-title regex). The runner resolves the host from `~/.grok/config.yaml`, exchanges the
dev key for a token, and exports `DATAGROK_URL` + `DATAGROK_AUTH_TOKEN` to Playwright.

> The globally-installed `grok` may predate `--skip-puppeteer` (added in datagrok-tools
> 6.2.x). If it runs the Puppeteer pass instead, drive Playwright directly:

### Directly via Playwright (fast local loop)

```powershell
cd C:\DataGrok\NewBitbucket\reddata\public\playwright-public
$env:DATAGROK_URL = "https://dev.datagrok.ai"
$env:DATAGROK_AUTH_TOKEN = (Invoke-RestMethod -Method Post `
  "https://dev.datagrok.ai/api/users/login/dev/<devKey>").token
npx playwright test EDA --reporter=list                  # all EDA tests
npx playwright test EDA/anova.test.ts --reporter=list    # single file
npx playwright test EDA --grep "PCA" --reporter=list     # filter by title
npx playwright test EDA/pls.test.ts --headed             # visible browser (debug)
```

`e2e/global-setup.ts` turns `DATAGROK_URL` + `DATAGROK_AUTH_TOKEN` into browser storage
state at `e2e/.auth.json` (gitignored). Delete that file to force a fresh login.

### On CI

The Jenkins **Test-Playwright** job runs the whole `public/playwright-public` suite on a
fresh `ui_tests` stack. To exercise the EDA suite there:

* `PREREQ_PACKAGES` must include `EDA` — the ML menu (ANOVA/PCA/PLS/MVA/Train Model/Pareto)
  is provided by the `@datagrok/eda` package and won't exist on the stack otherwise.
* `TEST_GREP=EDA /` selects exactly these tests. The runner's grep is case-insensitive, so
  a bare `EDA` also matches unrelated titles (e.g. `…TestUpdateData…`); the `/` keeps it to
  the `EDA / *` describes.

A pipeline `SUCCESS` alone does **not** prove the tests ran — a `global-setup` failure
collects 0 tests yet still reports "passed". Read the real result from the build artifact
`public/playwright-public/result/playwright-public.jest` (`Running N tests`, `N passed`).

### Reports

`playwright.config.ts` uses the `list` reporter, plus a `json` reporter when
`PLAYWRIGHT_JSON_OUTPUT_NAME` is set (the grok runner sets it and converts the JSON to the
platform CSV). Don't reintroduce `expect.poll` for readiness waits: each failed poll retry
is recorded as a captured error event, so a green suite could still show errors fluctuating
run to run. These tests wait via `page.waitForFunction` (polls in-page, rejects only once
on timeout) plus a single trailing `expect`, keeping the error count at `0`.

## Last results

13 passed — **~47 s locally (dev), ~1.5 min on CI** (`ui_tests`). Per-test (CI build #56):

| # | Scenario | Result | Time |
|---|----------|--------|------|
| 1 | ANOVA on demog.csv (Box plot + Analysis + F-test tabs)                          | PASS | 2.6s |
| 2 | Linear Regression on cars.csv predicting price                                  | PASS | 2.7s |
| 3 | PLS Regression on cars.csv with 3 components predicting price                   | PASS | 2.6s |
| 4 | Softmax on iris.csv predicting Species (numeric features only)                  | PASS | 2.6s |
| 5 | XGBoost classification on iris.csv predicting Species                           | PASS | 2.7s |
| 6 | XGBoost regression on cars.csv predicting price                                 | PASS | 2.6s |
| 7 | MVA on cars.csv (Grid + 3 Scatter + 2 Bar)                                      | PASS | 3.9s |
| 8 | Pareto cars-with-missing — empty + string columns excluded from Min/Max         | PASS | 2.7s |
| 9 | Pareto cars.csv — auto-select `model` as Label                                  | PASS | 2.6s |
| 10 | Pareto demog.csv — auto-select `USUBJID` as Label                              | PASS | 3.0s |
| 11 | Pareto cars.csv — Description/Objectives/Axes/Labels/Legend categories present | PASS | 2.7s |
| 12 | PCA on cars.csv adds PC1/PC2/PC3, then PC1 (2)/PC2 (2)/PC3 (2) with Center+Scale | PASS | 5.7s |
| 13 | PLS dialog opens on cars.csv with all expected inputs and accepts Components=3  | PASS | 2.8s |

## UI vs API split per scenario

Every test loads its CSV through `openDemoCsv`, which calls `grok.dapi.files.readCsv` +
`grok.shell.addTableView` on the shared page — no navigation, no dblclick (see "Shared page"
above). The *analysis* itself stays in the UI (menu → dialog → RUN) wherever the platform
allows it; the table below lists where each scenario falls back to the API and why.

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
[pls-run.md](../../packages/UsageAnalysis/files/TestTrack/EDA/pls-run.md).

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
See [softmax-run.md](../../packages/UsageAnalysis/files/TestTrack/EDA/MLMethods/softmax-run.md).

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

See [pareto-front-viewer-run.md](../../packages/UsageAnalysis/files/TestTrack/EDA/pareto-front-viewer-run.md).
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
| `test` (custom) | Extends `@playwright/test` with a worker-scoped shared logged-in page; overrides the built-in `page` fixture. Import `test`/`expect` from `./helpers`, not `@playwright/test` |
| `BASE` | `process.env.DATAGROK_URL` |
| `inputHost(safeName)` | `[name="input-host-…"]` selector |
| `inputEditor(safeName)` | `${inputHost(safeName)} .ui-input-editor` |
| `openDemoCsv(page, fileName)` | Load CSV via `grok.dapi.files.readCsv` + `addTableView` (no navigation) |
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
| `resetShell(page)` | Escape → `grok.shell.closeAll()` (no navigation; resets the shared page in place) |
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

* Screenshot is in `test-output/<test-name>-chromium/test-failed-1.png` after a run (the
  config sets `outputDir: 'test-output'`; failures keep traces + screenshots there).
* Add an ad-hoc screenshot: `await page.screenshot({ path: 'test-output/probe.png', fullPage: true });`

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
$env:DATAGROK_AUTH_TOKEN = (Invoke-RestMethod -Method Post `
  "https://public.datagrok.ai/api/users/login/dev/<devKey>").token
del e2e\.auth.json
npx playwright test EDA --reporter=list
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
