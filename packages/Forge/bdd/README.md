# Forge — BDD features

The manual scenarios of Forge (`core/docs/features/platform/predictive-modeling/forge/manual-checks.md`) as features
on `@datagrok-libraries/bdd`. Each feature names in its description the scenarios it translates and what it does not
claim.

| Feature | Scenarios | What it walks |
|---|---|---|
| `train/train-iris.feature` | 2.1, 6.2 (1-3), 2.5, 5.1, 5.2 (1-2), 5.13 (1-2), 2.2, 2.3, 2.8, 4.7 (1), 5.5, 5.6, 5.7 | the view on iris and its first training, the Runs grid and a restored run, a feature the method cannot use, the methods of a regression and an invalid setting, every method and the user's choice, a classifier's results and its wrong predictions, the Save dialog's Name checks, Save without data and what reached the server, the catalog row and its Details, Reference and Copy, SVM's training rows, the two balloons' links |
| `train/train-winequality.feature` | 2.4, 5.2 (4), 5.3, 5.12 (1-3, 5) | the type classifier's panels, the quality regressor and the Target/Features conflict, every regression method on 6497 rows, a change that supersedes a training, a regressor's results, a hovered point as the table's row |
| `apply/apply-iris.feature` | 3.1 (3-5), 3.4, 4.10, 3.3, 3.13, 3.7 (1, 2, 4), 3.8, 3.14 (2-5), 3.9, 3.5 (1-2) | the Apply dialog's presets and the predictions against the stored Accuracy, the catalog's Apply model icon and the column's tag and Predicted by pane, close, far and repeated column names, missing values skipped and imputed, a row without a model file, Forge:applyModel and its refusal, every red row of the mapping, the old tool's menu and dialog, the trash-icon delete and the orphaned column; every application counted on the server |
| `apply/apply-winequality.feature` | 3.2, 3.6, 4.4 (3), 3.5 (1-2) | the quality model preset with every column matched, whole-number predictions with the stored MAE, a cancelled application (a batch of one row) and the model's Activity and Last run |
| `apply/apply-demog.feature` | 3.15 (2-4) | the Apply dialog on demog's natural gaps: the gapped columns and their counts, the rows without a prediction by default and the skipped count, Impute predicting every row on a copy |
| `catalog/catalog.feature` | 4.0 (3), 4.1, 4.2 (1-2), 4.3, 4.4 (1-2), 4.5, 4.6, 4.7 (2-3), 4.8, 3.7 (3), 4.11, 4.12 (1-3), 5.15, 2.6, 4.9 | the model as a platform entity, its Details, tooltip, Performance and Activity, two models compared, the row menu and the preset table, tags and a rename, Download, Applicable to, a share with the second account taken back, Run Evaluation on the stored data and on an open table, the deletes from a second catalog and from the menu |
| `app/app-and-menus.feature` | 1.1, 1.2, 6.4, 4.13 (1-2), 1.3, 1.4, 1.6, 2.9 | the app from Apps and from ML, the menus listed exactly, the Browse nodes, the Methods grid and its tooltips, the catalog's captions and icons, the old tool's Train view and Apply dialog |

All are `@journey`: one test each, the Background once, the scenarios in order, each putting back what it changed.
The six that save models are `@serial` (one at a time; the catalog, the Apply dialog and Applicable to read every
model the account sees); `app-and-menus` saves nothing and runs in the parallel lane. The apply and catalog features
train their models in the Background through `Forge:trainModel` (no data stored unless a step says so), so every claim
names the feature's own models. One iris scenario is `@known-failure`: the disabled OK of
Save model shows no reason on hover, a platform regression since 2026-10-02/03, confirmed by hand on 2026-10-15. A
probe (a plain `ui.dialog` with a `ui.setDisabled` button) found the tooltip's hover catcher under the dialog (z-index
3000) and the pointer on the dialog's command bar, while the same button outside a dialog shows its tooltip. While
tagged, the scenario passes by failing and fails the run once the tooltip shows: drop the tag once the platform is
fixed. A scenario stops at its first failing step, so its last step (CANCEL) keeps the journey clean only once the
tooltip shows; until then the harness closes the dialog with Escape after the failure.

## Running

From `packages/Forge` (Git Bash; in PowerShell set the variables with `$env:NAME = 'value'` first):

```bash
cd ../../libraries/bdd && npm run build && cd -          # only when the library reports a stale build
npx playwright install chromium                            # once per machine
DATAGROK_URL=http://localhost:8080 DATAGROK_LOGIN=admin DATAGROK_PASSWORD=admin DATAGROK_SHARING_LOGIN=forge.tester \
  npx grok-bdd run --reporter=list                         # compile --check, then Playwright
npx grok-bdd run --reporter=list generated/train/train-iris.test.ts   # one feature (the same variables)
npx grok-bdd run --reporter=list generated/apply generated/catalog    # one or more folders
```

- The session signs in through the login form (`DATAGROK_LOGIN`, `DATAGROK_PASSWORD`; no dev key is needed).
  `DATAGROK_SHARING_LOGIN` names the second account, so the setup creates no `bddsecond` user.
- The server claims run the `grok s` CLI of the workspace (`datagrok-tools`) from Node, against the alias
  `FORGE_BDD_GROK_HOST` (`local` by default, `~/.grok/config.yaml`). The CLI signs in with that alias's keypair, and
  its user must be the page's: every server step compares `grok s raw GET users/current` (read once) with the page's
  current login and fails when they differ, since a user's queries read only the rows granted to that user.
- A `grok s` call costs about 0.2 s on the local stand; a full run makes 198 (about 40 s; app-and-menus none): a model
  row is read once for the fields a model keeps, and a model sweep that finds nothing reads nothing back (a run sweep
  always reads its count back).
  `FORGE_BDD_GROK_LOG=<file>` appends one line per call with its time.
- `DATAGROK_SHARING_LOGIN` is the account `catalog.feature` shares a model with (the stand's `forge.tester`); its
  name in the Share dialog and the Sharing pane is the account's friendly name, read from the server.

The stand needs: Forge published by the session's user (`grok publish local`: the debug build is the admin's), EDA
with its five methods and `knnImpute` (Impute in `apply-iris` and `apply-demog`), PowerGrid (the Forms viewer of
Compare), the second account, and `iris.csv`, `winequality.csv` and `demog.csv` in `System:DemoFiles`. A run uses at
most two lanes, also with `--workers 4`: the `@serial` features share one, one after another, and app-and-menus runs
beside them. A full run took 2.9-3.7 minutes with `--workers 4` and 3.1-3.3 with `BDD_FRESH_PAGE=1` (2026-10-16, five
runs; the serial lane: train-iris 48-66 s, catalog 28-49 s, apply-demog 20-44 s, train-winequality 22-28 s,
apply-iris 22-25 s, apply-winequality 8-9 s; app-and-menus 16-40 s beside it). `--headed` runs two windows, and the one in
the background runs its timers slowly (apply-iris took 2.9 minutes there).

## What a feature leaves on the server

| Left by | Cleanup, at feature start and end, read back |
|---|---|
| Saved models `forge-test-<purpose>-{time}` | `Given no Forge models named "…" are on the server`: the exact names and the family members older than an hour (a killed run), deleted in the page with `Forge:deleteModel` (the row, its model file folder, its uploaded training copy); read back by the exact names and the ids deleted (no row), the folders and the copies gone |
| Models of the API fixture `a Forge model "…" is trained on the current table through the API with:` | the same sweep, registered by the step for feature end; the Background's sweep lists every name a scenario may leave (a renamed model's new name too) |
| Hand-made rows `forge-test-<x>-hand-{time}` (`grok s domains insert`, no model file) | the same sweep (the row only: it has no folder) |
| Applications of a feature's model | deleted with the model (`Forge:deleteModel`); a scenario counts them for its own model only |
| Training runs of the tables `forge-test-iris-{time}`, `forge-test-wine-{time}`, `forge-test-demog-{time}` | `Given no Forge training runs of the table(s) "…" are on the server` (comma-separated): `grok s domains delete forge.training_run --filter dataset_name = …` for each table and each stale table of their families, at most 20 deletes of 1000 each (more fails); the count of every one of them read back 0 |
| The share of `catalog.feature` (4.12) with the second account | taken back in the same scenario through the Share dialog and read back gone with `grok s shares list`; if the scenario fails in between, taken back first at feature end with `grok.dapi.permissions.revoke` in the page and read back gone with `grok.dapi.permissions.get`, both of which go to the model's wrapper project, where the grant lives (the entity itself refuses a direct grant), also after the model is deleted; only a killed run (no feature end) leaves that grant on the model's entity, which outlives the model (below), granted to the second account only |
| The opened tables | the shell reset closes them; nothing is saved |

Deletes are soft and a saved model is an eager entity, so each saved model leaves an entity behind (an EMS issue).
A full run saves 15 models (train-iris 6, train-winequality 1, apply-iris 2, apply-winequality 1, apply-demog 1,
catalog 4) and so leaves 15 such entities, plus 5 platform projects named after models (the wrappers the platform keeps a row's shares
in). The
entities of `grok s raw GET 'entities?text=forge-test'` whose name starts with `forge-test` grew by 17 in a full run
(the search does not return every one: `entities?text=forge-test-cat-quality` lists some it misses). `training_run`
is eager too, so every recorded run leaves an entity as well; that search does not count them (a run's entity is not
named after its table).

## Forge's vocabulary (`bindings/`)

- `elements.ts`: the kinds `group` (a collapsible group by its caption, `aria-expanded` on its header), `results
  panel` (a panel of the results by its title) and `detail` (a caption-value row of a key-value table of the context
  panel, Details or Sharing: `value of "Applications" detail in context panel`), the elements `results area` (the
  Train view's readiness), `forge app` (the Forge view), `train form`, `form splitter`, `seed copy icon`, `tags box`,
  `tag chips`, `column search`.
- `results-state.ts` (no steps): the results area's `aria-busy`, `aria-invalid` and `data-forge-run-id`, and the seed
  shown (`shownSeed`: the results' Performance panel's, or a model's Performance pane's in the context panel when a
  step asks for it), read by `train.ts` and `server.ts`.
- `catalog-grid.ts` (no steps): the row of the catalog grid that shows a model, from the grid's drawn cells, the
  centre of an area of that row, and what a grid draws as text (its headers' captions and its cells).
- `grok-s.ts` (no steps): the CLI, the identity guard, the model, run, application and share reads, the second
  account and the take-back of its grant, the model and run sweeps, a step's comma-separated list, a dotted path into
  a row, a model file's folder. A model row read by its name (or created by a step) is remembered: the fields a model
  keeps (id, blob, metrics, seed, fingerprint, the training copy's table id) are read from it, also after the model is deleted, while a field that
  changes (tags, description) is read anew; a sweep forgets its names.
- `train.ts`: the shown training (`user remembers the shown training`, `a new training's result should be shown`
  [within N seconds]: ready and a run id other than the remembered one, which it then remembers; `the same training's
  result should be shown`), the demo-file opener (`OpenServerFile` unprocessed, so the table keeps its creation script
  and Save offers Reference), the API fixture `a Forge model "…" is trained on the current table through the API with:`
  (target, features, method, storage, tags, description), the task bar entries counted since `user counts the task bar
  entries` (with every text each showed): their count, the end of such an entry (it has left the page, whatever it
  shows by then), a text it showed or never showed (`Rows 1-1 of 6497`), its cancel (the cross of the counted entry
  that showed the text: a running job replaces it with its progress) and a change while a training runs, no balloon of
  any colour, the balloons closed, an element scrolled into view, the clipboard and a Runs cell against the shown seed,
  two readings against each other (both read until they agree), the catalog row of a model (listed, clicked, reopened:
  another row shown first, since the context panel is built again only for a model it does not show) and the context
  panel showing it or not, a balloon's link by its model, the sizes of every results panel kept.
- `server.ts` over `grok-s.ts`: the sweeps above, `the Forge model "…" is deleted through the API` (Forge:deleteModel
  in the page, the row read back gone), and the claims `… Forge model(s) named "…" should be on the server`,
  `the Forge model "…" on the server should have:` (a dotted path into a json column, `empty`, `present`, `number` (a
  finite number, not below 0) or the value as text), the fingerprint (its columns and no row values), the training copy (the table's id, rows and columns), the
  model file folder, server tables by name, the newest training run of a table, `the newest training run of the table
  "…" should be linked to the model "…"` (its `model_id` is the model's id), `the shown training should be the
  newest training run of the table "…", with:`, the run count against a remembered one (both read once the view is
  settled), the seed shown (in the results, or in the context panel) against the saved one, a grid cell against a stored metric to three decimals (the
  Performance pane); a hand-made row without a model file, a model's applications (their count, the newest one's
  fields and its time as the start of Last run, by the id read before the delete), `the downloaded file "…" should
  have the size of the model file of the Forge model "…"` (the newest such file this worker downloaded against the
  model file on the server, which is not empty), the folder of a deleted model gone, the model
  as an entity with its author's grants, a share with the second account (`DATAGROK_SHARING_LOGIN`) present or gone,
  and a column tagged with a model's id.
- `apply.ts`: the predictions against the stored metrics (the share of rows where two columns of a named table agree
  against the train Accuracy, their mean absolute difference against the train MAE), cells emptied, a value against a
  regular expression, a column's count of missing values, the rows with a gap in any of some columns, a column empty
  exactly in those rows, a named open table without a column matching a pattern, a named table closed
  (catalog 4.7's guard: 4.6 opens a second iris copy, and when it fails before closing it, 4.11 and 5.15 would read
  it), `Forge:applyModel` from the page and its refusal's message,
  a column without a tag, the option a Dart radio input has chosen, the columns a column picker offers exactly (opened,
  read, closed with Escape), a column against the predictions of the model on the features imputed with
  `Eda:knnImpute` (a copy predicted through `Forge:applyModel`, so one more application of the model: placed after the
  application counts), and the context panel's title.
- `catalog.ts`: the catalog's loads (`user remembers the catalog loads`, `every catalog should have reloaded`, `the
  shown catalog should have reloaded`, `the first Forge view should be shown`; every open Forge view's models pane,
  a hidden one read through `grok.shell.views`, in the order of the shell's views), a model's row selected with
  Control, right-clicked, hovered (entering from the left), its cell read, the grid not listing it (by its frame,
  after the grid settled), a grid's selected rows, a grid scrolled to a column, the captions a grid draws (the header
  cells' text, `showFriendlyName` off) and a text in none of its drawn cells and headers, the open views of a name, a
  top-menu group's items exactly, and the Share dialog's second account picked (its grant on the model in the context
  panel taken back first at feature end, should the scenario fail before removing it), removed and named by the name
  the platform shows.
- `columns.ts`: the user's table has no mouse-over row, then the one a viewer reports as hovered (the Train view is
  not a table view, so the library's "the table" is not the one it trains on).

Steps that another package could use carry a `LIB-PROMOTE:` comment.

## What the platform gave these features

- The Train view publishes `aria-busy`, `aria-invalid` and `data-forge-run-id` on its results (`ARCHITECTURE.md`,
  "Automation surface"), so "a new training's result should be shown" waits for the training that a change starts and
  names its run, and "should be invalid" waits for the problem it found; the groups publish `aria-expanded`.
- The scatter plot reports `hovered row` (1-based) and `marker of row <n>` areas; the grids of the results report
  `row order` and their cells by table row (`text of cell 71 of Sepal.Length` in Wrong predictions).
- Every grid of the results reports its cells (`text of cell 1 of Metric`), the scatter, bar and PC plots their rows,
  bars and axes; the Runs grid its `cell <r> of #` areas, clicked after the panel is scrolled into view.
- A red columns input marks its editor (`d4-invalid` on the `<div>` editor and the label, no `aria-invalid`), so the
  claim is `editor of Features input should be invalid`.
- The "Select columns..." dialog's grid reports only the rows it draws: a column below them is typed into its Search
  box first.
- The Forge view carries `forge-app`, and its models pane `.forge-catalog` publishes `aria-busy` and
  `data-forge-catalog-load` (`ARCHITECTURE.md`, "Automation surface"): a delete or an edit is claimed by a later load
  in every open catalog, a Refresh by a later load in the shown one, before the grid is read, so a row's absence is not
  read off the old frame.
- The catalog grid reports `row header <n>` and `cell <n> of <column>` areas by table row, and `header <column>`, where
  the grid's hit test finds the header cell whose text is the drawn caption (`{widget} should show the columns`). A
  grid in a narrow pane (Activity) reports only the
  columns it draws: the feature scrolls to a column first. A pointer that rests in a cell asks for no tooltip.
- A red column input of the Apply dialog marks its editor and label (`d4-invalid`), its Model choice the `<select>`; the
  summary of Columns carries `aria-invalid`; the disabled OK carries `aria-disabled` but shows no reason on hover in a
  dialog (a platform regression).
- A Dart radio input marks its option only by `checked` (no value attribute, every option reads "on").
- A task-bar entry shows its job's text first and then its progress ("Rows 1-1 of 6497").
- A vertical top-menu group holds its items in `.d4-menu-item-container-fixed`; the ML menu is shown over a table view
  only.
- The Share dialog's user typeahead finds an account by the start of a word and shows its friendly name ("Forge
  Tester"), as do the dialog's grant rows and the Sharing pane; `grok s shares list <id>` lists the grant made through
  the dialog, so it is the claim of record.
