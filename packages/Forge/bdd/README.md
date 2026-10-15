# Forge — BDD features

The manual scenarios of the Forge Train view as features on `@datagrok-libraries/bdd`. Each feature names in its
description the scenarios it translates and what it does not claim.

| Feature | Scenarios | What it walks |
|---|---|---|
| `train/train-iris.feature` | 2.1, 6.2 (1-3), 2.5, 5.1, 5.2 (1-2), 5.13 (1-2), 2.2, 2.3, 2.8, 4.7 (1), 5.5, 5.6, 5.7 | the view on iris and its first training, the Runs grid and a restored run, a feature the method cannot use, the methods of a regression and an invalid setting, every method and the user's choice, a classifier's results and its wrong predictions, the Save dialog's Name checks, Save without data and what reached the server, the catalog row and its Details, Reference and Copy, SVM's training rows, the two balloons' links |
| `train/train-winequality.feature` | 2.4, 5.2 (4), 5.3, 5.12 (1-3, 5) | the type classifier's panels, the quality regressor and the Target/Features conflict, every regression method on 6497 rows, a change that supersedes a training, a regressor's results, a hovered point as the table's row |

Both are `@journey @serial`: one test each, the Background once, the scenarios in order, each putting back what it
changed; the features that save models run one at a time. One iris scenario is `@known-failure`: the disabled OK of
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
```

- The session signs in through the login form (`DATAGROK_LOGIN`, `DATAGROK_PASSWORD`; no dev key is needed).
  `DATAGROK_SHARING_LOGIN` names the second account, so the setup creates no `bddsecond` user.
- The server claims run the `grok s` CLI of the workspace (`datagrok-tools`) from Node, against the alias
  `FORGE_BDD_GROK_HOST` (`local` by default, `~/.grok/config.yaml`). The CLI signs in with that alias's keypair, and
  its user must be the page's: every server step compares `grok s raw GET users/current` (read once) with the page's
  current login and fails when they differ, since a user's queries read only the rows granted to that user.
- A `grok s` call costs about 0.3 s on the local stand; a full run makes 80 (25 s: iris 43, winequality 19, the
  family sweeps 18). `FORGE_BDD_GROK_LOG=<file>` appends one line per call with its time.

The stand needs: Forge published by the session's user (`grok publish local`: the debug build is the admin's), EDA
with its five methods and `knnImpute`, `iris.csv` and `winequality.csv` in `System:DemoFiles`. A full run takes about
1.9 minutes (iris about 80 s, winequality about 30-40 s).

## What a feature leaves on the server

| Left by | Cleanup, at feature start and end, read back |
|---|---|
| Saved models `forge-test-<purpose>-{time}` | `Given no Forge models named "…" are on the server`: the exact names and the family members older than an hour (a killed run), deleted in the page with `Forge:deleteModel` (the row, its model file folder, its uploaded training copy); read back by the exact names and the ids deleted (no row), the folders and the copies gone |
| Training runs of the tables `forge-test-iris-{time}`, `forge-test-wine-{time}` | `Given no Forge training runs of the table "…" are on the server`: `grok s domains delete forge.training_run --filter dataset_name = …` for the table and each stale table of the family, at most 20 deletes of 1000 each (more fails); the count of every one of them read back 0 |
| The opened tables | the shell reset closes them; nothing is saved |

Deletes are soft and a saved model is an eager entity, so each saved model leaves an entity behind (an EMS issue).
Counted as the entities of `grok s raw GET 'entities?text=forge-test'` whose name starts with `forge-test`: about 7
per full run (the 7 saved models; one `grok test` plus one full run added 67). `training_run` is eager too, so every
recorded run leaves an entity as well; that search does not count them (a run's entity is not named after its table).

## Forge's vocabulary (`bindings/`)

- `elements.ts`: the kinds `group` (a collapsible group by its caption, `aria-expanded` on its header) and `results
  panel` (a panel of the results by its title), the elements `results area` (the Train view's readiness), `train
  form`, `form splitter`, `seed copy icon`, `tags box`, `tag chips`, `column search`.
- `results-state.ts` (no steps): the results area's `aria-busy`, `aria-invalid` and `data-forge-run-id`, and the seed
  the Performance panel shows (`shownSeed`), read by the steps of `train.ts` and `server.ts`.
- `grok-s.ts` (no steps): the CLI, the identity guard, the model and run reads and the model and run sweeps.
- `train.ts`: the shown training (`user remembers the shown training`, `a new training's result should be shown`
  [within N seconds]: ready and a run id other than the remembered one, which it then remembers; `the same training's
  result should be shown`), the demo-file opener (`OpenServerFile` unprocessed, so the table keeps its creation script
  and Save offers Reference), the API fixture `a Forge model "…" is trained on the current table through the API`
  (for the catalog and apply features), the task-bar cancel and the task-bar count (every entry added since `user
  counts the task bar entries`, with every text it showed), a change while a training runs, no balloon of any colour,
  the balloons closed, an element scrolled into view, the clipboard and a Runs cell against the shown seed, two
  readings against each other, the catalog row of a model (listed, clicked) and the context panel showing it, a
  balloon's link by its model, the sizes of every results panel kept.
- `server.ts` over `grok-s.ts`: the sweeps above, `the Forge model "…" is deleted through the API` (Forge:deleteModel
  in the page, the row read back gone), and the claims `… Forge model(s) named "…" should be on the server`,
  `the Forge model "…" on the server should have:` (a dotted path into a json column, `empty`, `present` or the
  value), the fingerprint (its columns and no row values), the training copy (the table's id, rows and columns), the
  model file folder, server tables by name, the newest training run of a table, `the newest training run of the table
  "…" should be linked to the model "…"` (its `model_id` is the model's id), `the shown training should be the
  newest training run of the table "…", with:`, the run count against a remembered one (both read once the view is
  settled), the seed shown against the saved one.
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
