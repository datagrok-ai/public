# u2 agent evaluation

UI work an agent is asked to do in a plugin with u2: twelve single jobs (`t01`…`t12`) and three whole features
(`f01`…`f03`). Each is built around misuse that really happened, and in each the trap is incidental: the agent gets
a product request, never a word about the hazard, the checks or being evaluated. Each task's machine grader checks
the behaviour asked for and the absence of the misuse. The pass rate of fresh agents' solutions is u2's
**AI-friendliness score**.

## The tasks

This table is the private key. It never goes to an agent, and neither does anything else under `eval/`.

| Task | Entry | The job | The real misuse |
|---|---|---|---|
| `t01` | `submit-dialog.ts` | a dialog of read-only rows, with pencil-to-edit and download icons | B-1: dead rail icons on read-only rows; Enter/Escape reaching the dialog |
| `t02` | `compound-summary.ts` | a summary form whose Structure row is a molecule drawing | N-1: a 160×100 drawing on the 28px options rail |
| `t03` | `merge-dialog.ts` | a file field that also takes a Datagrok Files file or path | B-2: "pick from Files" lost |
| `t04` | `molecule-column-input.ts` | a picker for a table's molecule column | B-4: the hidden `~` column offered and defaulted |
| `t05` | `settings-dialog.ts` | a function's settings form | B-3/N-2: crash with a nameless table open |
| `t06` | `campaigns-panel.ts` | a panel that builds its rows after `await load()` | unowned builds after `await`; builds into a disposed panel |
| `t07` | `new-template-form.ts` | a SAVE that is off until the form is right | reasonless (boolean) or "Disabled" reasons |
| `t08` | `campaigns-table.ts` | a grouped table with row hover | B-6: hovered row = group header colour |
| `t09` | `campaign-ribbon.ts` | TableView ribbon actions, one adding a `DG.Row` | C-7: a dropped promise; N-3 (a `DG.Row` returned to the button runner) is absorbed by u2, which takes a value a `run` answers as a plain answer |
| `t10` | `campaign-catalog.ts` | a catalog grouped by a date and by an object | N-5: group header ≠ cell text; `[object Object]` groups |
| `t11` | `template-editor.ts` | an editor laid out as designed | package CSS reaching into `.u2-*`, `!important`, inline style hacks |
| `t12` | `current-view-badge.ts` | a badge following `onCurrentViewChanged` | a platform subscription nobody ends |
| `f01` | `campaign-actions.ts` | a campaign's ribbon, Submit and Merge dialogs in a TableView of a busy session | t01, t03, t04, t09, and a dependent input built after an `await` |
| `f02` | `compute-settings.ts` | a compute settings pane that loads, lists functions with their parameters, and saves | t05 (here the nameless table breaks the whole pane, as B-3 did), t06, t07 (a reasonless SAVE and an inert pane), t11 |
| `f03` | `campaigns-home.ts` | an app home: a loading, grouped catalog with a details pane | t02, t06, t08, t10, t11 |

`task.json` keeps the rest of the key: `name` (what the task was called before the ids went neutral), `source`
(the bug), `hazards` (for a feature task, the single tasks it combines), `bad` (the checks each bad solution
must fail) and `absorbed` (bad solutions of a bug u2 now contains, which must fail nothing, each with why: if one
starts failing again, the guard that absorbed it is gone).

## Keeping the trap incidental

- **Neutral ids and entries.** Folders are `tNN`/`fNN`; the file an agent writes is named for the feature
  (`merge-dialog.ts`), never for the hazard.
- **Prompts are product requests.** A `task.md` says what the feature does and the API surface it is built on (the
  function the plugin calls, the platform calls the feature makes), as a developer gets it. It never names the
  hazard or what the grader looks at ("must not get lost", "whatever tables happen to be open").
- **One environment for every task.** What the test page offers is listed once, in `AGENT-PROMPT.md`, the same for
  every task, so no list singles out a trap's ingredient.
- **The hazard lives in the fixtures.** The graders' data is what users have: a campaign table read back with
  Chem's `~Molecule.canonicalSmiles` first, a session with a nameless table and another view current
  (`lib/platform.mjs` `busySession`), a refresh that fails, a pane closed while it loads.
- **Features, not puzzles.** `f01`…`f03` are whole features in which several hazards meet the way they met in
  HitTriage, none of them the point of the request.

## Layout

```
eval/
  AGENT-PROMPT.md  the instruction an agent gets, with the task's task.md
  tasks/<id>/
    task.md        the product request
    task.json      the private side (above), and `entry`: the file the agent writes
    grader.mjs     the checks: {id: async (mod, ctx) => …}, failing with a message
    good.ts        a reference solution (+ good.css where it has one)
    bad-*.ts       solutions reproducing the real bugs (+ their .css)
  lib/             the grading substrate (below)
  run-graders.mjs  proves the graders valid; writes results/baseline-graders.json
  grade-run.mjs    scores a folder of agent outputs
  results/         results a grading run writes (local, git-ignored)
  runs/<date>/     agent outputs: <id>/<k>/ per sample, RUN.md, score.json (local, git-ignored)
```

## How a solution is graded

Each solution runs in its own node process on u2's test substrate: `tests/dom-shim.js`, the platform stub
`tests/dg-stub.mjs` and the doubles in `tests/platform-doubles.mjs`, as they are. `lib/platform.mjs` adds what the
doubles lack: `columns.bySemType`, `col`, `rows.addNew` answering a `DG.Row`-like proxy, every `on…` event js-api's
DataFrame declares, a TableView with a ribbon, `DG.ViewBase` / `DG.View` / `DG.TableView` classes that `instanceof`
the page's views (served through the page's `datagrok-api/dg`, `lib/hooks.mjs`), and the busy session;
`lib/harness.mjs` stands in for `grok.chem.drawMolecule` and `ui.pickFileFromFiles`.
`@datagrok-libraries/u2` resolves to u2's test build (`.test-build`, compiled by `npm test`), so **run `npm test`
in `libraries/u2` first**; the graders refuse a build older than its sources or CSS. A `.ts` solution is
transpiled the way u2 compiles its own code, and a `.css` import is installed in the page as package CSS, after
u2's sheet.

Graders act as the user does and read only what u2 exposes publicly. A click is refused when the element couldn't
be operated: hidden, inert, disabled, or `pointer-events: none`. A choice is picked by the text the user reads.
Reads go through `Control.forElement(…)` `.getState()` and `getItems()`, `data-u2-part` / `data-u2-item`, ARIA, the
DOM, and `lib/css.mjs`. That module is a small cascade over the page's sheets, so it can see what the shim cannot:
hover colours, `pointer-events`, widths, gaps. Every solution also gets two checks: **loads** (it imports) and
**typecheck** (one TypeScript program over all the solutions, typed against u2's sources).

A solution scores the fraction of its task's checks it passed. u2 may change, and the graders keep working as long
as the behaviour holds. The reference solutions use today's API, though: when the API changes, they need the
codemod like any other consumer.

## Grader validity

```
npm run eval -- --repeat 3         # = node eval/run-graders.mjs --repeat 3; one task: node eval/run-graders.mjs f02
```

First `lib/promises.mjs` looks up, in the page, every platform API `AGENT-PROMPT.md` promises (the DataFrame
events read off js-api's own declaration); a missing one makes the run invalid, since an agent that used it would
fail for the page, not for u2. Then every `good.ts` must pass all of its checks, and every `bad-*.ts` must fail
**exactly** the checks its `task.json` lists (an `absorbed` one: none), in every round; a check whose verdict
changes between rounds is reported as flaky. A bad solution that
failed one check fewer or more would be failing for some other reason. The command exits 1 otherwise. A run over
all the tasks writes `results/baseline-graders.json`: the u2 version, each round's verdict, each check's id, and
each bad solution's failures with their messages. The harness is not part of u2's `npm test` or CI. It grades u2
strict, as u2's node runner tests it (`U2_STRICT=1` unless set): a misuse throws, so a bad solution's ownership bug
fails every check that needs what it never built; a check calls the solution (`need(mod, name)`) under a scope of
its own, as a page would.

## Running an agent: the protocol

1. **A fresh agent per task and sample.** It starts with no conversation or memory and gets
   [`AGENT-PROMPT.md`](AGENT-PROMPT.md) with its slots filled and the task's `task.md` pasted in: a neutral "you
   are adding this to a Datagrok plugin" instruction, what it may read (u2's docs and source, the JS API source),
   and the test page's environment. It writes into an empty folder outside the repo, so neither `eval` nor a task
   id is in any path it sees. The plan asks for 3 samples per task.
2. **Copy each output folder to `eval/runs/<date>/<id>/<k>/`**, `k` = 1, 2, 3…, unchanged (the entry file, any
   other `.ts`, the `.css`). Write down the model, the doc set and the prompt's revision in
   `eval/runs/<date>/RUN.md`.
3. **Grade it:** `node eval/grade-run.mjs eval/runs/<date> [label]`. This prints each task's score, pass@1 and
   failed checks, writes `eval/runs/<date>/score.json` and, given a label, `results/<date>-<label>.json`.
4. **The numbers.** A sample scores its fraction of checks passed. A task scores the mean over its samples, and its
   **pass@1** is the share of its samples that passed every check. A set (`single` = t01–t12, `feature` = f01–f03,
   `all`) scores the mean over its tasks; a task with no output counts 0. Sample *i* of every task makes the run
   one agent per task would have scored: the min–max and sd of those runs is the **spread**, how far one run's
   number can sit from the mean. The plan's target is a mean of ≥ 0.9.

Runs made before the ids went neutral (2026-10-27: `<NN-name>/solution.ts`) are still read: `grade-run.mjs` maps
`04-hidden-column` to `t04` through `task.json`'s `name`, and takes `solution.ts` where there is no entry file.

## Adding a task

A newly found misuse becomes a task, or a hazard in a feature task. Take the next neutral id (`t13`, `f04`). Write
`task.md` as the product request a developer would get, keeping the hazard in the fixtures and out of the text;
extend `AGENT-PROMPT.md`'s environment if the feature needs a platform call it does not list. Write `grader.mjs`
with checks that act and read through public surfaces. Add `good.ts`, at least one `bad-*.ts` reproducing the real
bug, and `task.json` with `name`, `entry`, `title`, `source` and the checks each bad solution must fail. Then
`npm run eval -- --repeat 3` must pass.
