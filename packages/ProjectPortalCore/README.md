# Project Portal Core

Several portfolios on one Datagrok instance. One organization can run small-molecule
discovery, bioprotein discovery and bioprotein CMC side by side. Each portfolio has its own
program registry feed, its own stage-gate model, and its own access groups.

## The model

**The portfolio is the unit of configuration.** Everything that varies between portfolios
hangs off the `portfolio` row. The alternative, one profile per instance, breaks as soon as
one customer runs discovery and CMC together.

| Table | What it holds | Security |
|---|---|---|
| `stage_model` | A named, versioned stage-gate process. Portfolios may share one. | table, all users view |
| `stage` | A stage or gate of a model. Keyed by `(stage_model_id, code)`, so `LO` exists in both discovery models. | table |
| `gate_criterion` | What must hold to enter a stage: `document`, `data_condition` or `decision`. A data condition may carry a structured `condition` the portal evaluates. | table |
| `portfolio` | Setting, stage model, **stage authority** (`feed` or `portal`), groups, program noun. | table |
| `registry_feed` | SQL against a named connection that yields the portfolio's programs. | table |
| `program` | Keyed by `(portfolio_id, code)`. `predecessor_id` links a program that continues in another portfolio. | **row, private by default** |
| `stage_transition` | The gate log: every stage change, including `sync` changes observed in a registry. Written once. | master → program |
| `compound` | A projection of the compound registry, keyed `(system, code)`. Never minted, never promoted, never a folder. | table |
| `program_compound` | A compound in a program, with its role: explored, hit, lead, backup, candidate, tool, discontinued. Many-to-many. | master → program |
| `compound_role_change` | Every role change: who, when, why, from the feed or by hand. Written once. | master → program |
| `campaign`, `series`, `study`, `process_stage` | Work units under a program, keyed `(program_id, code)`. A campaign is open-ended and reaches an `outcome`; a study is closed-ended and completes. | master → program |
| `criterion_status` | Met / waived / not met per program and criterion. | master → program |

### Rules that make several portfolios work

- **Stage authority.**
  - Where the registry owns the stage (typically development), the portal mirrors it. Every change it observes is logged as a `sync` transition, and the portal refuses to record decisions.
  - Where the portal owns the stage (discovery), decisions are made here. Programs advance one stage at a time, and only when the next stage's required criteria are met or waived.
  - The registry's own spelling is always kept in `source_stage`.
- **One source, several portfolios.**
  - Two portfolios may read the same registry; each feed's SQL takes its own slice.
  - The **routing audit** (`portfolioRoutingAudit`) runs each source's full key list. It reports keys that no portfolio claims and keys that two portfolios claim.
- **Lineage, not moves.**
  - A molecule that leaves a discovery portfolio and appears in a development feed becomes a new development program, with `predecessor_id` pointing to the discovery program.
  - The discovery program is marked `transitioned`.
  - If both feeds deliver the key at the same time, the sync reports an overlap instead.
- **Missing is not deleted.** A program that drops out of its feed is flagged `missing_from_feed`, never removed.
- **Fail closed.**
  - `program` rows are private to their author until they are shared.
  - Sync shares each new program with its portfolio's groups: approvers ⊂ contributors ⊂ viewers, with View to viewers and Edit to contributors.
- **A workspace per portfolio.**
  - `buildWorkspace` creates a Spaces tree: the portfolio's `workspaceRoot`, a folder per program named `CODE · Name`, and the portfolio's `workstreams` under each program.
  - Compounds, plates, results and predictions are never folders.
  - It is additive: it creates what is missing and never renames, moves or deletes. It records each program's folder on the program row (`space_id`), so a program renamed in its registry keeps its folder, and the report names the difference.
  - The root is shared View to viewers and Edit to contributors, and its folders inherit that. Unlike program rows, spaces accept the grant, so the workspace is isolated per portfolio today.
  - **It grows with the stage (PORT-35).** A stage declares `folders`, paths under the program's folder (`Candidate Package`, `DMPK/In vivo PK`). A program has the folders of every stage it has reached: the builder adds them, and a gate decision that advances a program adds the new stage's folders at once when its workspace exists. A program that recycles keeps them.
  - **Orphans are reported; the builder never touches them (FR-22).** Every full build, and `workspaceOrphans`, lists root folders no program claims and the folders of programs that are terminated, transitioned, completed or missing from their feed.
  - **Parking is a separate, explicit step.** `parkOrphans(code, dryRun, folders)` moves the listed folders (or the ones named, by path or space id) into `_Parked/<date>` under the root. A dry run is the default. Applying needs one of the portfolio's approvers. A move keeps the folder's id, files, subfolders and grants. Each move is recorded in `parked_folder` (from where, to where, why, who, when). `unparkFolder(id)` moves one back, by the record's id or the folder's, unless a folder of the same name took its place; an emptied `_Parked/<date>`, then `_Parked`, is removed. `parkedFolders(code, all)` lists them. `_Parked` is never itself an orphan, and the builder neither rebuilds a parked program folder nor makes a new one ("folder parked; not rebuilt").
  - **Deleting is out of scope.** Nothing here deletes a workspace folder; removing a parked folder is a person's call, by hand.
  - **In the app.** The ribbon's **Workspace** (Advanced mode) shows the orphan report with **Preview park**, **Park** on each folder and **Park all** for approvers, and the parked folders with **Unpark**.

### Work units

A portfolio enables the work units it uses (`workUnits`): discovery typically campaign and series,
development study and process stage. `workUnitFeeds` deliver them from a source system, one feed per
kind, each row naming its `program` by code. They sync after the programs, by the same rules:
- keyed by program and code, so a re-run changes nothing;
- the kind's own columns are mapped and checked against their choices; anything else lands in `details`;
- a row for a program the portfolio does not hold is rejected and reported;
- a unit its feed stops delivering is flagged `missing_from_feed`, never deleted; a unit added by hand is left alone;
- series maturity is checked against the portfolio's `seriesMaturity`, and drift is warned.

A work unit is never a folder of compounds (FR-26). The program page lists the units by kind.

### Compounds in programs

A portfolio's `compoundFeed` (registry `system`, `connection`, `sql`) returns one row per program and
compound: `program`, `compound` (the registration code), `role`, and `name`, `smiles`, `aliases`,
`batches` when present. Sync runs after the programs and work units:
- one `compound` row per registry code, however many programs hold it; a tool compound serves several programs with a different role in each;
- the registry moves a role, and the change is logged in `compound_role_change`;
- `setCompoundRole` changes a role by hand, with a rationale; that link becomes manual, and later syncs report the registry's role without overwriting it (the `program_flow` rule);
- a link the feed drops is flagged `missing_from_feed`, never deleted.

Which programs hold a compound is under each program's security; the compound rows themselves are
readable by all users, as a registry usually is. Rows are written one by one, which suits a
portfolio's lead and candidate sets, not a full 10⁶-compound registry (NFR-01): bulk binding waits
on EMS over external databases (PLAT-13).

### Data conditions and readiness

A data condition is evaluated, not ticked. The portfolio declares a **program dataset**
(`programData`): one SQL on a named connection, with `@program` standing for the program's code,
returning one row per compound, lot or whatever unit the criteria count. A criterion's `condition`
is a declarative predicate over those rows, the one rule language of 07 §2.10:

```json
{"scope": "lead_set", "measured": ["solubility_um", "papp", "clint"], "minCount": 10}
{"scope": "series", "where": [{"column": "series_maturity", "op": "in", "value": ["validated", "lead"]}], "minCount": 3}
{"where": [{"column": "purpose", "op": "=", "value": "tox"}, {"column": "released", "op": "=", "value": true}]}
```

- **Scope.** `program` counts every row. `lead_set` counts rows whose `role` is a lead role (`leadRoles`,
  default lead, backup, candidate). `series` groups rows by `series`, and needs `minSeries` series
  (default 1) to reach `minCount` each.
- **Operators.** `<` `<=` `>` `>=` `=` `!=` `in` `not_in` `present` `absent`. A blank never compares.
- **Measured.** A row qualifies only with a value in every measured column; the rows that miss one
  are named in the gap list.
- **Readiness** (`programReadiness`) evaluates every criterion of a stage: data conditions from the
  data, documents and decisions from `criterion_status`. A recorded waiver wins; a recorded "met"
  does not outvote the data. It returns met over required, a score, and the gap list.
- **The gate uses it.** Advancing needs the evaluated conditions met (or waived). On advance the
  portal writes what the data said into `criterion_status.evidence`, beside the decision.
- A data condition with only a free-text `rule` is still recorded by hand.

### The portfolio table

The Table button on the board opens the portfolio as a Datagrok table view, at
`/apps/ProjectPortalCore/Portal/<code>/table`. All portfolios opens every portfolio you can see as
one table, at `/apps/ProjectPortalCore/Portal/table`; there, the portfolio column is a filter.
- **Rows.** One per program: stage (a category in the model's order), days in stage, next stage, required criteria met over total, a readiness score and a ready flag, compounds by role, work units by kind, whether a workspace exists. Unstaged programs read `Unstaged`.
- **Days in stage** count from the last transition into the current stage (`stage_transitions`), or from the program's creation when it never moved.
- **Readiness** is the gate panel's own computation (`readinessFor`).
- **Fast to read.** The data is read in a fixed number of batched queries. Only the program dataset is read per program, a few at a time, because its SQL selects by `@program`.
- **Filters and viewers.** The view opens with the filter panel, programs by stage and modality, days in stage, the readiness matrix (programs × the criteria of their next stage, state as color) and compounds by role.
  - The matrix and the compound roles are tables of their own, linked to the programs by id, so the filters drive every viewer.
- **Clicking a row** opens the program's page beside the grid, the same page the board shows.
- `portfolioTable(portfolio)` returns the same programs table to scripts.

### Model runs

A model run is saved to a program with `saveModelRun(call, portfolio, program, ...)`. The package
ships one model to try it on: `FedBatchTiter` (`scripts/fed-batch-titer.js`, `meta.role: model`),
a synthetic fed-batch titer model for EX-BIOC. It takes feed rate, temperature, duration and seeding
density, and returns the final titer, peak cell density, final viability and a daily time course.
The numbers are invented for the example.

Saving a run does four things:
- **Annotates the call.** The call's options carry `portalPortfolio`, `portalProgram`, `portalStage`
  (the stage at run), `portalUnitOperation`, `portalRunType`, `portalPurpose`, and the `title`,
  `description` and `tags` the Compute run history shows. Options persist with the call and can be
  queried (`options.portalProgram = "EXB-203"`).
- **Saves it.** Dataframe inputs and outputs are uploaded first and shared with the portfolio's
  viewers. The call itself has no ACL.
- **Links it into the workspace.** The run is linked, not moved, into the folder the program's
  current stage declares (a run at D5 goes to `Drug Substance/Process Lock`), into a folder named in
  `folder`, or into the program folder.
- **Records a `model_run`.** The row carries the stage at run, unit operation, run type, purpose
  and the catalogued `model` (by the function's nqName). A gate criterion selects it with
  `selectRunForCriterion`: the run becomes the `selected` artifact, an earlier one `superseded`,
  and the criterion is met unless it was waived.

### Model catalogue and model cards

A portfolio file declares its model catalogue under `catalogue` (F-12, FR-50):

```json
"catalogue": {
  "unitOperations": [{"code": "EX-PROD", "name": "Production bioreactor", "ordinal": 2}],
  "runTypes": [{"code": "EX-SIM", "name": "Simulated", "definition": "..."}],
  "models": [{"code": "ProjectPortalCore:FedBatchTiter", "name": "Fed-batch titer (example)",
    "unitOperations": ["EX-PROD"], "defaultRunType": "EX-SIM",
    "card": {"version": "1.0", "question": "What final titer ...?", "intendedUse": "...",
      "inputs": [{"name": "feedRate", "units": "% v/v per day", "min": 0, "max": 6}],
      "applicability": "...", "performance": "...", "lineage": "...", "owner": "...", "approval": "draft"}}]
}
```

- **A model is its function.** Its code is the function's nqName (`Package:Name`), the handle that survives a release.
- **Checked with the set.** A model may name a unit operation or run type that any portfolio of the set declares. A code two portfolios declare must be declared the same way; a model two portfolios declare lists both.
- **Loaded with the portfolios, idempotently.** Rows are matched by code. A model's portfolios and unit operations only grow, so loading one portfolio never drops another's claim. Runs saved before the model was catalogued are linked to it.
- **One card per model version** (`model_card`); `model.current_version` names the current one, and older cards stay as lineage.
- **The card is written around the question** (PORT-41): question, intended use, inputs with their validity range, outputs, applicability, performance, lineage, limitations, owner, approval (`draft`, `in_review`, `approved`, `retired`). Inputs and outputs the card does not list come from the function's signature; a listed one takes what it leaves out (caption, units) from the function.
- **Where it shows.** The program page lists the program's saved runs; a run's model opens its card. `modelCard(model)` returns the card as a widget.
- A package with a catalogue but no portfolio hands it over with `loadModelCatalogue(catalogue, portfolio)`.

The example catalogue (EX-BIOC) declares five unit operations, three run types and `FedBatchTiter`
with a draft card. The codes carry the `EX-` prefix because unit operations and run types are global.

### Model recommendation

A card may say where its model applies (F-13, FR-51). Every key is optional; a key left out matches anything:

```json
"appliesTo": {"tasks": ["physchem"], "stages": ["H2L", "LO"], "modalities": ["Small Molecule"],
  "workUnits": ["series"], "predicts": ["solubility_um"], "rank": 1}
```

- **Missing assays.** For the gate the program page shows (the next stage, or the stage clicked on the track), each open data condition lists the columns it needs that rows lack: `measured` columns on the rows that pass `where`, and columns compared with a number in `where` on the rows in scope. A column matched against a value or list (`series_maturity in [...]`) classifies a row and is not an assay. A criterion with only a rule in words is listed without a column. With a series, only its rows count.
- **Suggested models.** Candidates are the models the portfolio's catalogue declares, and any model whose card declares `appliesTo`. A deprecated model or a retired card is never suggested. Each declared key must match: stages (the program's stage or the one it works towards), modality, work-unit kind (the portfolio's, or `series` for a series), and the task when one is asked for. Order: models predicting a missing column, then the more keys matched, then approved, in review, draft, then `rank`, then name.
- **Program page.** *Missing assays for <stage>* says which rows lack what, and which suggested model predicts it meanwhile. *Suggested models* lists each with why; the name opens its card, **Run** opens its dialog filled in from the program: a dataframe input gets the program dataset (opened as a table), inputs named `program`, `portfolio` or `series` get those codes.
- `recommendModels(portfolio, program, stage, series, task)` returns one row per suggested model (model, name, rank, predicts, matched, reasons, approval, installed). `missingAssays(portfolio, program, stage, series)` returns one row per missing column.

EX-SMD declares two synthetic models, `SolubilityEstimate` and `HergRisk`: tiny scripts that compute an
invented number from potency. They exist to show the flow; their cards say so.

### Predictions

What a catalogued predictive model said about each compound of a discovery program, kept per model
version (F-22, FR-58/103).

- **The record.** A `prediction` row: program, compound (code, and the registry compound when the code is one), model (catalogue code), model version (the card's version), the program-data column it stands in for (`predicts`), the value, units (read off the column name: `_um` is uM), and the scoring run (`model_run_id`, and its `run_id`). Master-secured through the program. Removing the run removes its predictions.
- **Batch scoring.** `scoreCompounds(portfolio, program, model)` runs the model once over the program dataset, its dialog's own prefill (the table, the program's codes). The card's `appliesTo.predicts` names the columns; the output column holding each is `predicted_<col>`, `estimated_<col>`, `<col>_predicted`, `<col>_pred`, `<col>_estimate` or `pred_<col>`, never the measured copy. The call is saved with `saveModelRun` (purpose screening, title `Scoring · <model> <version> · <program>`), so provenance reads prediction → run → model.
- **Idempotent per version.** When every compound already has a prediction from this model version, nothing runs. Otherwise one run, and only the compounds not yet scored are stored; what is stored is never rewritten. A new card version scores again.
- **Predicted vs measured.** Per model and column, the latest prediction per compound joined with the program dataset's measured value of the same column: how many were measured, the mean absolute error, Pearson's r.
- **Program page** (discovery portfolios, or wherever predictions exist; behind More in Simple mode). **Score** beside **Run** on a suggested model whose card names what it predicts. *Predictions (n)*: per model and column the fit, the first 12 compounds predicted and measured, **Predicted vs measured** (a table view with a scatter plot), and Provenance of the scoring run. Each compound listed shows its predictions. A missing assay shows the stored predictions of the compounds lacking it ("Predicted: EXC-1006 34.1 uM (...)") before the models that could predict it; `missingAssays` has a `predicted` column.
- **Elsewhere.** The activity feed has a `prediction` entry per scoring run; provenance has a *Predictions* node under the scoring run.
- **Functions.** `scoreCompounds(portfolio, program, model)` and `programPredictions(portfolio, program, model)` (compound, predicts, model, version, predicted, measured, units, model run).

### Program data

Each portfolio declares its program dataset, `programData`: one SQL on a named connection, with
`@program` for the program's code (F-09). Data conditions are evaluated over it, and the program
page opens it with **Open data** as a table view named after the program (`EXS-KIN1 · Kinase
inhibitor`). A portfolio may name it and further datasets beside it:

```json
"programData": {"name": "SAR table", "connection": "System:Datagrok", "sql": "... WHERE program = @program"},
"programDatasets": [{"code": "assays", "name": "Assay results", "connection": "System:Datagrok",
  "sql": "... WHERE program = @program",
  "viewers": [{"type": "Box plot", "options": {"categoryColumnNames": ["assay"], "valueColumnName": "value"}}]}]
```

- **One button per dataset:** *Open data* for the program dataset, *Open Assay results* for a further one.
- **Viewers** are the platform's, with its property names. Without any, the defaults follow the portfolio's setting. Discovery: a scatter plot of the first two numeric columns, coloured by the first category, and a histogram. Development: a bar chart over the first category, and a histogram. The program column is never used. A viewer naming a column the table lacks is skipped.
- **Governed:** the query runs as the user on a connection shared with them, and the program is read under row security first.
- `programData(portfolio, program, dataset)` returns the table; `openProgramData(portfolio, program, dataset)` opens it.

### Compounds needing attention

A portfolio may declare `attention` rules over its program dataset (F-17). A rule uses the
condition language of the criteria (scope, `where`, `measured`), but it flags rows instead of
counting them:

```json
"attention": [
  {"code": "herg_unmeasured", "name": "Potent, hERG not measured", "severity": "warning",
   "condition": {"where": [{"column": "ic50_um", "op": "<", "value": 0.5}], "measured": ["herg_ic50_um"]},
   "message": "{compound}: IC50 {ic50_um} uM, no hERG IC50"},
  {"code": "lead_no_solubility", "name": "Lead without solubility",
   "condition": {"scope": "lead_set", "measured": ["solubility_um"]}}
]
```

- **What trips a rule.** A row trips it when the scope admits it and every `where` predicate holds. If the rule names `measured` columns, the row must also lack at least one of them. `minCount` and `minSeries` do not apply.
- **Severity** is `critical`, `warning` (the default) or `info`. **Message** fills `{column}` with the row's value and `{missing}` with what it lacks. Without a message, a line names the rule, the values it tested and what is missing.
- **An item** is the row's `compound`, `id`, `code`, `lot` or `name`. It counts once, however many rules it trips.
- **Program page.** *Needs attention (n)* shows a count per rule that trips, then each flagged item with its rules and messages. A click opens the program dataset with that row current and every flagged row selected.
- **Portfolio table.** An `attention` column (flagged items, a filter) and `attention rules` (the count per rule). Both are hidden for a portfolio without rules. A program whose data cannot be read has no count.
- `programAttention(portfolio, program)` returns one row per item and rule: item, rule, rule name, severity, message.
- **Examples.** EX-SMD: potent but hERG not measured, a lead without solubility, hERG IC50 under 10 µM. EX-BIOC: a tox lot not released.

### Scheduled sync

A sync runs on call: the board's **Sync** button, `syncPortfolio`, `syncAllPortfolios`. A portfolio
may also declare a `schedule` and `owners` (F-01):

```json
"schedule": "every 6h",
"owners": ["admin"]
```

- **Schedule.** An interval (`every 6h`, `30m`, `1d`, at least 15 minutes) or a five-field cron line in UTC (`0 6 * * 1-5`, `@daily`). When both day fields are set, either matching is enough, as in cron.
- **`syncDuePortfolios(codes)`** syncs the portfolios whose schedule is due, discovery first. A portfolio without a schedule is never touched; `codes` (comma-separated) narrows the check further. A sync is due when the portfolio has never synced, or the next slot after the last run's start has come. Any run counts, by hand or failed.
- **Sync history.** Every sync, by schedule or by hand, is a `sync_run` row: started, finished, trigger, who, counts of created, updated and flagged rows, unresolved stage spellings, error, the report, and the per-program changes. A dry run is not recorded, and a sync never fails because its history could not be written. The board's **Sync history** shows the schedule, the next slot, the owners and the last 20 runs. `syncHistory(portfolio)` returns the runs as a table.
- **Owners are notified** in the platform's notifications (the bell) after a scheduled run that changed or flagged something, failed, or found different unresolved spellings or overlaps than the run before. The server lets a user notify another user only with the Edit users permission; otherwise the run records who could not be reached.
- **Per-program changes** (new, fields changed, stage moved, no longer in the feed, work units, compounds added, role changes, advanced by data) are on the result and in `sync_run.changes`.
- **New data for followers.** After the feeds, every sync checks the followed programs for new data (see below).
- **Nothing runs by itself yet.** The platform schedules only server-side scripts and queries (`#schedule:`), not a package's JavaScript functions, so the package cannot register its own schedule. Something must call `syncDuePortfolios` (see Known platform issues).
- The examples declare schedules: EX-SMD every 6 h, EX-BIOD daily at 06:00 UTC, EX-BIOC weekdays at 06:00 UTC. Owner: `admin`.

### Following a program: new-data notification

A user follows a program, or one of its series, and hears when new results land (F-11, FR-83).

- **Follows** are rows of `follow`: user, program, optional work unit (kind and code), and two snapshots of the program dataset: as last checked for news, and as at the user's last visit. Each user writes their own. Unfollowing turns the row off (`active`); following again turns it back on with a fresh baseline.
- **The snapshot** is cheap and conservative: rows, filled cells (so a value landing in a blank cell of a wide SAR table counts), and, when `programData.watermark` names a column that only grows (a run date, an increasing id), its largest value. With a watermark, rows past the last value are new, so a replaced row counts; without one, only growth in the counts is news. Rows or values that disappear, or a value overwritten by another, are never reported. A series follow looks only at rows whose `series` column is that series.
- **When.** After every sync, by hand or by schedule, and on `checkNewData(code)`. Only followed programs are read, one query each. Following takes the data as it stands as the baseline, so what is already there is never news.
- **Telling.** Each follower with news gets one notification per portfolio in the platform's notifications (the bell), the same route the sync uses for owners: `EXS-KIN1: 1 new row, 4 new values`. Nothing external. A sync run by a user without Edit users can notify only themselves.
- **Program page.** **Follow**/**Unfollow** beside Open workspace, and Follow/Unfollow on each series row. A follower sees *What's new since your last visit* once per open of the program; the visit is then recorded.
- `followProgram(portfolio, program, series)` and `unfollowProgram(...)` do the same by code.

### Pinned programs: My programs on Home

Pinning a program pins its workspace folder. The pin is a platform favorite on the folder, the user's own or a
group's; the portal stores nothing. A pin made or removed anywhere (the program page, the star, the Workspace) is
the same pin.

- **Program page.** **Pin to my Workspace** / **Unpin** beside Follow, in the app and in the portfolio table's
  program panel. A program without a built folder cannot be pinned; the button says why. An approver also gets
  **Pin for <the portfolio's viewers>** (a group favorite) and, when the folder is missing, **Build its folder**
  (this program's folders only).
- **Board and table.** A pinned program's card says *pinned*; the portfolio table has a `pinned` column.
- **Home.** The **My programs** widget: a tile per pinned program the user can see, across portfolios. A tile shows
  code and name, portfolio, stage and its place on the track, status, the next gate's readiness (or *ready*), the
  attention count and *new data* when a followed program has news the user has not seen. A group's pin is marked
  *for <group>*. A click opens the program's page. With nothing pinned, the widget says how to pin.
- **How it is read.** One call for the user's favorites, one for their groups' (`getFavoritesForGroups`), the
  programs by `space_id`, then the portfolio table's facts for those programs only, so readiness and attention are
  the table's (`readinessFor`). A program row the user cannot see gives no tile.
- **The widget** is a `meta.role: dashboard` function (`myProgramsWidget`, `meta.order: 2`), which Home lists.
  The user removes it with its close icon and brings it back under *Customize widgets…*; Home has no reordering.
  It reads again each time Home comes back into view.
- **Address.** `/apps/ProjectPortalCore/Portal/<portfolio>/<program>` opens a program's page; the app keeps the
  address current as programs are opened.
- `pinProgram(portfolio, program, group)`, `unpinProgram(...)` and `pinnedPrograms()` do the same by code. Pure
  rules in `src/model/pinning.ts`.

### Activity feed

What happened to a program, newest first, read from the portal's own records (F-21, FR-81). No
event store: each record already says who did what and when.

- **Kinds.** Gate decisions (by an approver, the registry or the data), evidence, runs saved, runs selected for a criterion, documents linked and approved, snapshots, compound role changes, registry changes and work-unit changes the sync saw, and new data. Evidence attached within two minutes of its decision is folded into the decision's entry ("with 2 pieces of evidence"); evidence added later is its own entry.
- **Sources.** `stage_transition`, `evidence`, `model_run`, `criterion_artifact`, `document`, `program_snapshot`, `compound_role_change`, the portfolio's last 50 `sync_run` rows (each run's change for the program, and the new-data lines followers were told, kept in `sync_run.changes.newData`), and each follow's last news with `follow.news_on`. The same news line told twice within 15 minutes (two followers, or a sync and a check) is one entry. Stage moves and role changes are read from the gate log and the role history, not from the sync's changes, so nothing shows twice.
- **A series.** Narrowed to a series, the feed keeps what is tagged with it (`series KIN-A: new`, `EXS-KIN1 series KIN-A: 2 new values`) and what is about its compounds (from the program dataset's `series` and `compound` columns): the series' "what is new".
- **Program page.** *Activity (n)* with a kind filter (each kind with its count), a series filter when the program has series, eight entries a page, Newer/Older. An entry over a record has **Provenance**. Each series row has **Activity**, which narrows the section to it. In Simple mode the section is behind More.
- **My feed.** The ribbon's **My feed** opens everything the user follows: each followed program's feed, or for a program followed only by series, the series' part. A click on the program code opens it.
- **Functions.** `programFeed(portfolio, program, series, kinds, limit)` and `myFeed(kinds, limit)` return the feed as a table (at, portfolio, program, kind, what, who, work unit, record); `kinds` is comma-separated.

### Modes and quick actions per persona

The board has a **Mode** and a **Persona** select (F-08, FR-30–32). Both are the user's own and are
remembered in the platform's user settings (`grok.userSettings`, entry `ProjectPortalCore`), so they
follow the user to another browser.

- **Advanced** (the default) shows the whole program page, as before.
- **Simple** keeps the heavier sections behind **More**: missing assays and suggested models, needs attention, saved runs, history. The gate, the work units and the compounds stay. The ribbon drops Sync history, Export and All portfolios.
- **Quick actions** head the program page: one card per action the persona uses most, each opening something the page already does. A card that cannot be used now is dimmed and says why on hover.

| Persona | Cards |
|---|---|
| Chemist | Open data (the SAR table), needs attention, missing assays, suggested models, follow |
| Biologist | Open further data (assay results), open data, missing assays, follow, open workspace |
| Modeler | Suggested models, save a run, choose a run, missing assays, open data |
| Process scientist | Save a run, choose a run, open data, open workspace, follow |
| Lead | Gate decision, portfolio table, needs attention, follow, open workspace |

- **Which persona.** The one the user chose. Until they choose, the portfolio's mapping of groups to personas decides (the first group the user is in), then its default, then the setting: chemist in discovery, process scientist in development.
- **A portfolio may declare its own** (optional):

```json
"personas": {"default": "chemist",
  "groups": [{"group": "Portal EX-SMD Approvers", "persona": "lead"}],
  "cards": {"modeler": ["suggested_models", "save_run",
    {"name": "Portfolio as a table", "function": "ProjectPortalCore:portfolioTable"}]}}
```

  A card is an action code (`open_data`, `open_other_data`, `needs_attention`, `missing_assays`, `suggested_models`, `save_run`, `choose_run`, `follow`, `open_workspace`, `gate_decision`, `portfolio_table`) or a function. A function gets `portfolio`, `program` and `stage` when it has inputs by those names; a view it returns opens, a table opens as a table view. A persona the file leaves out keeps the core's cards.
- EX-SMD and EX-BIOC map their approvers to the lead persona; their defaults are chemist and process scientist.

### Program roles

A program may have its own viewers, contributors and approvers, on top of the portfolio's three
groups (F-16, FR-90). A role is held by a group or by a user (the user's personal group). Roles nest:
an approver may do what a contributor may. A program role only adds; it never takes away what the
portfolio's groups give.

- **Stored** in `program_role`: program, group, role, who granted it and when. One row per program and group; removing a role turns the row off (`active`).
- **Managed** on the program page (*Program roles*, behind More in Simple mode) and with `setProgramRole(portfolio, program, principal, role)`, by the portfolio's approvers only. A program approver cannot change roles. `programRoles(portfolio, program)` lists them.
- **Enforced in the app** where the portfolio's groups are: a program approver records gate decisions (and Met / Waive), a program contributor or approver selects runs. `recordDecision` and `selectRunForCriterion` check again before they write. Letting data advance a whole portfolio stays with the portfolio's approvers.
- **The workspace folder** of the program is shared to match: View to viewers, Edit to contributors and approvers. Removing a role removes the grant. The portfolio's own groups are never touched. A folder built later gets the roles from the builder.
- **The program row is not shared yet.** Program rows refuse grants (GROK-21041, see Tests), so a role does not make the row visible to someone outside the portfolio's groups. `syncRowGrants` in `src/program-roles.ts` applies the same grants to the row; it runs when `ROW_GRANTS` is set to true. Once GROK-21041 is fixed, set it, and the same table drives the row grants. The two isolation tests stay unchanged as the acceptance check.
- `program_role` is table-secured with View and Edit for all users, like `follow`: the app enforces who may change it, the server does not.

### Decisions with evidence

A gate decision rests on evidence (F-25, FR-71/74). The decision is the gate-log entry
(`stage_transition`: program, decision, from, to, who, when, rationale, source); the `evidence` table
says what it rested on. There is no separate `decision` table.

- **What a piece is.** `kind` (`model_run`, `snapshot`, `file`, `url`, `criterion`), `target`, a soft reference (a run's `run_id`, a file path, a URL, a gate criterion's id, a snapshot's id) as plain text with no foreign key, and `name` and `context` frozen when attached (a run's model, the criterion it was selected for, who chose it and why; a criterion's state, basis and what the data said). `role` is `readiness` for a criterion as the gate saw it, `supporting` for the rest.
- **It survives its target.** The history looks each target up again: a run renamed shows `(now <new title>)`, a run deleted shows `(gone; as recorded)` and a click shows what was recorded. URLs are not fetched.
- **When deciding.** The gate panel suggests the gate's readiness (one tick for its required criteria) and the run selected for each of its criteria, all ticked. Untick what the decision does not rest on; type a link or a file path and **Attach**. Advance, Hold and Terminate write the ticked pieces beside the gate-log entry.
- **By API and by data.** `recordDecision` without a list writes the suggestion, so an advance by data records its readiness the same way (the free-text `stage_transition.evidence` is kept). An empty list writes none. A piece that fails its check is reported and left out; the decision stands.
- **Afterwards.** The program's history lists each decision with its evidence, each piece clickable: a run shows its saved call in the context panel, a criterion shows its stage's gate, a URL opens, a CSV file opens as a table. `addEvidence(portfolio, program, kind, target, name, decision)` adds a piece to a recorded decision (default the latest), dated when added; the program's approvers only. `decisionEvidence(portfolio, program, decision)` lists decisions and evidence, with each target's state (`live`, `renamed`, `gone`, `unchecked`).
- `evidence` is master-secured through the program, like the gate log.

### Program snapshots

A snapshot freezes a program at one time (F-32, FR-15), so a decision can point at what the team
saw, not at what the data says today.

- **What it holds.** `program_snapshot`: when and by whom, a label and note, the stage and status, every criterion that applies with its state (the next gate's data conditions as evaluated, with what the data said), the run selected for each criterion, and the program dataset. Every column is write-once; nothing in the package updates or deletes a snapshot.
- **The table** is written as a CSV file to the program's workspace folder (`snapshot-<program>-<time>Z.csv`), and its SHA-256 is kept on the row. Opening the snapshot reads the file back and says whether it still holds what was frozen (intact, changed, missing). A program without a workspace folder keeps the table's size and fingerprint only, and the note says so. A failed write of the row removes the file it wrote.
- **Taking one.** **Snapshot now** on the program page (*Snapshots*, behind More in Simple mode), `snapshotProgram(portfolio, program, label, note)`, or **Snapshot now** in the gate panel when deciding: taken before the decision changes anything and attached as evidence. Contributors and approvers of the portfolio or the program.
- **Reading one.** Each snapshot on the program page, and each snapshot in a decision's evidence, opens read-only: stage, criteria, runs, the table's state, and **Open table**, a table view with editing off. `programSnapshots(portfolio, program)` lists them.

### Documents

A document is a record over one file (F-27, FR-70), so a gate's document criteria point at real
files with a version and a status, not at a tick.

- **What it holds.** `document`: the file's full path as the platform reads it (`System:AppData/...`, or a file of the program's workspace folder, `Space:Program:Files/name`), the folder when it is the program's, its kind, title, version, status (`draft`, `approved`, `superseded`), who wrote and linked it and when, who approved it and when, and the document criterion it answers.
- **Kinds** come from the portfolio's optional `documentKinds` (`[{code, name, description}]`). Without it, each document criterion's code is a kind. `other` is always there.
- **Approval meets the criterion.** Approving a document meets the criterion it answers unless the criterion was waived, the same way selecting a run does. A draft puts a criterion that was not started in progress. Approving a newer version supersedes the older approved document of the same kind and criterion. A superseded or withdrawn document never un-meets a criterion. Contributors and approvers link drafts; only approvers approve.
- **Document cards** on the program page (*Documents*): one per document criterion of the current gate (the program's stage) and the next, each with its state, its documents (approved first, newest version first; each opens: text in a dialog, CSV as a table, anything else downloads), what is missing, the run selected for the criterion, and the decisions that cite one of its documents. **Link document…** picks a file of the program folder, takes a path, or uploads a file into the folder; a file of the same name is never overwritten. Documents no card shows (no criterion, or another stage's) follow the cards.
- **As evidence.** A document is evidence of kind `file`, its path the target and its record in the context. The gate panel suggests the approved documents of the gate's criteria, ticked; `addEvidence` takes a document by id or path; the history finds a document by its record (renamed, gone).
- **Functions.** `programDocuments(portfolio, program)` lists a program's documents; `linkDocument(portfolio, program, path, criterion, kind, title, version, status)` links a file or changes the document over it. `loadExampleDocuments()` links the synthetic EX-BIOC documents in `files/documents/` (EX-* portfolios only).
- `document` is master-secured through the program, like the gate log.

### Evidence tables

The evidence behind a document criterion, generated into a document (F-40, FR-72), so a written
document can say which runs, data and decisions it rests on.

- **What goes in.** For one document criterion of a program: the run selected for it, with the saved call's inputs and outputs and their units; the program's newest snapshot, with the criterion's state as frozen and up to 25 rows of its table; the data conditions of the criterion's gate, as readiness evaluates them now; the decisions whose evidence cites the criterion, its run, or a document over it that a person wrote. A part with nothing in it says so in one row.
- **What comes out.** One table (`section, source, ref, item, value, units, state, on`), written twice into the program's workspace folder: a CSV and a Markdown file with a table per section (`evidence-<criterion>-<program>-v<n>-<time>Z.md`). The Markdown file is linked as a new draft document over the criterion, kind `evidence_table` when the portfolio declares it, else `other`, version 1, 2, 3 per criterion. It goes on like any document: approved directly, or submitted where the portfolio asks for approval.
- **What it records.** `document.generated`: the criterion, when and by whom, each source (`run`, `snapshot`, `condition:<code>`, `decision:<id>`) with its soft reference and a content hash, one hash over all of them, and the CSV's path and SHA-256.
- **Again.** Generating again makes the next version. Older generated drafts of the criterion are superseded at once; an approved one stays approved until the new version is approved, which supersedes it the usual way.
- **In the app.** **Generate evidence table** on each document card (contributors and approvers); **Table** on a generated document opens its CSV as a table view. `generateEvidenceTable(portfolio, program, criterion)` does the same and returns the table.
- **Evidence changed since written** (F-41, FR-73). The same sources are read again and compared with the record: another run selected, a newer snapshot, a data verdict that changed, a decision added or no longer citing it. Any difference makes the document stale. The card shows **stale** and what changed, with **Regenerate** (the next version; older drafts superseded, an approved one kept until the new one is approved). The program page stores what it finds on `document.evidence_check`; the activity feed shows a stale document as an entry.
- **Who is told.** After every sync (as for new data) and on `checkEvidenceChanges(code)`, each stale document's author and the program's followers are told in the platform's notifications, once per state of the evidence. Superseded documents are not checked.
- Pure assembly in `src/model/evidence-tables.ts`; the platform side in `src/evidence-tables.ts`. Needs the program's workspace folder.

### Approval workflow

A portfolio may ask for a second person (F-30, FR-57). It is the ArtifactAlignment pattern
(publish, pending, approved) rebuilt on the portal's own tables. Off by default: without `approval`
every decision and approval takes effect at once, as before.

- **Opting in.** `approval: {decisions: ["advance", "terminate"], artifacts: ["document", "model_card"], approvers: 1}`. `decisions` are gate decisions (advance, recycle, hold, terminate; portal authority only). `artifacts` are documents and model cards. `approvers` is how many approvers other than the proposer must approve (1 to 5, default 1).
- **A decision becomes a proposal.** An approver's advance or terminate is stored in `approval_request` as pending, with its rationale, the stage it goes to, and the evidence frozen now (the gate panel's ticks, a snapshot taken now). The program does not move. Only one decision per program waits at a time. `recordDecision` refuses such a decision unless it carries out the pending request.
- **Review.** Another approver (a portfolio approver or a program approver) approves or rejects it with a rationale. The proposer cannot; they may withdraw it. On the last approval it goes through `recordDecision` as the reviewer: gate log (source as proposed), evidence, workspace folders. A rejection is logged on the request; nothing changes. The proposer is told in the platform's notifications. If the decision cannot take effect any more (the gate is no longer met), the review is not kept and the request stays pending.
- **Advance by data.** In a portfolio with both `autoAdvance` and approval of advance, data proposes the move instead of making it (source `data`), once per program. Any approver may approve a proposal by data, the one who ran the sync included. The sync report says "waiting for approval".
- **Documents.** A draft or rejected document is **Submitted for approval** (contributors and approvers): draft → submitted → approved or rejected. Approving goes through the same path as before: the criterion is met, the older version superseded. A submitted document cannot be changed or approved directly.
- **Model cards.** `submitModelCard(portfolio, model, version)`: draft → in review → approved or rejected, for a model in the portfolio's catalogue. Reloading the catalogue keeps an approval the workflow set.
- **In the app.** The gate panel's buttons say **Propose advance to …** and **Propose termination**, and are off while a decision waits. *Waiting for approval* on the program page (both modes) lists the program's requests: Approve and Reject with a rationale for an approver who did not propose it, Withdraw for the proposer, the last three closed. The history says who proposed and who approved a decision. The ribbon's **Waiting for me** lists every pending request the user may review, across portfolios.
- **Example.** EX-BIOD asks for one approver on advance, terminate, documents and model cards. No existing test or scenario decides in EX-BIOD.
- `approval_request` is table-secured, View and Edit for all users, like `program_role`: the app checks who proposes and reviews. Master security through the program would hide a request from every approver but the program's author until GROK-21041 is fixed.

### Provenance

From any portal record, what it rests on and what rests on it (F-28, FR-75): data → model run →
document → decision, and back. The graph is the portal's own records; nothing new is stored.

- **Nodes.** The program, its decisions (gate-log entries), model runs, the tables their saved calls read and wrote, the catalogued model a run used, documents, snapshots, gate criteria, the program dataset, a compound, and the files and links a decision cites.
- **Edges** run from what something rests on to what rests on it: a compound is in the program dataset; the dataset is frozen in a snapshot; a run reads its input tables and its model and writes its output tables; the run selected for a criterion feeds the documents answering it (or the criterion when none does); a document answers its criterion; each piece of a decision's evidence points at the decision; a decision moves the program. Evidence is an edge, not a node: starting from a piece of evidence starts from its target. A cited record deleted since is still a node, named as the evidence froze it.
- **The walk** goes `depth` steps (default 4) each way, breadth first. A compound's programs (up to 10) are read when the walk starts from a compound; otherwise one program's records. Saved calls are read for their tables for the 40 newest runs.
- **In the app.** **Provenance** on the program page, beside each decision and each piece of evidence in the history, and beside each saved run, document and snapshot. The dialog shows *Rests on* and *Used by* as indented trees; a click opens the record (a run's call in the context panel, a snapshot, a document's file, the program data, a model card, a table, another program on its page), and **Trace** starts the walk again from it.
- **Functions.** `provenance(kind, id, depth)` returns `{start, nodes, edges}` (kind: program, decision, evidence, model_run, snapshot, document, compound; a run's `run_id` and a compound's code work as ids). `showProvenance(kind, id)` opens the dialog.

### Stage advance driven by data

Off by default. A portfolio whose stages the portal owns may let data move its programs, up to a
limit stage, so the top gates still need an approver:

```json
"autoAdvance": {"through": "H2L", "countWaivers": false}
```

- **When a program moves.** It is active, already in its model, and its next stage is not past `through`. That stage has required criteria, and every one was evaluated met from the program's data. It moves one stage per pass.
- **What never counts.** A waiver, unless `countWaivers` is true: a waiver is a person's call. A criterion recorded met by hand, a document or a decision: not data. A stage with no required criteria: nothing for data to decide.
- **What never moves.** Held, terminated, completed, transitioned or missing programs. Unstaged programs: entering the model is a decision. Programs whose data cannot be read.
- **How it is recorded.** Like an approver's advance (`recordDecision`): a `stage_transition` with `decision` advance and `source` data, `evidence` naming each criterion and what the data said, the data conditions' results kept on `criterion_status`, and the new stage's folders added to the workspace. Approvers' decisions are logged with source `approver`, registry changes with `registry`.
- **When it runs.** After every sync of the portfolio, by hand or by schedule; the run's report and history name each program it moved (`EXS-GPC3: advanced by data HF → H2L`), and owners hear of it. Or on call: `autoAdvancePortfolio(code, dryRun)`. A dry run (the default) lists what would move and why each other program stays. Applying needs the caller to be an approver, as a manual advance does.
- **In the app.** Sync history says how far data advances stages and has **Preview advance by data**. The program's history marks a move `(by data)`; hover for the evidence.
- **It is a field, not a third stage authority.** Stage authority says who owns the stage (registry or portal); auto-advance is a way the portal decides, so the gate rules, the approver groups and the feed's handling stay as they are.
- The example EX-SMD lets data advance through Hit-to-Lead. Its programs are unstaged, so nothing moves until an approver places one in the model; a program at Hit Finding whose data shows a validated series and a sub-micromolar compound then enters Hit-to-Lead on the next sync.

### Choosing the run a criterion rests on

A criterion can rest on a model run (F-15). It says so in the portfolio file:

```json
{"code": "process_description", "name": "Locked process description", "kind": "document",
 "artifact": "model_run", "runFilter": {"unitOperations": ["EX-PROD"]}}
```

- **The filter** narrows the program's runs by `unitOperations`, `runTypes` and `models` (catalogue codes; a model by its nqName). Every list it has must match; a run without the code does not pass. No filter offers every run.
- **The gate panel** shows the run selected for the criterion (who, when, why) and the runs it superseded. **Choose run…** opens the dialog.
- **The dialog** lists the runs the filter admits and says how many it leaves out. Tick two or more to compare them: inputs, then outputs, labelled from the function with units, what differs in bold, a table output by its row count. Pick one and give a rationale.
- **Who may choose:** the portfolio's contributors and approvers, also where the registry owns the stage. Choosing a run is evidence, not a stage decision. Not on a program that has ended. `selectRunForCriterion` asks again before it writes.
- **What a choice does** is unchanged: the run becomes `selected`, the one before `superseded` (kept, shown as history), and the criterion is met unless it was waived.
- `criterionRuns(programId, criterionId)` and `runComparison(items)` return what the panel and the dialog show.

### Saving a run from the app

The **Save to program** dialog (F-14, FR-53) asks for portfolio, program, unit operation, run type,
purpose and an optional title, then calls `saveModelRun` and says where the run was linked
(`Saved to EXB-203 at D5, linked into …/Drug Substance/Process Lock`).
- **Filled in before you touch it.** From the program page: its portfolio and program. Otherwise from the call's own annotations when it was saved before, else from the model's catalogue: its first portfolio you can see, its first unit operation, its default run type. A code the lists do not hold is left blank.
- **Programs you can see.** The program list is the portfolio's programs under row security, running ones first.
- **Validity.** An input outside the card's validity range is a warning, not a refusal.
- **Where it opens.**
  - After every run of a catalogued model, a balloon offers *Save to program…* (for 20 s). It follows runs from any run view, the function dialog and Compute's included, once the package has loaded in the session (opening the app loads it).
  - A run's context menu (`Save to program…`), wherever the platform shows one for a FuncCall.
  - The program page's **Save a run…** button, over this session's runs of catalogued models.
  - `saveRunToProgram(call)` from a script.
- Compute2's RichFunctionView has no extension point for a third-party action, so the dialog is not a button in its ribbon.

`model_run.func_call`, the soft reference to the saved call, cannot be written today (platform
gap, see Tests). Every row also holds `run_id`, the call's id as plain text. The insert tries the
soft reference first and writes the row without it when the server refuses, so the reference is
written once the platform is fixed. **Migration debt:** rows without `func_call` need it backfilled
from `run_id` then (`func_call` is write-once through the API, so this is a SQL migration), after
which `func_call` can be made required again.

## Configuration

Portfolios are data, not code. Each portfolio is a JSON file: its groups, stage model (stages and
gate criteria), stage authority, and registry feeds. The package ships three examples in
`files/portfolios/`, on synthetic registries written as `VALUES` on `System:Datagrok`, so every stand
can run them:

| File | Registry | Stage model | Authority |
|---|---|---|---|
| `ex-sm-discovery.json` | example small-molecule registry, 5 projects | EX-SMD-V1: TV → HF → H2L → LO → CS → PCD | portal |
| `ex-bio-discovery.json` | example biologics registry, research slice | EX-BIOD-V1: TI → LG → LO → DA → CN (hands off to EX-BIOC) | portal |
| `ex-bio-development.json` | example biologics registry, development and commercial slice | EX-BIOC-V1: D1 → D7 | feed |

The example data covers the cases a real registry raises: a key no portfolio claims, a molecule
handed on from discovery (lineage), one still active in discovery while in development (overlap), a
development gate code on a research project, and a stage spelling that resolves to nothing.

Customer portfolios live in customer packages, not here. A customer package hands its files to the
core with `loadPortfolioConfigs` and audits its sources with `routingAuditOf`.

### Profile tooling

Portfolio files load and validate as a set. Four functions work on them as files (F-18):

- **`exportPortfolio(code)`** reads a loaded portfolio back into the file format: stage model, criteria, feeds, groups by name, datasets, rules, schedule, auto-advance, handoff, and its catalogue (the unit operations and run types it declared, the models that name it, each with its current card). Loading the export gives the same portfolio. Rows a newer file no longer declares are still exported, because loading never deletes.
- **`diffPortfolios(a, b)`** compares two portfolios, each a loaded code or a file's JSON. Both are normalized first (defaults written out, empty lists dropped, keys sorted, cards completed from the function as loading does). The lines come in three groups: vocabulary (stages, criteria, work units, series maturity, catalogue codes), workspace tree (root, workstreams, stage folders), other fields (feeds, groups, datasets, rules, schedule). Long text such as SQL is reported by length, not printed.
- **`exportProfile(codes)`** bundles loaded portfolios (default: all) as `{format: "project-portal-profile", version, exportedOn, portfolios}`.
- **`importProfile(bundle, dryRun)`** checks the bundle's portfolios as a set (a handoff may name a portfolio already loaded), diffs each against the loaded one (new, changed, unchanged), and unless `dryRun` (the default) loads them. A bundle that does not validate is refused whole. Loading adds and updates; removals in the diff are reported, not applied.
- **In the app.** The board's **Export** shows the open portfolio's file, with a download.

`files/sources.json` lists each source system's full key query for the routing audit.
`validatePortfolioSet` checks the files as a set before anything is written. It checks:
- duplicate codes;
- a stage spelling that could resolve to two stages;
- a shared stage model declared two different ways;
- a handoff that doesn't leave from a terminal stage, or goes to a portfolio that isn't configured.

## Functions

| Function | What it does |
|---|---|
| `Project Portal` (app) | Portfolios, the selected portfolio's own board, and the gate for the selected program, with a link to the program's workspace. `/apps/ProjectPortalCore/Portal/<code>` opens one portfolio, `.../<code>/<program>` one program's page; `.../<code>/table` its table, `.../Portal/table` every portfolio's. |
| `loadPortfolios` | Loads this package's portfolio files: groups, stage models, portfolios, feeds. Idempotent. |
| `loadPortfolioConfigs(configs)` | Loads portfolios handed over by another package, as a JSON array, checked as a set |
| `portfolioTable(portfolio)` | The programs of a portfolio, or without one of every portfolio you can see, as a table: stage, days in stage, readiness for the next gate, compounds by role, work units, workspace |
| `programData(portfolio, program, dataset)` | One program's rows of the program dataset, or of a further dataset by code, as a table named after the program. Before 0.3.0 it was registered as `programDataTable`; that name still works |
| `openProgramData(portfolio, program, dataset)` | Opens that table as a table view with the portfolio's viewers |
| `recommendModels(portfolio, program, stage, series, task)` | Catalogued models that fit the program, best first, a model predicting an assay the gate lacks first |
| `missingAssays(portfolio, program, stage, series)` | The columns the gate's open data conditions need and rows lack, with the suggested models that predict each |
| `setProgramRole(portfolio, program, principal, role)` | Gives a group (by name) or a user (by login) a role on one program: viewer, contributor, approver, or none to remove it. Portfolio approvers only; shares the program's folder to match |
| `programRoles(portfolio, program)` | The program's roles as a table |
| `decisionEvidence(portfolio, program, decision)` | The program's gate decisions with their evidence, one row per piece, and whether each target is live, renamed or gone |
| `snapshotProgram(portfolio, program, label, note)` | Freezes the program now: stage, criteria, selected runs, and the program dataset as a file in its workspace folder |
| `programSnapshots(portfolio, program)` | The program's snapshots, newest first |
| `programDocuments(portfolio, program)` | The program's documents: criterion, kind, title, version, status, path |
| `proposeDecision(portfolio, program, decision, rationale, stage)` | Proposes a gate decision in a portfolio that asks for approval of it; approvers only |
| `submitDocument(portfolio, program, path, rationale)` | Submits a document for approval (draft or rejected → submitted); contributors and approvers |
| `submitModelCard(portfolio, model, version, rationale)` | Submits a model card for approval (draft or rejected → in review) |
| `reviewRequest(request, verdict, rationale)` | Approves or rejects a request; on the last approval it takes effect. Never the proposer |
| `withdrawRequest(request)` | Withdraws a pending request; its proposer only |
| `approvalRequests(portfolio, program, status)` | A portfolio's or program's requests: what, status, progress, proposer, reviews, result |
| `waitingForMe()` | The pending requests the current user may review, oldest first |
| `linkDocument(portfolio, program, path, criterion, kind, title, version, status)` | Links a file as a document of the program, or changes the document over it; approving meets its criterion unless waived |
| `generateEvidenceTable(portfolio, program, criterion)` | Generates the evidence table of a document criterion into a new draft document version in the program folder (CSV and Markdown) and returns the table |
| `checkEvidenceChanges(code)` | Marks each generated evidence table of the portfolio stale or current against the evidence now; tells a stale one's author and the program's followers once. The sync runs it too |
| `provenance(kind, id, depth)` | What a portal record rests on and what rests on it, `depth` steps each way, as `{start, nodes, edges}` |
| `showProvenance(kind, id)` | The same in a dialog, each record clickable |
| `loadExampleDocuments()` | Links the synthetic EX-BIOC documents the package ships. EX-* only; idempotent |
| `addEvidence(portfolio, program, kind, target, name, decision)` | Adds a run (by `run_id` or title), a criterion (by code), a URL or a file path to a decision (default the latest), its name frozen now. Approvers only |
| `checkNewData(code)` | Checks every followed program of a portfolio for new data and notifies the followers; also runs after every sync |
| `followProgram(portfolio, program, series)` | Follows a program, or one of its series, as you; the data now is the baseline |
| `unfollowProgram(portfolio, program, series)` | Stops following |
| `scoreCompounds(portfolio, program, model)` | Runs a catalogued predictive model once over the program's compounds, saves the run, stores one prediction per compound and predicted column; idempotent per model version |
| `programPredictions(portfolio, program, model)` | The program's latest predictions with the measured value beside each |
| `programFeed(portfolio, program, series, kinds, limit)` | The program's activity, newest first, as a table; a series narrows it to that series, `kinds` (comma-separated) to those kinds |
| `myFeed(kinds, limit)` | The activity of everything you follow, newest first |
| `programAttention(portfolio, program)` | The items of one program that trip the portfolio's attention rules, one row per item and rule |
| `programReadiness(portfolio, program, stage)` | Every criterion of a stage (default: the next one) met, waived or open, with data conditions evaluated, and the gap list |
| `syncPortfolio(code, dryRun)` | Runs one portfolio's program feeds, then its work-unit and compound feeds. With `dryRun` it prints the plan without applying it. |
| `syncDuePortfolios(codes)` | Syncs only the portfolios whose schedule is due, records each run, notifies owners of a run that changed something |
| `syncHistory(portfolio)` | A portfolio's sync runs as a table, newest first |
| `exportPortfolio(code)` | A loaded portfolio as its portfolio file (JSON) |
| `diffPortfolios(a, b)` | Differences from `a` to `b` (each a loaded code or a file's JSON): vocabulary, workspace tree, other fields |
| `exportProfile(codes)` | Loaded portfolios as one bundle, to move them to another instance |
| `importProfile(bundle, dryRun)` | Validates a bundle as a set, diffs it against what is loaded, and unless `dryRun` loads it; refused when it does not validate |
| `autoAdvancePortfolio(code, dryRun)` | Advances each program whose next gate the data alone meets, one stage, never past `autoAdvance.through`. `dryRun` (default) says what would move and why the rest stay |
| `syncAllPortfolios` | Syncs discovery portfolios before development ones, so lineage resolves in one pass |
| `portfolioRoutingAudit` | For each source: keys no portfolio claims, and keys claimed twice |
| `routingAuditOf(sources)` | The same audit for sources handed over by another package |
| `buildWorkspace(code, dryRun, program)` | Builds one portfolio's workspace, or with `program` one program's folders. With `dryRun` it lists what it would create. Refused while another build of the same root runs. |
| `workspaceOrphans(code)` | The folders a prune-or-park pass would look at; a report, nothing is changed |
| `parkOrphans(code, dryRun, folders)` | Moves those folders (or the ones named) into `_Parked/<date>` under the root and records each move; never deletes. Dry run by default; approvers apply |
| `unparkFolder(id)` | Moves a parked folder back where it was; approvers only |
| `parkedFolders(code, all)` | The portfolio's parked folders, and with `all` those restored since |
| `pinProgram(portfolio, program, group)` | Pins the program's folder in your Workspace (a favorite), or a group's (approvers only), so it gets a tile in My programs |
| `unpinProgram(portfolio, program, group)` | Takes the pin away |
| `pinnedPrograms()` | Your pinned programs across portfolios: stage and step, status, next gate's readiness, attention, new data |
| `myProgramsWidget` | The My programs widget on Home |
| `setCompoundRole(portfolio, program, compound, role, rationale)` | Changes a compound's role in a program by hand, logged |
| `saveModelRun(call, portfolio, program, unitOperation, runType, purpose, title, folder)` | Annotates, saves and links a run into the program's workspace, and records it as a `model_run` |
| `saveRunToProgram(call)` | Opens the save-to-program dialog for a run |
| `modelCard(model)` | The card of a catalogued model (by function nqName), as a widget |
| `loadModelCatalogue(catalogue, portfolio)` | Loads a catalogue handed over without a portfolio file, checked first. Idempotent |
| `buildAllWorkspaces` | Builds the workspace of every portfolio that declares a `workspaceRoot` |

## Tests

- **`grok test`: package tests.**
  - Configuration rules, the sync planner, gate rules, data conditions and readiness. These are pure functions.
  - Compounds on the example data: roles, a compound in two programs, history, a manual change, idempotence.
  - Work units on the example data: series, campaigns, studies, process stages; a rejected row; idempotence.
  - Readiness on the example data: hit-to-lead and lead-optimization conditions, development lots, the gate.
  - The portfolio table: one row per program, stage order, readiness and days in stage against `programReadiness` and the transition log, compounds and work units, filters, the linked readiness matrix, and a restricted user's rows.
  - Integration against the example portfolios: feeds, authority, lineage, overlap, routing, gate decisions.
  - Schema round-trips: the run-selection loop and the workspace vocabulary.
  - The workspace builder: tree, idempotence, preview, a renamed program, concurrent builds, stage folders from the builder and on advance, the orphan report.
  - Pinned programs: the rules (what can be pinned, who pins for a group, pins from the user's and the groups' favorites, a tile's track, readiness and address, tiles that follow the pins, unseen news), and on EX-SMD a pin and unpin of EXS-PRO2 through favorites (its tile, `pinnedPrograms`, the table's `pinned` column), an unpin in the platform taking the tile away, a pin for the viewers group, a program without a folder refused, and a program's address opening its page. The admin's favorites are put back.
  - Parking: the rules (where each orphan goes, already parked, only some, who may, unparking over a folder of the same name, the parking folder never an orphan, the summary), and on EX-SMD a throwaway folder previewed, parked, left out of the report and the build, unparked with `_Parked` removed and the tree as before; a terminated scratch program's folder parked and not rebuilt; a non-approver refused.
  - Isolation between portfolios, of program rows and of the workspace, run as a throwaway restricted user.
  - Model runs: the example model, a run saved at process lock (annotations, the stage folder, the reload, the `model_run` row by run id), and a gate document criterion selecting it.
  - Saving from the app: the dialog's rules (what it fills in, the program list, what it refuses, the result line) and the dialog itself, driven in the test: filled in for a Fed-batch titer run, a validity warning, saved to a program.
  - Scheduled sync: the rules (intervals and cron in UTC, when a sync is due, a run's counts and per-program changes, when owners hear of it), and on EX-SMD and EX-BIOD a recorded run, a due and a not-yet-due schedule, a platform notification, the function and the history panel.
  - Auto-advance: the rules (one stage, never past `through`, waivers, recorded criteria, held and unstaged programs, the setting), and on EX-SMD a dry run, a real pass with its gate-log entry and evidence, a sync that advances and its history, and the example put back.
  - Profiles: the rules (normal form, readable diff lines by kind, bundle checks, refusal, a handoff to a loaded portfolio), and on the examples the export of each EX file equal to the file, the functions on an edited EX-SMD, and a bundle dry-run, refused, imported and put back.
  - Attention: the rules (what trips one, severity, messages, the lead set, what a rule needs), and on EX-SMD the flags of EXS-KIN1 and EXS-GPC3, an EX-BIOC lot, the function, the portfolio table's count, and a flagged compound opened.
  - Modes and personas: the rules (the persona from groups, default and setting; cards from the defaults, the portfolio and functions; the config check; what Simple hides), and on EXS-KIN1 the settings remembered, the program page's cards for a chemist and a modeler, Simple with More and back, and a function card that opens a table.
  - Decisions with evidence: the rules (what an item needs, a typed link, readiness as evidence, the suggestion, renamed and gone targets), and on EX-SMD a scratch program advanced with its readiness, with an explicit list and a refused item, a run added later then renamed and deleted, the gate panel attaching a link, and an advance by data with its readiness evidence.
  - Program snapshots: the rules (criteria frozen in stage order, recorded or evaluated; the file name; the table check; the summary), and on EX-SMD a snapshot of EXS-KIN1 with its SAR table in the program folder, read back intact and then noticed changed, the functions and the dialog, and a decision that takes a snapshot of a program without a workspace.
  - Documents: the rules (the vocabulary, what a document needs, who may link and approve, what a status does to the criterion, superseding, order, the cards, a document as evidence, the config check), and on EX-BIOC the example documents and EXB-201's card, a scratch project's draft, approval, a newer version, a waived criterion, the refusals and the functions, and a file uploaded into EXB-201's folder.
  - Provenance: the rules (the graph's edges, a cited record deleted since, a walk from a decision up to the data and down to the program, depth, the tree, starts from evidence, a run id and a compound code, the text), and on EX-BIOC a document resting on the run selected for its criterion and its model, the run's output table, and on EX-SMD a decision's readiness and snapshot up to the program data, a start from evidence, the function, the dialog and the program page's links.
  - Approval workflow: the rules (the config check, what needs approval, who may propose and review, N approvers and one rejection, withdrawal, Waiting for me, a proposal by data), and on EX-BIOD a scratch program's advance proposed, refused to its proposer and to a non-approver, approved by another person and recorded; a rejection; a hold recorded at once; a document submitted and approved, meeting its criterion; a model card approved and kept on reload; and on EX-SMD an advance by data proposed once and recorded on approval.
  - Program roles: the rules (the strongest role, roles only add, who may change them, grants and what changes), and on EXS-KIN1 a throwaway user made viewer, contributor and approver (access in the app, the folder's grant, the folder seen as that user), a downgrade, removal; a non-approver refused; the functions.
  - New-data notification: the rules (snapshot, baseline, growth only, a value in a blank cell, a watermark that sees a replaced row, the config check), and on EXS-KIN1 a follow and a series follow with their baselines, a quiet check, news told once to the follower, a series follower told only about the series, what's new since the last visit, the sync hook, and the functions.
  - Predictions: the rules (output columns, units, what a scoring stores and keeps, the latest version, predicted against measured, missing assays), and on EXS-KIN1 a solubility scoring stored against its saved run, the same version again running nothing, predicted vs measured and the missing assay's stored prediction, the functions, and the provenance from the run to the model and the predictions.
  - Activity feed: the rules (entries from each kind of record, evidence folded into its decision, news told twice counted once, newest first, a series' part, kinds, pages), and on EXS-KIN1 new data a check found in the program's feed and the follower's, a series follow that sees only its series, and the functions.
  - Model recommendation: the rules (applicability keys, missing assays per criterion and series, matching and order, a dialog's prefill), and on EX-SMD the gaps of EXS-KIN1 towards Lead Optimization and the model that predicts one, the functions, a stage no model is for, a task filter, an EX-BIOC program offered only its own catalogue, and the example models run on the SAR table.
  - Program data: the rules (the dataset list, what a further dataset needs, the table's name, default viewers), and on EX-SMD the SAR table and assay results of EXS-KIN1, what is refused, and the table view.
  - Choosing a run: the rules (what a criterion takes, the filter, who may choose, the comparison, the history), and on EX-BIOC the candidates of D5's process description, two runs compared, a choice superseded, and the dialog.
  - The model catalogue: its rules (codes, cards, merging across portfolios, validity ranges), the example catalogue loaded twice with nothing rewritten, the card joined with the function, the card rendered, and a run saved with the catalogue's codes.
- **`npm run test:bdd`: Gherkin scenarios** in `bdd/features/`, compiled to Playwright by
  `@datagrok-libraries/bdd`. Features address controls only by their `data-u2-name`.

On a bleeding-edge stand dated 2026-10-02, 280 of the 283 package tests and all 33 BDD scenarios pass.
The two isolation tests fail because of a platform gap, not a portal bug:
- The server refuses `permissions.grant` on a promoted domain row, even from the row's author with Share (`~can_share` true).
- `privileges_service.savePermission` checks Share against the projects repository, which never holds domain rows.
- Programs therefore stay private to their author, which fails closed. The tests stay in place as the acceptance check for the fix.
- Program roles wait on the same fix for row grants; the hook is `ROW_GRANTS` in `src/program-roles.ts` (see Program roles).

The third, saving a model run, fails on a second platform gap: the server's existence check for a
soft reference to `Core.func_calls` aliases its target `p`, and the permission subquery it embeds
declares its own `permissions p`, so `p.func_id` resolves against the wrong table. The test stays as
the acceptance check; `saveModelRun` falls back to `run_id` meanwhile (see Model runs).

## Known platform issues found building this

1. **Sharing a domain row is refused.** A promoted domain row cannot be shared through `grok.dapi.permissions.grant`. See above.
2. **`string_list` columns can't be written from JS.**
   - `Types.fromValue` reports any `List` as `list`, so a JS array never validates against a `string_list` column.
   - The workaround is to hold stage spellings in `stage.spellings` as a `|`-separated string.
3. **The npm-installed `grok` lags the repo.** It rejects current `domain-schema.schema.json` manifests, Grit's included. Use the workspace build at `node_modules/.bin/grok`.
4. **`dapi.connections.filter('name = "X"')` matches the friendly name.** `findConnection` in `src/portal.ts` falls back to an exact match.
5. **A range filter over an all-zero numeric column drops every row.** A column of `[0, 0, 0]` in the filter panel leaves 0 rows passing; `[5, 5, 5]` keeps all three. The portfolio table leaves a numeric column out of the filter panel until its values vary.
6. **The bar chart sorts categories by name, not by the column's category order.** Sorting by category ignores `setCategoryOrder`, while the box plot respects it. The table's bar chart splits on a hidden `stage in order` column (`01 D1`, `02 D2`).
7. **Listing a space's linked runs failed once (GROK-21043, closed as cannot reproduce).** In a test run on 2026-10-01,
   `GET /spaces/{id}/children?includeLinked=true` failed with `Converting object to an encodable object failed` once
   the space linked a run. Rechecks the same day did not reproduce it: existing program and stage folders holding linked
   runs, and a fresh space with a run just linked, list them (with and without a `FuncCall` types filter), and the
   folder view shows them. The ticket was closed; reopen it with a repro if it recurs. The tests still read
   `grok.dapi.projects.find(id).links`.
8. **A panel that takes a function is not offered for one.** A `meta.role: panel` function whose input is a `func`
   (or an `object`), with a condition that holds, does not appear in the context panel of a model script. The
   model card is therefore reached from the program page and through `modelCard(model)`, not as an info panel.
9. **A package cannot schedule its own functions.** The server's scheduler (`#schedule:` with an optional
   `#schedule.runAs:` group) runs only server-side scripts and queries, and the JS API has no way to create a
   schedule. A package's JavaScript function does not run headless either. So `syncDuePortfolios` exists but is not
   registered; it runs when someone, or something with a browser session, calls it.
10. **No JS API to send a notification.** `POST /users/notifications` creates one, but only for the caller unless they
   hold Edit users, and the stand's `grok.dapi` has no notifications client at all. The portal calls the route with
   `fetch`.
11. **A column named `user` breaks the schema deploy.** The deployer does not quote column names in its DDL, so
   `"user": {"type": "user"}` fails with `syntax error at or near "user"`. `follow` names it `user_id`.
12. **A row of a table without soft delete cannot be deleted.** On a table declared `"softDelete": false`, both
   `delete(id)` and a transaction's delete answer `not-found` for a row that `query` returns. Changing `softDelete`
   afterwards is refused without a migration. `follow` keeps rows and turns them off instead.
13. **`FuncCall.edit()` runs a call whose inputs are all set** instead of showing its dialog. The portal's Run builds
   the dialog with `buildEditor`.
14. **A package function registers under its TypeScript name.** A different `//name:` header is ignored. The
   function documented as `programData` was therefore registered as `programDataTable`; since 0.3.0 the TypeScript
   function is named `programData`, and `programDataTable` stays as a second function for scripts that use it.
15. **Deleting a user setting that is not there throws.** `grok.userSettings.delete(name, key)` fails with
   `Cannot read properties of null` when the key is absent. The portal deletes only a key it has read.
16. **`in (...)` over a string column of uuid-shaped values fails.** `space_id in ("…", "…")` on `program.space_id`
   (a string) answers `operator does not exist: text = uuid`; `space_id = "…"` works. The pinned programs are read
   with `=` joined by `or`, 50 at a time.
