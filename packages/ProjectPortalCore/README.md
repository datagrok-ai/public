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
  - **Orphans are reported, never touched (FR-22).** Every full build, and `workspaceOrphans`, lists root folders no program claims and the folders of programs that are terminated, transitioned, completed or missing from their feed. Pruning or parking them is a separate, explicit step this package does not take.

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
- **Per-program changes** (new, fields changed, stage moved, no longer in the feed, work units, compounds added, role changes) are on the result and in `sync_run.changes`. Data-driven stage advance and new-data notification will read them.
- **Nothing runs by itself yet.** The platform schedules only server-side scripts and queries (`#schedule:`), not a package's JavaScript functions, so the package cannot register its own schedule. Something must call `syncDuePortfolios` (see Known platform issues).
- The examples declare schedules: EX-SMD every 6 h, EX-BIOD daily at 06:00 UTC, EX-BIOC weekdays at 06:00 UTC. Owner: `admin`.

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

`files/sources.json` lists each source system's full key query for the routing audit.
`validatePortfolioSet` checks the files as a set before anything is written. It checks:
- duplicate codes;
- a stage spelling that could resolve to two stages;
- a shared stage model declared two different ways;
- a handoff that doesn't leave from a terminal stage, or goes to a portfolio that isn't configured.

## Functions

| Function | What it does |
|---|---|
| `Project Portal` (app) | Portfolios, the selected portfolio's own board, and the gate for the selected program, with a link to the program's workspace. `/apps/ProjectPortalCore/Portal/<code>` opens one portfolio; `.../<code>/table` its table, `.../Portal/table` every portfolio's. |
| `loadPortfolios` | Loads this package's portfolio files: groups, stage models, portfolios, feeds. Idempotent. |
| `loadPortfolioConfigs(configs)` | Loads portfolios handed over by another package, as a JSON array, checked as a set |
| `portfolioTable(portfolio)` | The programs of a portfolio, or without one of every portfolio you can see, as a table: stage, days in stage, readiness for the next gate, compounds by role, work units, workspace |
| `programData(portfolio, program, dataset)` | One program's rows of the program dataset, or of a further dataset by code, as a table named after the program |
| `openProgramData(portfolio, program, dataset)` | Opens that table as a table view with the portfolio's viewers |
| `programAttention(portfolio, program)` | The items of one program that trip the portfolio's attention rules, one row per item and rule |
| `programReadiness(portfolio, program, stage)` | Every criterion of a stage (default: the next one) met, waived or open, with data conditions evaluated, and the gap list |
| `syncPortfolio(code, dryRun)` | Runs one portfolio's program feeds, then its work-unit and compound feeds. With `dryRun` it prints the plan without applying it. |
| `syncDuePortfolios(codes)` | Syncs only the portfolios whose schedule is due, records each run, notifies owners of a run that changed something |
| `syncHistory(portfolio)` | A portfolio's sync runs as a table, newest first |
| `syncAllPortfolios` | Syncs discovery portfolios before development ones, so lineage resolves in one pass |
| `portfolioRoutingAudit` | For each source: keys no portfolio claims, and keys claimed twice |
| `routingAuditOf(sources)` | The same audit for sources handed over by another package |
| `buildWorkspace(code, dryRun, program)` | Builds one portfolio's workspace, or with `program` one program's folders. With `dryRun` it lists what it would create. Refused while another build of the same root runs. |
| `workspaceOrphans(code)` | The folders a prune-or-park pass would look at; a report, nothing is changed |
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
  - Isolation between portfolios, of program rows and of the workspace, run as a throwaway restricted user.
  - Model runs: the example model, a run saved at process lock (annotations, the stage folder, the reload, the `model_run` row by run id), and a gate document criterion selecting it.
  - Saving from the app: the dialog's rules (what it fills in, the program list, what it refuses, the result line) and the dialog itself, driven in the test: filled in for a Fed-batch titer run, a validity warning, saved to a program.
  - Scheduled sync: the rules (intervals and cron in UTC, when a sync is due, a run's counts and per-program changes, when owners hear of it), and on EX-SMD and EX-BIOD a recorded run, a due and a not-yet-due schedule, a platform notification, the function and the history panel.
  - Attention: the rules (what trips one, severity, messages, the lead set, what a rule needs), and on EX-SMD the flags of EXS-KIN1 and EXS-GPC3, an EX-BIOC lot, the function, the portfolio table's count, and a flagged compound opened.
  - Program data: the rules (the dataset list, what a further dataset needs, the table's name, default viewers), and on EX-SMD the SAR table and assay results of EXS-KIN1, what is refused, and the table view.
  - Choosing a run: the rules (what a criterion takes, the filter, who may choose, the comparison, the history), and on EX-BIOC the candidates of D5's process description, two runs compared, a choice superseded, and the dialog.
  - The model catalogue: its rules (codes, cards, merging across portfolios, validity ranges), the example catalogue loaded twice with nothing rewritten, the card joined with the function, the card rendered, and a run saved with the catalogue's codes.
- **`npm run test:bdd`: Gherkin scenarios** in `bdd/features/`, compiled to Playwright by
  `@datagrok-libraries/bdd`. Features address controls only by their `data-u2-name`.

On a bleeding-edge stand dated 2026-10-02, 162 of the 165 package tests and all 19 BDD scenarios pass.
The two isolation tests fail because of a platform gap, not a portal bug:
- The server refuses `permissions.grant` on a promoted domain row, even from the row's author with Share (`~can_share` true).
- `privileges_service.savePermission` checks Share against the projects repository, which never holds domain rows.
- Programs therefore stay private to their author, which fails closed. The tests stay in place as the acceptance check for the fix.

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
