# Project Portal Core

Several portfolios on one Datagrok instance. One organization can run small-molecule
discovery, biologics discovery and biologics development side by side. Each portfolio has its own
program registry feed, its own stage-gate model, and its own access groups.

## The model

**The portfolio is the unit of configuration.** Everything that varies between portfolios
hangs off the `portfolio` row. The alternative, one profile per instance, breaks as soon as
one organization runs discovery and development together.

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
  - `buildWorkspace` creates a Spaces tree: the portfolio's `workspaceRoot`, a folder per program named `CODE · Name`, and the portfolio's `workstreams` under each program. Stages add their folders as a program reaches them; a portfolio may add unit-operation folders from its process flows.
  - Compounds, plates, results and predictions are never folders. The builder is additive: it never renames, moves or deletes. Orphan folders are reported, and parked only when an approver asks.

## Features

Each feature is described in detail in [docs/features.md](docs/features.md):

- **Portfolios and sync:** [Workspace per portfolio](docs/features.md#workspace-per-portfolio), [Work units](docs/features.md#work-units), [Compounds in programs](docs/features.md#compounds-in-programs), [Scheduled sync](docs/features.md#scheduled-sync), [The portfolio table](docs/features.md#the-portfolio-table).
- **Gates and evidence:** [Data conditions and readiness](docs/features.md#data-conditions-and-readiness), [Stage advance driven by data](docs/features.md#stage-advance-driven-by-data), [Decisions with evidence](docs/features.md#decisions-with-evidence), [Program snapshots](docs/features.md#program-snapshots), [Documents](docs/features.md#documents), [Evidence tables](docs/features.md#evidence-tables), [Approval workflow](docs/features.md#approval-workflow), [Provenance](docs/features.md#provenance), [Data package manifest](docs/features.md#data-package-manifest).
- **Models and runs:** [Model runs](docs/features.md#model-runs), [Saving a run from the app](docs/features.md#saving-a-run-from-the-app), [Model catalogue and model cards](docs/features.md#model-catalogue-and-model-cards), [Choosing the run a criterion rests on](docs/features.md#choosing-the-run-a-criterion-rests-on), [Comparing a criterion's runs in Compute2](docs/features.md#comparing-a-criterions-runs-in-compute2), [Model recommendation](docs/features.md#model-recommendation), [Predictions](docs/features.md#predictions).
- **Discovery:** [Program data](docs/features.md#program-data), [Compounds needing attention](docs/features.md#compounds-needing-attention), [Observations](docs/features.md#observations), [Design sets and selection sets](docs/features.md#design-sets-and-selection-sets), [DMTA status](docs/features.md#dmta-status).
- **People:** [Modes and quick actions per persona](docs/features.md#modes-and-quick-actions-per-persona), [Program roles](docs/features.md#program-roles), [Write permissions on portal tables](docs/features.md#write-permissions-on-portal-tables), [Following a program: new-data notification](docs/features.md#following-a-program-new-data-notification), [Pinned programs: My programs on Home](docs/features.md#pinned-programs-my-programs-on-home), [Activity feed](docs/features.md#activity-feed).

## Design docs

The design behind the package is in [docs/](docs/README.md): why the portfolio is the unit of configuration
([multi-portfolio-design.md](docs/multi-portfolio-design.md)), the [architecture](docs/architecture.md), the
[schema notes](docs/schema.md), the [build log](docs/build-log.md) with the decisions open to revisit, and the
[platform findings](docs/platform-findings.md).

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
- a handoff that doesn't leave from a terminal stage, or goes to a portfolio that isn't configured;
- a flow step naming a unit operation that no catalogue of the set declares.

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
| `securePortalTables(dryRun)` | Server write permissions on `program_role`, `parked_folder`, `approval_request`: Edit to each portfolio's approvers (contributors for requests), none to All users. Admins; idempotent; preview by default |
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
| `dataPackage(portfolio, program)` | A program's data package manifest: criteria to the next gate, documents, evidence tables, runs, datasets, snapshots, with status and source |
| `dataPackageCsv(portfolio, program)` | The same as CSV text |
| `addObservation(portfolio, program, subject, code, title, finding, evidenceKind, evidenceTarget)` | Records an observation on the program, a series, campaign, study or compound, with one piece of evidence; contributors and approvers |
| `retractObservation(id, reason)` | Retracts an observation, kept and marked; its author or an approver |
| `observations(portfolio, program, subject, code)` | A program's observations as a table, current first |
| `createDesignSet(portfolio, program, series, name, source, proposals, scoreColumn, scoreName, lowerIsBetter)` | Proposed structures for a series from a table, ranked by score; contributors and approvers |
| `designSets(portfolio, program)` | The program's design proposals by set and rank, with what was chosen for synthesis |
| `chooseForSynthesis(portfolio, program, designSet, labels, chosen)` | Marks proposals chosen for synthesis, or not |
| `createSelectionSet(portfolio, program, name, campaign, rationale, compounds)` | A hit list of the program's compounds |
| `selectionSets(portfolio, program)` | The program's hit lists |
| `loadExampleDesignSets()` | The synthetic KIN-A design sets and K1 hit list on EXS-KIN1. EX-* only; idempotent |
| `dmtaStatus(portfolio, program)` | Each compound's DMTA stage with dates, design label, series, cycle days and sources |
| `dmtaCycles(portfolio, program)` | Per series: compounds at each DMTA stage and the median cycle time |
| `criterionRunsToCompare(portfolio, program, criterion)` | A run criterion's candidate runs by model, with run ids: what to add in Compute2's Compare Runs |
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
| `workspaceModelCard(model)` | The card beside a model's preview in the Workspace tab on Home (`meta.role: workspaceModelPanel`); nothing for an uncatalogued function |
| `loadModelCatalogue(catalogue, portfolio)` | Loads a catalogue handed over without a portfolio file, checked first. Idempotent |
| `buildAllWorkspaces` | Builds the workspace of every portfolio that declares a `workspaceRoot` |

## Tests

- **`grok test`: package tests.** Pure rules for every feature, and integration tests on the example portfolios.
  Coverage is listed in [docs/testing.md](docs/testing.md).
- **`npm run test:bdd`: Gherkin scenarios** in `bdd/features/`, compiled to Playwright by
  `@datagrok-libraries/bdd`. Features address controls only by their `data-u2-name`. Run with `--workers=1`.

Three package tests fail until the platform is fixed, and stay as the acceptance checks: the two isolation tests
(GROK-21041, a domain row cannot be shared) and saving a model run with its soft reference (GROK-21042).
`saveModelRun` falls back to `run_id` meanwhile.

## Known platform issues

The platform issues found building this, and the workaround for each, are listed in
[docs/platform-findings.md](docs/platform-findings.md).
