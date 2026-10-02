# Project Portal Core changelog

## 0.3.0 (unreleased)

- Activity feed (F-21, FR-80/81/83): a program's activity, newest first, read from the portal's own records: gate decisions (by an approver, the registry or the data) with their evidence, evidence added later, runs saved and selected, documents linked and approved, snapshots, compound role changes, registry and work-unit changes the sync saw, and new data. No event store. Pure entries, merge, filter and pages in `src/model/feed.ts`
- The program page's *Activity*: a kind filter, a series filter (the series' "what is new": its work units, its news, its compounds), pages of eight, Provenance on an entry over a record; **Activity** on each series row. The ribbon's **My feed**: everything the user follows, a followed series narrowing its program to the series
- `programFeed(portfolio, program, series, kinds, limit)` and `myFeed(kinds, limit)` return it as a table
- The sync keeps what followers were told in `sync_run.changes.newData`; a follow keeps when its last news was told (`follow.news_on`)
- `portal` schema 0.2.21: `follow.news_on`. Additive, no migration script

- Provenance (F-28, FR-75): from a program, decision, piece of evidence, model run, snapshot, document or compound, the portal's records as a graph walked both ways: what it rests on (the program data and the compound in it, a snapshot of the data, a run's input tables and model, the run selected for a criterion, the document answering it, the criterion, every piece of a decision's evidence) and what rests on it (a run's output tables, documents, decisions, the program a decision moved). Pure graph and walk in `src/model/provenance.ts`
- The provenance dialog: *Rests on* and *Used by* as indented trees; each record opens (a run's call, a snapshot, a document's file, the program data, a model card, a table, a program on its page) and **Trace** starts again from it. **Provenance** on the program page, on each decision and piece of evidence in the history, and on each saved run, document and snapshot
- `provenance(kind, id, depth)` returns `{start, nodes, edges}`; `showProvenance(kind, id)` opens the dialog

- Documents over files (F-27, FR-70): a `document` record over one file, in the program's workspace folder or at any file path the platform reads, with its kind, title, version, status (draft, approved, superseded), author, and the gate criterion of kind document it answers. Approving one meets that criterion unless it was waived; linking a draft puts a criterion not started in progress; a newer approved version of the same kind and criterion supersedes the older
- Document cards on the program page (*Documents*): one per document criterion of the current and the next gate, with its state, the documents over it (each opens; a draft has **Approve** for approvers), what is missing, the run selected for it and the decisions that cite it. **Link document…** links a file of the program folder, a path, or an upload into the folder (no overwrite). Documents no card shows are listed below
- Documents are evidence of kind `file` with the record in the context: the gate panel suggests the approved documents of its criteria, `addEvidence` takes a document by id or path, and the history looks a document up by its record
- `programDocuments(portfolio, program)`, `linkDocument(portfolio, program, path, criterion, kind, title, version, status)` and `loadExampleDocuments()`. Optional portfolio field `documentKinds`; without it, one kind per document criterion, and `other`. `exportPortfolio` writes it back; `diffPortfolios` reports it as vocabulary
- EX-BIOC declares four document kinds and ships four synthetic documents (process description draft for EXB-201, an approved development report, an IND quality draft for EXB-202, an approved transfer package for EXB-104)
- `portal` schema 0.2.20: table `document`, `portfolio.document_kinds`. Additive, no migration script

- Dated program snapshot (F-32, FR-15): `snapshotProgram(portfolio, program, label, note)` freezes the program's stage and status, the state of every criterion that applies (the next gate's data conditions as evaluated), the runs selected, and the program dataset (the SAR table on EX-SMD) as a CSV file in the program's workspace folder, its SHA-256 on the record. Written once; a file changed afterwards shows as changed. Without a workspace folder the table's size and fingerprint are kept, not the table
- The program page lists snapshots (*Snapshots*, behind More in Simple mode) with **Snapshot now**; each opens read-only: stage, criteria, runs, and the table in a view with editing off. `programSnapshots(portfolio, program)` lists them
- A decision may take a snapshot now as evidence (**Snapshot now** in the gate panel, `snapshot` on `recordDecision`); the history's snapshot link opens it. `addEvidence` takes a snapshot by id
- `portal` schema 0.2.19: table `program_snapshot`. Additive, no migration script

- Decisions with evidence (F-25, FR-71/74): a gate decision (the gate-log entry) carries evidence rows: a model run by `run_id`, a file path, a URL, a criterion's state at the gate, and a program snapshot. Each points at its target by a soft reference and freezes the name and context when attached, so the history reads the same after the target is renamed or deleted, and says so (`now …`, `gone; as recorded`)
- The gate panel suggests the gate's readiness and the runs selected for its criteria, all ticked; the approver may untick any and attach a link or file path. The program's history lists each decision with its evidence, each piece clickable (a run opens its call, a criterion its stage, a URL its page)
- `recordDecision` writes the suggestion when no list is given, so an advance by data now records its readiness evidence the same way as an approver's. `stage_transition.evidence` is still written
- `decisionEvidence(portfolio, program, decision)` lists decisions with their evidence and whether each target is live, renamed or gone; `addEvidence(portfolio, program, kind, target, name, decision)` adds a piece to a decision already recorded (default the latest), approvers only
- `portal` schema 0.2.18: table `evidence`. Additive, no migration script

- Program roles (F-16, app side): viewers, contributors and approvers per program, on top of the portfolio's groups, held by a group or a user. Managed by the portfolio's approvers on the program page (*Program roles*) and with `setProgramRole`; listed by `programRoles`. A program approver records gate decisions and a program contributor selects runs, checked where the portfolio's groups are checked. The program's workspace folder is shared to match (View to viewers, Edit to contributors and approvers), also by the builder. Row grants wait on GROK-21041; `ROW_GRANTS` in `src/program-roles.ts` turns them on
- `portal` schema 0.2.17: table `program_role`. Additive, no migration script
- Simple mode keeps Program roles behind More

- Simple and Advanced mode, quick actions per persona (F-08, FR-30–32): the board's Mode and Persona selects, remembered per user in the platform's user settings. Simple keeps missing assays, suggested models, attention, saved runs and history behind More, and drops Sync history, Export and All portfolios from the ribbon; Advanced (the default) shows everything. Quick-action cards for chemist, biologist, modeler, process scientist and lead head the program page, each opening an existing action (open data, needs attention, missing assays, suggested models, save a run, choose a run, follow, open workspace, gate decision, portfolio table); one that cannot be used says why
- `personas` (optional) in a portfolio file: groups mapped to personas, a default persona, and cards per persona, including function cards (`{name, function}`) that get the program's portfolio, program and stage. Without it the persona comes from the user's choice, else the setting
- `portal` schema 0.2.16: `portfolio.personas`. Additive, no migration script
- EX-SMD and EX-BIOC map their approvers to the lead persona
- `programData(portfolio, program, dataset)` is now registered under that name, as documented; `programDataTable`, the name it was registered under by mistake, stays as a second function

- New-data notification (F-11, FR-83): a user follows a program, or one of its series, from the program page (**Follow**, and Follow on each series row) or with `followProgram`/`unfollowProgram`. After every sync, and on `checkNewData(code)`, each followed program's dataset is read once and compared with each follow's snapshot (rows, filled values, and the largest value of `programData.watermark` when declared); followers with news get one notification per portfolio in the platform's notifications. Only growth is news; following takes a baseline. The program page shows *What's new since your last visit* to a follower
- `programData.watermark` (optional): a column that only grows, so a replaced row counts as new
- `portal` schema 0.2.15: table `follow`. Additive, no migration script

- Model recommendation v1 (F-13, FR-51): a card may declare `appliesTo` (tasks, stages, modalities, work-unit kinds, the program-data columns it predicts, rank; a key left out matches any). The program page shows *Missing assays* for the gate it shows (the columns open data conditions need and rows lack, and which suggested model predicts each) and *Suggested models*, best first, each with **Run**, which opens the model's dialog filled in from the program (the program dataset, program and series codes). `recommendModels(portfolio, program, stage, series, task)` and `missingAssays(portfolio, program, stage, series)` return them
- `portal` schema 0.2.14: `model_card.applies_to`. Additive, no migration script. `exportPortfolio` writes it back
- EX-SMD declares two synthetic predictive models, `SolubilityEstimate` and `HergRisk` (tiny scripts, invented numbers)

- Profile tooling (F-18): `exportPortfolio(code)` reads a loaded portfolio back into the file format (round-trips the three example files); `diffPortfolios(a, b)` compares loaded portfolios and files in readable lines (vocabulary, workspace tree, other fields); `exportProfile(codes)` and `importProfile(bundle, dryRun)` move a set of portfolios between instances: validated as a set, diffed against what is loaded, refused when invalid, loaded unless dry run. The board's **Export** shows the open portfolio's file
- `validatePortfolioSet` and `loadPortfolioSet` take the codes already loaded, so a handoff may name one of them
- `portal` schema 0.2.13: `portfolio.handoff` (was validated, not stored), `unit_operation.portfolios`, `run_type.portfolios`. Additive, no migration script

- Stage advance driven by data (opt-in, off by default): a portal-authority portfolio may declare `autoAdvance: {through, countWaivers}`. After each sync, and on `autoAdvancePortfolio(code, dryRun)`, an active program whose next gate is met from its data alone moves one stage, never past `through`; a waiver counts only with `countWaivers`, a criterion recorded met by hand never does, held, ended and unstaged programs stay. The move is recorded like an approver's advance (gate log, criterion evidence, workspace folders) with source `data` and the evidence. The sync's report, history and owner notification name each move; Sync history has a preview
- The gate log records who moved a program: `approver`, `registry` or `data`
- `portal` schema 0.2.12: `portfolio.auto_advance`, `stage_transition.source`, `stage_transition.evidence`. Additive, no migration script
- EX-SMD lets data advance through Hit-to-Lead

- Scheduled sync (F-01): a portfolio may declare a `schedule` (an interval such as `every 6h`, or a cron line in UTC) and `owners` (user logins). `syncDuePortfolios(codes)` syncs only the portfolios whose schedule is due; a portfolio without one is never touched
- Sync history: every sync, by schedule or by hand, is recorded (started, finished, trigger, who, counts of created, updated and flagged rows, unresolved spellings, error, report, per-program changes). The board's **Sync history** shows schedule, next slot, owners and runs; `syncHistory(portfolio)` returns them. `syncPortfolio` returns the same report as before
- Owners are notified in the platform's notifications after a scheduled run that changed something, failed, or found new standing issues
- `portal` schema 0.2.11: `portfolio.schedule`, `portfolio.owners`, table `sync_run`. Additive, no migration script
- The example portfolios declare schedules and an owner. The package cannot register the schedule itself (the platform schedules only server-side scripts and queries)
- Compounds needing attention (F-17, FR-104): a portfolio declares `attention` rules (code, name, severity, condition, message) in the criteria's condition language. A row trips a rule when the scope admits it, `where` holds and, if the rule names `measured` columns, it lacks one. The program page shows *Needs attention*: a count per rule and the flagged compounds (or lots); a click opens the program data at that row. The portfolio table has an `attention` column (a filter) and `attention rules`. `programAttention(portfolio, program)` returns the flags
- `portal` schema 0.2.10: `portfolio.attention`. Additive, no migration script
- EX-SMD flags potent compounds without hERG, leads without solubility and a hERG IC50 under 10 µM; EX-BIOC flags a tox lot not released
- Program data one click away (F-09): the program page's **Open data** opens the portfolio's program dataset for the program as a table view named after it, with the viewers the portfolio declares or the defaults for its setting (discovery: two measures plotted against each other, coloured by the first category, and a histogram; development: a bar chart over the first category, and a histogram). `programData(portfolio, program, dataset)` returns the table; `openProgramData(...)` opens it. The query runs as the user, on a connection shared with them
- `programData` takes an optional `name` and `viewers`; a portfolio may declare further datasets per program (`programDatasets`: code, name, connection, SQL with `@program`, viewers), each with its own button
- `portal` schema 0.2.9: `portfolio.program_datasets`. Additive, no migration script
- EX-SMD calls its program dataset the SAR table and declares Assay results beside it
- Choose the run a criterion rests on (F-15): a criterion may declare `artifact: "model_run"` and a `runFilter` (unit operations, run types, models by code). The gate panel shows the run selected for it and the ones it superseded, and **Choose run…** opens a dialog over the program's runs the filter admits: tick two or more to see their inputs and outputs side by side (what differs in bold), pick one, say why. `selectRunForCriterion` keeps its behaviour (a later run supersedes; the criterion is met unless waived) and now refuses a user who is neither contributor nor approver, a program that has ended, and a blank rationale. `criterionRuns` and `runComparison` return what the panel and the dialog show
- `portal` schema 0.2.8: `gate_criterion.artifact`, `gate_criterion.run_filter`. Additive, no migration script
- EX-BIOC: D5 `process_description` rests on a production bioreactor run (`EX-PROD`)
- Save to program from the app (F-14, FR-53): a dialog asks for portfolio, program (the ones you can see in that portfolio), unit operation, run type, purpose and title, fills in what it can (the program page's portfolio and program, the call's earlier annotations, the model's catalogue portfolio, unit operation and run type), warns when an input is outside the card's validity range, then calls `saveModelRun` and says where the run was linked. It opens from the balloon that follows each run of a catalogued model, from a run's context menu, from the program page's Save a run button (this session's runs), and as `saveRunToProgram(call)`
- `buildWorkspace` takes an optional `program`, to build one program's folders
- Model catalogue and model cards (F-12, FR-50): a portfolio file declares its `catalogue` (unit operations, run types, models with their cards), checked with the set and loaded with the portfolios, idempotently; `loadModelCatalogue` loads one handed over without a portfolio file. Runs saved before their model was catalogued are linked to it on load
- Model card: the question the model answers, intended use, inputs with their validity range, outputs, applicability, performance, lineage, limitations, owner, approval state. Inputs and outputs the card does not list come from the function's signature. `modelCard(model)` returns it; the program page opens it from a saved run
- Program page: the program's saved runs, newest first, with model, stage at run, unit operation, run type and purpose
- `portal` schema 0.2.7: `model_card` (one per model version); `model.current_version`, `model.portfolios`, `model.unit_operations`, `model.default_run_type_id`. Additive, no migration script
- Example catalogue for EX-BIOC: unit operations EX-SEED, EX-PROD, EX-HARV, EX-PROA, EX-UFDF; run types EX-SIM, EX-TWIN, EX-BENCH; `FedBatchTiter` with a draft card
- `saveModelRun` records the uploaded tables' ids in the call's options (`portalTables`), so removing a run can remove its tables
- Loading configuration compares json values with their keys sorted, so a reload no longer rewrites json columns the server returns in its own key order
- A program has a workspace only when its folder still exists: the portfolio table's workspace flag and the program page's Open workspace button check the server, not the recorded `space_id`, so folders deleted by hand no longer show as present
- Portfolio table: one row per program in a Datagrok table view with the filter panel (Table and All portfolios on the board; `/Portal/<code>/table`, `/Portal/table`). Stage in model order, days in stage, readiness for the next gate, compounds by role, work units, workspace. Default viewers: programs by stage and modality, days in stage, a readiness matrix and compounds by role, all driven by the filters. A row opens the program's page. `portfolioTable(portfolio)` returns the table
- Model runs: `saveModelRun` annotates a call (portfolio, program, stage at run, unit operation, run type, purpose) in its options, saves it, links it into the program's workspace (the current stage's folder when there is one) and records a `model_run`; `selectRunForCriterion` selects a run for a gate criterion and supersedes the one selected before
- Example model `FedBatchTiter` (JavaScript, `meta.role: model`): a synthetic fed-batch titer model for EX-BIOC
- `portal` schema 0.2.6: `model_run.run_id`, the saved call's id as text, always written; `model_run.func_call` is no longer required, because the server refuses the soft reference today. The insert tries it first and falls back to `run_id` alone (migration debt, see the README). Migration script `databases/portal/0001_model_run_run_id.sql` and its `down/` twin
- Data conditions are evaluated: a portfolio declares a program dataset (`programData`, SQL with `@program`), and a criterion a structured `condition` (scope program / lead set / series, `where`, `measured`, `minCount`, `minSeries`)
- Readiness: `programReadiness` evaluates every criterion of a stage and returns a score and a gap list; the gate uses the evaluated conditions and records what the data said on advance
- Criteria may be limited to a modality (`modality`)
- App: the gate panel shows readiness and the gaps; any stage on the track shows what entering it takes; the program page lists its work units and its compounds by role
- Work units: `campaign`, `series`, `study`, `process_stage` under a program; a portfolio enables a subset (`workUnits`) and feeds them (`workUnitFeeds`), synced after the programs; flagged, never deleted, when a feed drops them
- Workspace grows with the stage: stages declare `folders`; the builder adds those of every stage a program has reached, and an advance adds the new stage's at once (PORT-35)
- Orphan report: `workspaceOrphans`, and every full build, list root folders no program claims and folders of ended programs; nothing is pruned or parked (FR-22)
- Compounds in programs: `compound` (a registry projection, never minted or promoted), `program_compound` with a role, `compound_role_change` history; a portfolio's `compoundFeed` syncs them; `setCompoundRole` changes a role by hand and the registry no longer overwrites it
- `portal` schema 0.2.2: `gate_criterion.condition`, `portfolio.program_data`; 0.2.3: the work-unit tables, `portfolio.work_units`, `portfolio.series_maturity`, `registry_feed.target`; 0.2.4: `stage.folders`; 0.2.5: `compound`, `program_compound`, `compound_role_change`
- Examples: EX-SMD compounds with roles, batches and a tool compound in two projects; stage folders (EX-SMD: Biology/Assay Development, DMPK/In vivo PK, Candidate Package; EX-BIOC: Regulatory/IND, Drug Substance/Process Lock, Regulatory/BLA); series and campaigns (EX-SMD), campaigns (EX-BIOD), studies and process stages (EX-BIOC); EX-SMD compounds by series and role (hit-to-lead and lead-optimization conditions, an optional hERG margin); EX-BIOC lots for the tox-lot condition

## 0.2.0 (2026-09-30)

- Workspace builder: a Spaces tree per portfolio (root, a folder per program, its workstreams), additive and idempotent, with a preview; the root is shared with the portfolio's groups
- App: the address carries the open portfolio's code, and the gate panel opens the program's workspace
- `portal` schema 0.2.1: `portfolio.workstreams`, `program.space_id`
- Portfolio files take `workstreams`. The examples declare Chemistry, Biology, DMPK, CompChem (EX-SMD); Protein Engineering, Biology, Biophysics, Developability (EX-BIOD); Drug Substance, Drug Product, Analytical, Regulatory (EX-BIOC)

## 0.1.0 (2026-09-30)

- `portal` domain schema: stage models, stages, gate criteria, portfolios, registry feeds, programs, stage transitions, criterion status
- Portfolio files as configuration: groups, stage model, stage authority, registry feeds
- Sync planner with stage authority, lineage across portfolios, overlap and missing detection; routing audit
- Gate rules and the gate log
- Project Portal app (u2): portfolio list, per-model board, gate panel
- Package tests and BDD scenarios
