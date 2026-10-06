# Build log: decisions to revisit

The core was built feature by feature on the synthetic example portfolios. At each step the conservative choice
was made where the requirements left one open; this page lists those choices, so they can be confirmed or
changed. What each feature does is in [features.md](features.md); platform issues are in
[platform-findings.md](platform-findings.md).

| Feature | Ids | Schema |
|---|---|---|
| Readiness, work units, stage folders, compounds | F-02, F-04, F-20; FR-03/04, FR-06, FR-11, FR-13, FR-14, FR-22 (report) | 0.2.2 to 0.2.5 |
| Portfolio table, filters and analytics | F-06 | none |
| Model runs | F-14 (API), F-15 (API) | 0.2.6: `model_run.run_id`, migration `0001_model_run_run_id.sql` |
| Model catalogue and model cards | F-12, FR-50 | 0.2.7 |
| Save a run to a program from the app | F-14 (UI), FR-53 | none |
| Choose the run a criterion rests on | F-15 | 0.2.8: `gate_criterion.artifact`, `run_filter` |
| Program data one click away | F-09 | 0.2.9: `portfolio.program_datasets` |
| Compounds needing attention | F-17, FR-104 | 0.2.10 |
| Scheduled sync | F-01 | 0.2.11: `portfolio.schedule`, `portfolio.owners`, table `sync_run` |
| Stage advance driven by data | F-20 (auto-advance) | 0.2.12: `portfolio.auto_advance`, `stage_transition.source`, `stage_transition.evidence` |
| Profile tooling | F-18 | 0.2.13: `portfolio.handoff`, `unit_operation.portfolios`, `run_type.portfolios` |
| Model recommendation | F-13, FR-51 | 0.2.14: `model_card.applies_to` |
| New-data notification | F-11, FR-83 | 0.2.15: table `follow` |
| Modes and quick actions per persona | F-08, FR-30 to FR-32 | 0.2.16: `portfolio.personas` |
| Program roles | F-16 (app side), FR-90 | 0.2.17: table `program_role` |
| Decisions with evidence | F-25, FR-71, FR-74 (in part) | 0.2.18: table `evidence` |
| Dated program snapshot | F-32, FR-15 | 0.2.19: table `program_snapshot` |
| Documents over files | F-27, FR-70 | 0.2.20: table `document`, `portfolio.document_kinds` |
| Provenance panel | F-28, FR-75 | none |
| Activity feed | F-21, FR-80/81/83 | 0.2.21: `follow.news_on` |
| Predictions | F-22, FR-58/103 | 0.2.22: table `prediction` |
| Approval workflow | F-30, FR-57 | 0.2.23: table `approval_request`, `portfolio.approval` |
| Park orphan folders | FR-22 (action half) | 0.2.24: table `parked_folder` |
| Evidence tables generated into documents | F-40, FR-72 | 0.2.25: `document.generated` |
| Evidence changed since written | F-41, FR-73 | 0.2.26: `document.evidence_check` |
| Pinned programs on the landing page | pinning, My programs widget | none |
| Model card in the Workspace preview | Power Pack hook `workspaceModelPanel` | none |
| Quality pass | review of the features above | none |
| Unit-operation folders from a program's flow | F-05 leaf level | 0.2.27: `portfolio.flow_folders` |
| Data package manifest | F-43, FR-43, PORT-11 (as a view) | none |
| Cross-model run comparison | F-29, FR-55, PORT-48 | none |
| Observations | F-26, FR-76/106, PORT-91 | 0.2.28: table `observation` |
| Design sets and selection sets | F-23, FR-101/102, PORT-85/89 | 0.2.29: tables `design_set`, `design_proposal`, `selection_set` |
| DMTA status | F-24, FR-105, PORT-86 | 0.2.30: `portfolio.dmta` |
| Program-page panels per persona | F-31, FR-34, PORT-09 | none (`portfolio.personas`) |
| Write permissions on portal tables | security | 0.2.31: three tables View for all users; metadata only |
| Three small fixes | user names, no folders for ended programs | none |

## Readiness, work units, stage folders, compounds

*F-02, F-04, F-20; FR-03/04, FR-06, FR-11, FR-13, FR-14, FR-22 (report). Schema: 0.2.2 to 0.2.5.*

- **A compound's role:** the registry owns it until a person changes it. After that the link is manual and the registry's value is only reported. This mirrors `program_flow`.
- **A recorded "met" on a data condition** is overruled by the data. Only a waiver opens a gate the data closes.
- **Compound rows** are readable by all users. Which programs hold a compound is under each program's security.
- **Compound sync** writes row by row. That's fine for lead and candidate sets, but not for a 10⁶-compound registry (NFR-01); bulk binding waits on PLAT-13.
- **Series scope** reads `series` and `series_maturity` from the program dataset, not from the `series` table. Joining the two is open.

## Portfolio table, filters and analytics

*F-06. Schema: none.*

- **Isolation is untested in practice.** The restricted-user test passes, but because program rows cannot be shared (GROK-21041) a portfolio viewer sees no programs at all. It proves that an outsider sees nothing, not yet that a viewer sees their own portfolio's rows.
- **Workstreams are not a column.** They are a portfolio-level list, not a property of a program.
- **The table does not refresh by itself.** After a gate decision made from the row panel, reopen the table to see the change.
- **All portfolios shows every portfolio on the instance,** those loaded by other packages included.

## Model runs

*F-14 (API), F-15 (API). Schema: 0.2.6: `model_run.run_id`, migration `0001_model_run_run_id.sql`.*

- Selecting a run marks a document criterion met.
- The default folder is the current stage's first declared folder. A stage with no folders (for example D4) puts the run in the program folder, not the last folder a previous stage declared.
- Unknown unit operation or run type codes are an error. `model_id` is set only when the model is already in the `model` catalogue (keyed by nqName). Nothing is auto-registered.
- The FuncCall has no ACL (platform), so anyone with the id can load a run. Only its uploaded tables are shared with the portfolio's viewers.

## Model catalogue and model cards

*F-12, FR-50. Schema: 0.2.7.*

- The catalogue sits in the portfolio file, not in a file of its own. Unit operations and run types are global tables, so the example codes carry an `EX-` prefix.
- A model a second portfolio declares must be declared identically; there is no per-portfolio override of a card.
- Loading never removes a portfolio, unit operation or card from a model; retiring is done by setting `approval: retired` or `status: deprecated`.
- The example card is `draft` and says its numbers are invented.
- Earlier runs on EX-BIOC (EXB-201, EXB-203) were linked to the catalogued model when the catalogue loaded.

## Save a run to a program from the app

*F-14 (UI), FR-53. Schema: none.*

- The balloon appears after every run of a catalogued model, scripted runs included (package tests show it too). It needs the package loaded first, through an `init` function, not `autostart`, so the 1 MB bundle is not loaded for every user at start-up.
- The program page offers only this session's runs; runs from earlier sessions are not on the server unless saved.
- Purpose is optional, as in `saveModelRun`.
- An input outside the validity range warns, never blocks.

## Choose the run a criterion rests on

*F-15. Schema: 0.2.8: `gate_criterion.artifact`, `run_filter`.*

- Selecting a run still marks a document criterion met (kept as it was).
- Contributors may select, not only approvers, and on feed-authority portfolios too: a run is evidence, not a stage decision. The Met and Waive buttons still follow the old rule (approvers, portal authority only).
- The filter decides what the dialog offers. `selectRunForCriterion` does not enforce it, so a script can still select any run of the program.
- Filter codes are not checked against the catalogue; an unknown code just matches nothing.
- A saved call holds a table output as an id, so the comparison shows "saved table", not its size.

## Program data one click away

*F-09. Schema: 0.2.9: `portfolio.program_datasets`.*

- Default viewers depend only on the portfolio's setting, not its modality: no structure-aware viewers (no molecule column in the examples).
- A further dataset is per portfolio, the same for every program. No per-program or per-stage datasets.
- The view is not saved anywhere; each click runs the query again and opens a new view.
- The dataset code `program` is reserved for the program dataset.

## Compounds needing attention

*F-17, FR-104. Schema: 0.2.10.*

- "Opens the compound" opens the program dataset at the compound's row. There is no compound page yet, and compounds are never promoted to entities.
- An item is named by the row's `compound`, `id`, `code`, `lot` or `name` column, the same as the gap list. Two rows with the same compound count as one item.
- Rules are per portfolio, not per stage or modality. Every program of the portfolio is checked, ended ones too.
- No history: a flag is computed on each read, never stored. Nobody is notified when a new flag appears (that is the new-data notification, FR-83).
- Severity only orders and colours; nothing blocks on it.

## Scheduled sync

*F-01. Schema: 0.2.11: `portfolio.schedule`, `portfolio.owners`, table `sync_run`.*

- Cron is read in UTC, not the server's or the owner's time zone.
- A manual sync counts as the last run, so it pushes the next scheduled slot back. A failed run waits for the next slot instead of retrying.
- Only scheduled runs notify; a manual sync never does. A quiet run that finds the same unresolved spellings and overlaps as the last one does not notify.
- Owners are user logins, not groups. A run by a user without Edit users can notify only themselves; the others are listed as not notified on the run.
- `sync_run` is readable by all users and written by whoever syncs (today: admins). `syncAllPortfolios` records manual runs for every loaded portfolio, the customer ones too (history rows only).
- "Flagged" counts rows newly flagged as missing from their feed. Overlaps are standing issues, shown but not counted.

## Stage advance driven by data

*F-20 (auto-advance). Schema: 0.2.12: `portfolio.auto_advance`, `stage_transition.source`, `stage_transition.evidence`.*

- Entering the model (unstaged → first stage) is never by data; a person places the program.
- Only evaluated data conditions count. A data condition with only a free-text rule (EX-SMD's Hit Finding assay criterion) never advances by data, even when recorded met.
- A stage without required criteria is not entered by data ("nothing for data to decide"), although its gate is trivially met.
- One stage per sync, even when the data would carry the program further within `through`.
- The sync applies it as whoever runs the sync, and only if that user is an approver (admins are). Otherwise the run reports "not applied".
- `autoAdvancePortfolio` defaults to a dry run, unlike `syncPortfolio`.
- `through` is a stage code of the portfolio's model; no per-modality limit.

## Profile tooling

*F-18. Schema: 0.2.13: `portfolio.handoff`, `unit_operation.portfolios`, `run_type.portfolios`.*

- Loading never deletes, so an import cannot remove a stage, criterion, feed, folder or catalogue entry. The diff shows removals and the report says they are not applied. Retiring needs its own step.
- A model shared by two portfolios exports with the union of their unit operations, as stored.
- Rows loaded before 0.2.13 have no `portfolios` on unit operations and run types; export then falls back to those the portfolio's models use. A reload fills it.
- A bundle carries no instance identity or schema version beyond its own format version; connections are named, so the target instance must have connections of the same names.
- Groups are exported by name; import creates missing groups as loading always has.
- No guard limits import to EX-* codes in code; the restriction applied to testing on the stand only.

## Model recommendation

*F-13, FR-51. Schema: 0.2.14: `model_card.applies_to`.*

- `appliesTo` is on the card (per version), not a separate `model_applicability` table. One row per model version, keys as lists; no per-key weights.
- A model with no `appliesTo` is suggested for every program of the portfolios that catalogue it. A model of another portfolio is suggested only if its card declares where it applies.
- Stage matches the program's current stage or the stage the gate panel shows. An unstaged program shows nothing for stage-keyed models until a stage is clicked.
- Tasks are free words, not a checked vocabulary. Without a task asked, the task does not filter.
- The program page shows program-level gaps only; per-series gaps are in `missingAssays(..., series)`.
- A criterion already met is never a gap, even when some rows lack a column. A criterion not met with every value measured (potency too weak) is not a missing assay either.
- Run prefill is by input type and name only; it never maps card inputs to dataset columns.

## New-data notification

*F-11, FR-83. Schema: 0.2.15: table `follow`.*

- Counts, not a diff of rows: an edited value, a deleted row, or a row replaced by another without a watermark is not reported. A declared watermark is the robust option; none of the examples declares one.
- Only the program dataset is watched, not the further datasets (Assay results has a `run_date` that would make a good watermark).
- Manual syncs notify followers too (owners hear only of scheduled runs). A sync run by a user without Edit users reaches only that user; the rest are listed as not notified.
- Only a series follow narrows the rows (by the `series` column). A campaign, study or process-stage follow could be stored but hears about the whole program; the UI offers follow on series only.
- `follow` is table-secured: any user may read and edit any follow row. Row security would hide follows from the sync that reads them.
- Unfollowing turns a row off instead of deleting it (platform issue below).
- What's new is shown only to followers, and the visit is recorded on page open, not when the user has read it.

## Modes and quick actions per persona

*F-08, FR-30 to FR-32. Schema: 0.2.16: `portfolio.personas`.*

- Advanced is the default, so nothing changes for a user until they choose; FR-31 may want Simple first for new users.
- Mode and persona are per user, not per portfolio. A user in two portfolios with different jobs switches by hand; until they choose, each portfolio's group mapping decides.
- Simple hides sections in the app only. FR-31's permission trimming per persona group is not done: the platform's own tree is unchanged (PLAT-17).
- Work units and compounds stay visible in Simple; the program's saved runs are behind More.
- The core's cards per persona are a first guess; a portfolio overrides them with `personas.cards`.
- Function cards run as the user with no confirmation; the portfolio file is trusted configuration.

## Program roles

*F-16 (app side), FR-90. Schema: 0.2.17: table `program_role`.*

- `program_role` is table-secured, View and Edit for all users, like `follow`: the app checks who may change it, the server does not. Master security on `program` would hide roles from everyone but the program's author until GROK-21041 is fixed.
- One role per program and group. Removing a role turns its row off; the server could not delete rows of a table without soft delete (see new-data notification).
- Only portfolio approvers manage roles; a program approver does not, so a role cannot grant itself more.
- Program viewers get nothing in the app beyond the folder today: reading the program row needs the row grant.
- Folder sharing runs for any portfolio that uses roles, not only EX-*; on this stand it was exercised on EX-SMD only.
- Groups are not created per program (FR-90 says "created with the program"): a role names an existing group or user. Creating `<program> Viewers` etc. on sync is a later step.
- Compound role changes (`setCompoundRole`) have no access check today, at either level; not added here.

## Decisions with evidence

*F-25, FR-71, FR-74 (in part). Schema: 0.2.18: table `evidence`.*

- No `decision` table: FR-74's `progress / drop / nominate / select` types and the work-unit link are not modelled; only gate decisions carry evidence. A document (FR-71's other parent) has none either.
- The target is plain text, not a soft-reference column: the server refuses soft references to `Core` tables today (the `model_run.func_call` bug), and a URL or file path has no table.
- Evidence can be added to a decision afterwards (`addEvidence`), dated when added. The gate log itself stays write-once.
- Readiness is one tick in the panel for all required criteria of the gate, one row per criterion when written. Optional criteria are left out.
- A decision without a list (API, data) gets the suggestion; an empty list gets none.
- A file target is checked with `dapi.files.exists`; a space's own files may report as unchecked.
- `evidence` is master-secured through the program, like the gate log, so it waits on GROK-21041 for anyone but the program's author too.

## Dated program snapshot

*F-32, FR-15. Schema: 0.2.19: table `program_snapshot`.*

- The table is a file in the program folder, not an uploaded table entity. Anyone with Edit on the folder can overwrite it; the SHA-256 catches that, it does not prevent it.
- Files sit at the program folder's top level. A `Snapshots` subfolder would show up in the workspace builder's orphan report.
- Only the program dataset is frozen. Further datasets (`programDatasets`, e.g. assay results) are not.
- Without a workspace folder the table is not kept anywhere. Embedding it in the row was ruled out because of size.
- Snapshots are master-secured through the program, like the gate log.
- The label defaults to `<program> <time>` (plus "at the advance decision" when a decision takes it). The page shows the label and a summary line.
- Viewers cannot take snapshots.

## Documents over files

*F-27, FR-70. Schema: 0.2.20: table `document`, `portfolio.document_kinds`.*

- Approving a document writes `criterion_status` met (unless waived), rather than readiness computing it from documents. Gate rules, snapshots and the portfolio table see it without change. A superseded document never un-meets the criterion; nor does deleting one.
- Only approvers approve; contributors link drafts. Approving is allowed on a feed-authority portfolio (EX-BIOC), where gate decisions are not.
- A document is evidence of kind `file` with `context.document`, not a new kind, so the `evidence.kind` choices did not change.
- One document per file path per program: linking the same path again changes that document. A new version is a new file.
- Uploads go to the top of the program's folder, never overwrite, and need the workspace folder. The file is not fingerprinted (snapshots are).
- The example documents are package files linked by path, not copies in the program folders.
- `document` is master-secured through the program, like the gate log, so it waits on GROK-21041 for anyone but the program's author.
- The column is `doc_version`: `version` is a system column of every domain table.

## Provenance panel

*F-28, FR-75. Schema: none.*

- Evidence is an edge, not a node. Starting from a piece of evidence starts from its target.
- FR-75's "prediction" is a run's output table. Runs are not linked to the program data: a saved call names its own input tables only.
- A compound is linked only to its program's dataset, and only when the walk starts from that compound, so a program's walk does not list every compound.
- The program is reached from its decisions only. Runs, documents and snapshots that no decision cites do not lead to the program.
- Saved calls are read for the 40 newest runs of the programs read; older runs show without their tables. A compound's walk reads at most 10 programs.
- Default depth 4 each way. Nothing is cached; each open reads the records again.
- A click opens the record. A decision or compound has no view of its own, so a click traces from it; a criterion click does nothing.

## Activity feed

*F-21, FR-80/81/83. Schema: 0.2.21: `follow.news_on`.*

- Read on demand, not stored. Each open of a program page reads ten record kinds; My feed reads them per followed program. Fine for tens of programs; a feed table or the platform's audit (PLAT-03/06) when it is not.
- Only the last 50 sync runs of a portfolio are read, and only the last news of each follow outside a sync. Older standalone checks are not in the feed.
- Attention count changes are not in the feed: the portal keeps no history of them. A sync could record the counts per program the way it records new data.
- Model-card approvals, program-role changes and follows are not in the feed (role rows have no history; a card has only its approval date).
- Entries are not filtered by what the reader may see beyond the program's own row security: master-secured records follow the program, `sync_run` and `follow` are table-secured, so a viewer of the program sees its sync and news lines.
- A series narrows by tag and compound only. A run, document or decision is never tagged with a series, so it drops out of a series' view.

## Predictions

*F-22, FR-58/103. Schema: 0.2.22: table `prediction`.*

- `prediction` is master-secured through the program, not a promoted entity with lazy promotion (FR-58). Enough while predictions are only read on the program page.
- One run per scoring, even when only a few compounds are new; the run scores the whole dataset and only new compounds are stored. Values already stored for the same version are never rewritten, even if the model would now say otherwise.
- The model version is the card's version, not the package version `model_run.model_version` records.
- The estimate column is found by name (`estimated_<col>` and five other spellings). A card cannot yet name its output column; a model whose output does not follow the names is refused.
- Units come from the column name's suffix (`_um` → uM), not from the card.
- Score is offered in discovery portfolios only; the function works anywhere a card names what it predicts.
- Predictions are not in the program-data conditions: a gate is never met by a prediction.

## Approval workflow

*F-30, FR-57. Schema: 0.2.23: table `approval_request`, `portfolio.approval`.*

- `approval_request` is table-secured, View and Edit for all users, like `program_role`: the app checks who proposes and reviews, the server does not. Master security through the program would hide requests from approvers who are not the program's author until GROK-21041 is fixed.
- A second approver who is not the program's author cannot review end to end today: reading the program row needs GROK-21041. The package test files the request as if a throwaway user proposed it and approves as admin; it also checks the throwaway user is refused as a non-approver. The BDD harness has one login, so its scenario checks that the proposer gets Withdraw, not Approve.
- Only approvers propose gate decisions. Contributors submit documents and cards but cannot propose a decision.
- The gate log's `decided_by` is the approver who made it take effect; the proposer is on the request (and shown in the history). No column changed.
- The evidence is frozen at the proposal; the readiness is checked again at approval. Nothing re-freezes it.
- Approvers are not notified of a new proposal (*Waiting for me* shows it); only the proposer hears the outcome. A non-admin reviewer may fail to notify someone else (the endpoint lets a user notify only themselves); the review stands and says so.
- `hold` and `recycle` may be listed; EX-BIOD lists only advance and terminate. A held program needs no second person there.
- A rejected document can be submitted again or linked again as a draft. A rejected model card can be submitted again.

## Park orphan folders

*FR-22 (action half). Schema: 0.2.24: table `parked_folder`.*

- Parking applies to whatever the report lists; there is no portfolio setting to turn it on (FR-22 says off by default). Off here means it never runs on its own: no sync or build parks anything.
- One `_Parked/<date>` per day, by UTC date. Two parks on one day share it.
- A parked program folder moves with its files, so a document linked by a `Space:...:Files/...` path inside it may not open until it is unparked. Only ended programs' folders are parked.
- Unpark refuses when a folder of the same name is back at the root; it does not merge.
- Removing an empty `_Parked` on unpark is the one deletion. A file or connection added to it by hand keeps it.
- `parked_folder` is table-secured, View and Edit for all users; the app checks the approver.
- On this stand EX-SMD's report lists EXS-ION4 (terminated). Nothing was left parked: the tests and the scenario park throwaway folders only.

## Evidence tables generated into documents

*F-40, FR-72. Schema: 0.2.25: `document.generated`.*

- The generated document is linked to the criterion like any other, so approving it meets the criterion (Process Lock's process description, for example). An evidence table could instead stay off the criterion and only be cited.
- The kind is `evidence_table` when the portfolio declares it, else `other`. No example portfolio declares it; approving a generated `other` document supersedes any other approved `other` document of the same criterion.
- Only the newest snapshot is read, and only 25 of its rows go into the table; the CSV and Markdown say how many more.
- Data conditions come from the criterion's own gate only. A gate without data conditions says so.
- A decision counts when its evidence cites the criterion, the run selected now, or a document a person wrote. Once another run is selected, decisions that cited only the old run drop out.
- The hashes are a 53-bit string hash (cyrb53), not SHA-256: they tell content apart, they are not a seal. A run and a snapshot are hashed by id and selection time or SHA-256, not by their values, because both are written once.
- Superseding older generated drafts skips the "only approvers mark superseded" rule: the generator replaces its own drafts. A submitted version is left to its review.
- Needs the program's workspace folder; without one, nothing is written.
- Provenance does not yet show a generated document's sources as edges.

## Evidence changed since written

*F-41, FR-73. Schema: 0.2.26: `document.evidence_check`.*

- Staleness is computed live on the program page (one read of selected runs, snapshots, evidence and the gate's readiness per criterion), and stored when it differs. Opening the page is then a write for users who may edit documents; others leave it to the sync.
- Told once per evidence hash. If the evidence changes back and forth, each new state is told again.
- Followers are told about every stale generated document of the program, not only of a series they follow.
- The sync checks every program of the portfolio that has a generated document, whatever the sync changed. Cheap while few documents are generated.
- "Revoked" is a decision whose evidence no longer cites the criterion, its run or its documents (the gate log is write-once). Selecting another run therefore also drops decisions that cited only the old run, and both show as changes.
- A notification to someone else needs a user who may edit users (as for followers, 9f10713145); otherwise the line says "not notified" and the document is still marked stale.

## Pinned programs on the landing page

*Pinning, the My programs widget. Schema: none.*

- A group's pin gives the group's members a tile too, marked *for <group>*. The program page's Pin/Unpin is about the user's own pin only.
- An approver pins for the portfolio's viewers group only, from the page; any other group goes through `pinProgram(..., group)`. Pinning for a group needs a portfolio approver (or an Administrator), not a group admin.
- **Build its folder** is offered to every approver of a portfolio with a workspace root, not only on EX portfolios.
- New data on a tile is the follow's `news_on` after its `seen_on`: what the sync found since the last visit. It is not computed live, so it lags the sync.
- The widget sits at `meta.order: 2` (after Spotlight) and spans two columns; tiles scroll inside its fixed 300 px.
- The tile opens the program in an already open Project Portal view rather than a second one.

## Model card in the Workspace preview

*Power Pack hook `workspaceModelPanel`. Schema: none.*

- The role name `workspaceModelPanel` and the side-by-side layout (card column 360 px) are Power Pack's to own; the Power Pack team should review the hook before it is released.
- The hook is for models only, not for apps, queries or other entities in the Workspace.
- The card shows in full; a compact Workspace version (approval and portfolios on top, the rest folded) may read better in the narrow column.
- The model's own preview form is narrower beside the card; its input values are cut off on a 1600 px screen.

## Quality pass

*A review of the features above. Schema: none.*

- `program_role`, `approval_request` and `parked_folder` grant Edit to all users. The app checks who may change a role or review a request; the server does not, so a user could make themselves a program approver from the console. Closed later by the grant change (see *Write permissions on portal tables*).
- A model card that declares only a modality is offered in every portfolio with that modality. That is what put another package's model into EX-SMD's suggestions. Keep, or require the portfolio too?
- The feed shows program roles as they stand. `program_role` keeps no history, so an earlier role of the same group is not an entry. A `program_role_change` table would give a full history.
- Approvers are told only when the person proposing may notify others (the server allows it to users who may edit users). For everyone else the request still shows under Waiting for me. A server-side notifier would fix this.
- The approval claim covers the double apply; a failure part-way through `recordDecision` can still leave a gate-log row and an open request.
- Program-role folder sync revokes a grant it did not make for a group that ever held a role.

## Unit-operation folders from a program's flow

*F-05 leaf level. Schema: 0.2.27: `portfolio.flow_folders`.*

- Declared flows live in the portfolio's JSON column, not in `process_flow`/`process_step`. Loading them into those tables (shared with a customer package's own loader) needs an owner per flow; not done.
- The pick is not written to `program_flow`, so other readers of that table (recommendation by unit operation) do not see EX flows.
- One flow per process type from the rules; parallel branches (`branch_group`) are folders like any step.
- Process-stage work units (EX-BIOC's USP/DSP/FORM/FILL) are not folders. They overlap with flows; one or the other should drive the tree.
- A step's unit operation must be in a catalogue of the set being loaded; a flow naming another, already loaded portfolio's unit operation is refused.

## Data package manifest

*F-43, FR-43, PORT-11 (as a view). Schema: none.*

- A view, not an authored table: no one can add a line by hand (an external data source, a lab notebook link). That is the "authored table" half of F-43; it needs a `data_package_item` table and an owner.
- Criteria are those of the stages reached and the next gate; later gates are left out.
- Documents of every status are listed, superseded ones too, so the history shows; a filter could hide them.
- Model inputs are not lines: a run's line references its saved call, which holds them. Provenance already walks a run to its input tables.
- Datasets are listed by connection, not by the query or its row count, so the manifest does not run the program's queries.
- No CSV-injection guard: a cell starting with `=` is written as is.

## Cross-model run comparison

*F-29, FR-55, PORT-48. Schema: none.*

- Seeding needs a Compute2 change: `CompareRuns(runIds: string[])` (or a `funcCalls` input) that adds each run with
  `entryFromFuncCall`. Small and additive; not made.
- A reviewer comparing other people's runs cannot use Compare Runs for them; Compare here (the portal's side by side)
  still works. A Compute2 option to list all authors would lift this.

## Observations

*F-26, FR-76/106, PORT-91. Schema: 0.2.28: table `observation`.*

- One subject per observation. A conclusion across two series is written on the program, or twice.
- No edit, only retract and write again. A typo costs a retraction.
- Evidence is runs, snapshots, links and files. Documents and decisions are not offered as evidence yet; a compound's measurements are cited through a snapshot.
- An observation on a compound is about the compound in this program; the same compound in another program does not see it.
- Not yet a node in Provenance, and not in the data package manifest.
- The evidence json is stored as `{items: [...]}`: a json column holds a map, not a list (the server refuses a list with "Cannot convert ... to map").

## Design sets and selection sets

*F-23, FR-101/102, PORT-85/89. Schema: 0.2.29: tables `design_set`, `design_proposal`, `selection_set`.*

- A design set is written once: no proposals added later, no rescoring. A new round is a new set.
- Ranking is by one score. Multi-parameter scoring (MPO) belongs to the tool that made the set; its single score is what is stored.
- Structures are stored as given (SMILES), not canonicalized; matching a proposal to a registered compound by structure is an exact string match (used by F-24).
- A selection set lists compounds the program already holds. Picking from a screen's full hit list, before the compounds are in the program, needs the registry to hold them first.
- `chooseForSynthesis` is the request, until request objects (F-49) exist.

## DMTA status

*F-24, FR-105, PORT-86. Schema: 0.2.30: `portfolio.dmta`.*

- *Made* has no date: the registry projection does not carry a registration or batch date. A compound feed column for it (or the batch's date) would give a design → make time.
- *Analyzed* is an observation on the compound or a role set by hand. A role the registry sets is not counted, and neither is an observation on the series.
- A registered compound never designed here starts at *made*; nothing is inferred backwards.
- *Tested* is any dated row of the compound in the chosen dataset, whatever the assay.
- Read on every open of a discovery program page (design sets, compounds, the assay dataset, observations, role changes, plug-ins). Fine on EX data; a large program may want it behind More or cached.
- Structure matching is exact SMILES. Canonicalizing needs RDKit in the browser (Chem package).

## Program-page panels per persona

*F-31, FR-34, PORT-09. Schema: none (`portfolio.personas`).*

- Unlisted panels go behind More rather than away: nothing is out of reach. A strict mode (hide them) is one flag if wanted.
- The gate can be put behind More for a persona. The Gate quick-action card then scrolls to nothing until More is opened.
- Panels are per persona per portfolio, not per user. A user's own layout (drag to reorder) is not built.
- Composition overrides the mode for that persona: Simple and Advanced look the same for a composed persona.
- Counting what is behind More builds those panels once more; cheap on EX data.

## Write permissions on portal tables

*Security. Schema: 0.2.31: three tables View for all users; metadata only.*

- `follow` is writable by every user, everyone's rows (low impact: follow state and news marks). Closing it needs row security, which the server refuses on a deployed table: a new row-mode table and a copy.
- Grants are per table: an approver of one portfolio can write another portfolio's roles and parked folders from the console; a contributor of any portfolio can write any request's status and reviews (a self-approval that looks approved, though the effect, the gate log and the program's stage, still needs Edit on the program row).
- When GROK-21041 is fixed and contributors get Edit on program rows, every master table opens to them too: a contributor could write the gate log, evidence, document status and criterion status directly, bypassing approval. Needs per-row write rules or a server-side write path.
- Functional, not security: non-admin approvers cannot carry out the effects that write master or reference tables (gate log, program stage, model card approval) until program rows can be shared.
- Program-role contributors and approvers outside the portfolio's groups get no table grant: they cannot propose or review from the app (they cannot see the program today either, GROK-21041).
- The rest needs the platform: a write rule per row for table-mode tables (e.g. Edit only where `portfolio_id` is one of the caller's portfolios), or a supported server-side write identity for a package. Not filed; worth raising alongside GROK-21041.

## Three small fixes

*User names instead of ids; no new folders for ended programs. Schema: none.*

- Ended means `completed` and `terminated`. `transitioned` and `missing_from_feed` still get new folders; both are park candidates already.
- Not configurable per portfolio: there is no portfolio-level settings column to hold it without a schema change.
- An ended program with no folder gets none, and its `space_id` stays empty.
- Names are kept for the session; a renamed user shows the old name until reload (`clearUserNames()` resets).

## Program rows visible through their folders

*Security, GROK-21041 workaround. Schema: none.*

- Verified on the stand with a plain user in one portfolio's viewers group: before, 0 programs on every board and a
  program's address fell back to the empty board; after, the portfolio's programs and their children (gate log,
  documents, evidence, runs), and none of another portfolio's.
- The row reaches whoever may view its folder. Anyone an admin shares a program folder with by hand also sees the
  program row. That is the intended rule (the folder is the program's), but it means folder sharing is now row
  sharing.
- View only. Contributors and approvers who are not the row's author still cannot write the master tables (gate
  decisions, evidence, documents, criterion status) on the server. That needs Edit on the row: GROK-21041.
- A program without a folder (a portfolio with no workspace, or before the builder runs) stays private to its
  author. Linking into the portfolio's root instead was not done: it would show every row in the root.
- `linkProgramRows` is per portfolio and not run on load; existing stands run it once per portfolio, or rebuild
  the portfolio's workspace.
- A plain user also needs the package released (`grok publish --release`): a debug publish is visible to its
  publisher only (platform finding 28), and the app answers "Application not found".

## Empty lists in filters

*Fix. Schema: none.*

- `queryIn` covers the reads that built `id in (...)` from a list (feed, provenance, evidence, run selections).
  Reads that build the list by hand and already skip an empty one (workspace flows, new-data, evidence tables,
  profile, portfolio table batches) were left as they are.
- My feed now says nothing for a follow whose program is gone. Such follows are not turned off; a clean-up of
  follows on deleted programs is not built.

