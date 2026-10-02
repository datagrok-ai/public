# Project Portal Core changelog

## 0.3.0 (unreleased)

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
