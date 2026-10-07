# Schema notes

The `portal` domain schema: `databases/portal/schema.json` (version 0.2.33, 41 tables). The design rules are in
[architecture.md](architecture.md#2-schema-design-rules). This page records how the tables group, who may write
them, the choices made while building, and how the schema grew.

## Tables

| Group | Tables |
|---|---|
| Spine | `portfolio`, `registry_feed`, `stage_model`, `stage`, `program`, `stage_transition` (the gate log) |
| Criteria | `gate_criterion`, `criterion_status`, `criterion_artifact` (which run answers a criterion) |
| Work units | `campaign`, `series`, `study`, `process_stage` |
| Compounds | `compound` (registry identity), `program_compound` (a compound in a program, with its values), `compound_role_change` |
| Models | `model`, `model_card`, `model_run`, `prediction`, `unit_operation`, `run_type` |
| Flows | `process_flow`, `process_step`, `flow_rule`, `program_flow` |
| Knowledge | `evidence`, `program_snapshot`, `document`, `observation` |
| Discovery | `design_set`, `design_proposal`, `selection_set` |
| Workflow | `approval`, `program_role`, `parked_folder`, `following`, `sync_log`; superseded and closed: `approval_request`, `follow`, `sync_run` |

## Security and who writes

The server enforces the table's security mode and its grants; checks in package functions are the app's (see
[platform-findings.md](platform-findings.md#what-the-server-enforces-on-writes)).

| Tables | Security | Who reads, enforced by the server | Who writes, enforced by the server |
|---|---|---|---|
| `approval` | row, linked into the program folder (a model card's into the portfolio root) | the portfolio's groups and the program's roles, through the folder (View only); the proposer | each portfolio's contributors group (Edit on the table; any row) |
| `sync_log` | row, linked into the portfolio root | the portfolio's groups; the run's author | each portfolio's approvers group |
| `following` | row, no row grants | its author only (admins included) | its author |
| `program_role` | table | all users | each portfolio's approvers group |
| `parked_folder` | table | each portfolio's approvers group | each portfolio's approvers group |
| `approval_request`, `sync_run`, `follow`, `compound` | table, closed | admins | admins |
| `stage_transition`, `criterion_status`, `evidence`, `document`, `program_snapshot`, `model_run`, `prediction`, `observation`, `criterion_artifact`, `campaign`, `series`, `study`, `process_stage`, `program_compound`, `compound_role_change`, `program_flow`, `design_set`, `design_proposal`, `selection_set` | master → program | Edit on the program row: its author (the sync) and admins, until program rows can be shared (GROK-21041) |
| `program` | row, private by default, linked into its folder | its folder's viewers (View only) | its author and admins |
| `stage_model`, `stage`, `gate_criterion`, `portfolio`, `unit_operation`, `run_type`, `process_flow`, `process_step`, `flow_rule`, `model`, `model_card` | table, View for all users | all users | admins |
| `registry_feed` | table, no grant | admins (it holds SQL) | admins |

The master tables (the row before `program`) are read by whoever reads the program row. A redeploy never revokes a
grant, so `loadPortfolio` (and `securePortalTables` for every portfolio) grants each portfolio's groups and takes
All users' grants away at runtime: Edit and Delete on the workflow tables, View as well on the closed ones and
`parked_folder`. Both need Share on the tables.

**Moving records (0.2.33).** Table security is one grant for the whole table, so the per-portfolio and per-user
records moved to row-secured tables: `approval` and `sync_log` (linked into a space, the pattern program rows
use) and `following`. A compound's name, structure, aliases and batches moved onto its `program_compound` links,
secured through the program; `compound` keeps the registry and code. `moveToSecuredTables(code, dryRun)` (admins)
copies a portfolio's requests and sync runs (same ids), links them and moves the compound values; compounds are
emptied once every link carries them. Follows are not copied: the server makes the caller a row's author, so a
copy would be the copier's. A redeploy cannot change a deployed table's security mode, hence new tables.

`approval` and `sync_log` use `defaultRowVisibility: table` so that Edit on the table lets the writer group write
any row (a reviewer is not the proposer). The flip side: a View grant on these tables would show every row, so no
group may hold one; `securePortalTables` takes All users' away if it appears.

**Still readable by all users:** `program_role` (group ids per program; everyone reads it to know their own
access), `portfolio` (configuration, including connection names and dataset SQL), and the catalogue (`model`,
`model_card`, `unit_operation`, `run_type`). The app scopes the catalogue to the user's portfolios; the server
does not.

**Remaining gaps.** Grants are per table, so an approver of one portfolio can write another portfolio's roles
from the console, and a contributor can write any request's status. Once program rows can be shared and
contributors get Edit on them, every master table opens to them too. Both need a per-row write rule or a
server-side write path from the platform.

## Choices made while building

- `criterion_status` is the per-program instance of every criterion kind, with owner, due date, document and
  notes. One status set: `not_started`, `in_progress`, `met`, `not_met`, `waived`. Data conditions use only `met`
  and `not_met`.
- `criterion_status.evidence` carries free-text notes rather than a second text column.
- `criterion_artifact` carries `program_id` as well as `criterion_status_id`: master security delegates one hop
  only.
- `gate_criterion.modality` (nullable) scopes a criterion to a modality; null means every modality.
- `flow_rule` is keyed by a `code` (`<modality or ANY>:<flow>`), because a nullable modality cannot serve in a
  business key.
- `process_flow.process_type` is free text, not choices, so the split into process types stays configuration.
- Tables are created in declaration order, so a table is declared after the tables it references.
- Lists that would be `string_list` are `|`-separated strings (`stage.spellings`, `portfolio.workstreams`), because
  a list cannot be written from JS ([platform-findings.md](platform-findings.md), item 2).
- A `json` column holds a map, not a list, so lists are stored under one key (`{datasets: [...]}`,
  `{items: [...]}`).
- The user column of `follow` is `user_id`: a column named `user` breaks the deploy.
- `follow` and `program_role` turn rows off (`active`) instead of deleting them.
- `document.doc_version`: `version` is a system column of every domain table.
- `model_run.run_id` is the call's id as text, always written; `func_call`, the soft reference, is written when the
  server accepts it (GROK-21042).

## Migrations

`databases/portal/0001_model_run_run_id.sql` (with `down/`) adds `model_run.run_id` and makes `func_call`
optional. **Debt:** once GROK-21042 is fixed, backfill `func_call` from `run_id` with SQL (the column is
write-once through the API), then make it required again.

Every other change so far has been additive and deployed from `schema.json` alone.

## Version history

| Version | Change |
|---|---|
| 0.2.0 | The first core: portfolios, stage models, programs, gate log, criteria, model catalogue, runs, flows |
| 0.2.1 | `portfolio.workstreams`, `program.space_id` (workspace per portfolio) |
| 0.2.2 to 0.2.5 | Data conditions, readiness, work units, stage folders, compounds |
| 0.2.6 | `model_run.run_id`; migration 0001 |
| 0.2.7 | Model catalogue and cards |
| 0.2.8 | `gate_criterion.artifact`, `run_filter` |
| 0.2.9 | `portfolio.program_datasets` |
| 0.2.10 | Attention rules |
| 0.2.11 | `portfolio.schedule`, `portfolio.owners`, `sync_run` |
| 0.2.12 | `portfolio.auto_advance`, `stage_transition.source`, `stage_transition.evidence` |
| 0.2.13 | `portfolio.handoff`, `unit_operation.portfolios`, `run_type.portfolios` |
| 0.2.14 | `model_card.applies_to` |
| 0.2.15 | `follow` |
| 0.2.16 | `portfolio.personas` |
| 0.2.17 | `program_role` |
| 0.2.18 | `evidence` |
| 0.2.19 | `program_snapshot` |
| 0.2.20 | `document`, `portfolio.document_kinds` |
| 0.2.21 | `follow.news_on` |
| 0.2.22 | `prediction` |
| 0.2.23 | `approval_request`, `portfolio.approval`; new choices on `document.status` and `model_card.approval` |
| 0.2.24 | `parked_folder` |
| 0.2.25 | `document.generated` |
| 0.2.26 | `document.evidence_check` |
| 0.2.27 | `portfolio.flow_folders` |
| 0.2.28 | `observation` |
| 0.2.29 | `design_set`, `design_proposal`, `selection_set` |
| 0.2.30 | `portfolio.dmta` |
| 0.2.31 | `program_role`, `parked_folder`, `approval_request`: View for all users (was View and Edit); metadata only |
| 0.2.32 | `stage_transition.decision` takes `record`; `stage_transition.outcome` (go, hold, no_go) for decision records; additive |
| 0.2.33 | `approval`, `sync_log`, `following` (row-secured replacements of `approval_request`, `sync_run`, `follow`); `program_compound.compound_code`, `compound_name`, `smiles`, `aliases`, `batches`, `details`; additive |
