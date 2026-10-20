# Architecture

The reference architecture Project Portal is built towards: one core package with one domain schema, and
customer configuration as data. The last section says where the build differs from the target.

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│  PRESENTATION  (core, configuration-composed; u2 over `portal` tables)            │
│  Portfolio board · Portfolio table · Program page · Gate panel                    │
│  Simple/Advanced · Quick actions and panels per persona · Model cards · Run save/ │
│  select/compare · Attention list · DMTA status · Approvals · Provenance · Feed    │
├──────────────────────────────────────────────────────────────────────────────────┤
│  WORKSPACE  (Spaces, generated)                                                  │
│  Root → program → {workstream | stage folder | process type → unit operation}     │
│  Rows and runs LINKED into spaces; compounds, results, predictions never folders │
├──────────────────────────────────────────────────────────────────────────────────┤
│  DOMAIN CORE  (one schema, `portal`)                                             │
│  spine      portfolio · program · compound · program_compound(role) ·            │
│             stage_model · stage                                                  │
│  work units campaign · series · study · process_stage (a portfolio enables some) │
│  criteria   gate_criterion(kind: document|data_condition|decision) ·             │
│             criterion_status · criterion_artifact(selection) · program_snapshot  │
│  models     model · model_card · model_run · prediction                          │
│  design     design_set · design_proposal · selection_set                         │
│  knowledge  document · evidence · observation · approval_request                 │
│  vocab      unit_operation · run_type · process_flow · process_step · flow_rule  │
├──────────────────────────────────────────────────────────────────────────────────┤
│  PROJECTIONS OVER PLATFORM AND PACKAGE ENTITIES                                  │
│  model_run → FuncCall · document → file · program data → query on a connection · │
│  compound → registry projection · evidence.target → soft ref (+ frozen name)     │
├──────────────────────────────────────────────────────────────────────────────────┤
│  FEEDS  (per portfolio)                                                          │
│  program registry · work units · compound registry · program datasets            │
├──────────────────────────────────────────────────────────────────────────────────┤
│  PLATFORM                                                                        │
│  Domain schemas · Spaces · FuncCall · row security · favorites · notifications · │
│  user settings · u2 · Power Pack (Home, Workspace) · Compute2                    │
└──────────────────────────────────────────────────────────────────────────────────┘
```

## 1. Packages

| Package | Contents |
|---|---|
| `ProjectPortal` (this package) | The `portal` schema; sync, compound binding, workspace builder, readiness, attention, recommendation, scoring, snapshots, documents, approvals, provenance; the app and its views; the portfolio loader and validator; three synthetic example portfolios |
| A customer package | Portfolio files (stage models, criteria, feeds, groups, catalogue, personas) handed to the core with `loadPortfolioConfigs`, plus any loaders or functions only that customer needs. Kept outside this repository |
| Packages the portal calls | Compute2 (run comparison), Power Pack (Home widgets, Workspace preview); discovery tools the portal can project over or open with the program's context |

The rule of the split: **the core has no customer nouns and no setting-specific noun in code**, and a portfolio
file has no code. If a second customer asks for the same thing, it moves into the core.

## 2. Schema design rules

1. **One schema.** A `ref` cannot cross plugins, so customer packages add rows and use `extensible` tables; they
   never declare a second schema.
2. **Security modes.** `program` is `row`, private by default. Work units, runs, criteria state, documents,
   evidence and the other per-program records are `master` → program. Vocabulary and configuration tables are
   `table`.
3. **Never promote** compounds, results, wells or predictions to entities.
4. **Business keys everywhere;** composite keys are handled explicitly, so every load and sync is idempotent.
5. **Reference tables over `choices`** where definitions are content or rows need external codes.
6. **Frozen beside live:** the stage at run, a role at decision, a name frozen beside every soft reference, a
   snapshot as frozen data.
7. **Migrations from day one** (`databases/portal/0001_*.sql`, with `down/`).
8. **The schema is named for the product** (`portal`), never for a customer or an upstream database.
9. **Setting-neutral columns.** No column in the spine or criteria tables presupposes a process or a discovery
   noun; those live in work-unit and vocabulary tables a portfolio enables.
10. **Data conditions are one rule language,** shared by gate criteria, attention rules and missing assays: a
    declarative predicate over the program dataset (column, operator, value, scope = lead set | series | program,
    minimum count).

## 3. Workspace builder contract

- A plan of nodes first, with a preview (`dryRun`) that lists every folder it would create.
- Additive by default: create what is missing; never rename, move or delete. Pruning is out of scope; parking is a
  separate, explicit step.
- One build per root at a time.
- A program's folder id is recorded on the program, so a rename in the registry keeps the folder.
- **A level may only be generated from `program`, a work unit, a declared workstream list, a stage's folders, or
  a process flow's steps.** Generating from `compound`, `plate`, `result` or `prediction` is refused.

## 4. Discovery assets: what goes where

| Concern | Where | Why |
|---|---|---|
| Structure, registration code, batches | The customer's registry, projected as `compound` | Never mint compounds |
| Structure rendering | The platform's Molecule semantic type | Platform, not portal |
| SAR table | The program dataset: a governed query opened as a table view | Live, not stored |
| Predictions | `prediction` rows from batch scoring | Predicted vs measured needs rows |
| Design set | `design_set` with structures as a Molecule column; chosen rows stand in for requests | Left half of DMTA |
| Requests | An external system, read through a pluggable DMTA source | Separate product |

## 5. Runs and predictions: what goes where

| Concern | Where |
|---|---|
| Execution record | FuncCall |
| Scientific meaning of a run | `model_run` |
| Per-compound outputs of a scoring run | `prediction` |
| Which artifact answers a criterion | `criterion_artifact` |
| Appearance in a program's space | A link into the stage folder (`spaces.addEntity`) |
| Access to the raw run | Follows the function's visibility (the call has no ACL); row security on `model_run` and sharing of the saved tables give the program scope |

## 6. Persona modes and pages

Simple and Advanced are app modes, not a separate shell. Pages are u2 components over `portal` tables. Quick
actions are cards declared per persona in the portfolio file; each runs an action the page already has, or a
function that receives the program's context (portfolio, program, stage). A portfolio can compose the program page
per persona from named panels.

## 7. As built

The build follows the target with these differences, each recorded in [build-log.md](build-log.md):

- **No `model_applicability` table.** Applicability is `model_card.applies_to`, per model version.
- **No `decision` table.** The gate-log entry (`stage_transition`) is the decision; `evidence` rows hang off it.
- **`prediction` is master-secured through the program,** not through the run.
- **No modality or other vocabulary tables yet;** `program.modality` is a string.
- **Approval is optional per portfolio** (`approval`); selection is the default.
- **Process flows** can be declared in the portfolio file (`flowFolders`) or loaded by a customer package into
  `process_flow`, `process_step` and `program_flow`.
- **Nothing runs on a schedule by itself:** a package cannot register a scheduled function
  ([platform-findings.md](platform-findings.md), item 9).
