# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

TypeScript Datagrok app for exploring and validating CDISC [**SEND**](https://www.cdisc.org/standards/foundational/send)
(Standard for Exchange of Nonclinical Data — preclinical / animal study data). Registered as the **Preclinical Case**
app (`#app` in `package.ts`); surfaces a browse-tree of studies where each study exposes per-domain
table views plus cross-domain analysis views (Validation, Matrix, Measurements, Microscopic
Findings, Observation Timelines) and a Study Summary. Loads `.xpt` files (SAS transport) via the
`readSas.py` script, caches parsed DataFrames to `.d42`, and validates through a Celery /
CDISC-CORE Docker container.

## Architecture

Start here — the rest can be read off these:

- **`src/package.ts`** — app entry, `#app` / `#fileHandler` / `#appTreeBrowser` registrations,
  and the two dispatch maps: `SUPPORTED_VIEWS` (names) and `VIEW_CREATE_FUNC` (studyId → view).
  Adding a view means editing both.
- **`src/preclinical-study.ts`** — `PreclinicalStudy` is the in-memory model of a loaded study.
  Owns `domains: PreclinicalDomains` (one nullable `DG.DataFrame` per SEND domain code, plus
  `supp: DG.DataFrame[]`), lifecycle flags (`initCompleted`, `validated`, `loadingStudyData`),
  and the `validationCompleted` RxJS `Subject`. `init()` → `process()` (DM joins, baseline calc) →
  `validate()` (async CDISC-CORE call). `ensureValidationColumnsAdded()` waits for validation,
  then decorates each domain DataFrame with error columns.
- **`src/utils/app-utils.ts`** — the study-loading pipeline and the `studies` registry
  (`{[studyId]: PreclinicalStudy}`, module-level). `createStudiesFromAppData` walks `SEND/`,
  `createStudyWithConfig` builds one study, `initStudyData` / `readClinicalData` do the lazy
  file-to-DataFrame load (with `.d42` fast path), `openStudy` / `loadView` drive navigation.
- **`src/data-preparation/data-preparation.ts`** — `calculateLBBaselineColumns` (per-subject/test
  baseline, change, pct-change, max/min post-baseline on the `lb` domain) and
  `createAllMeasurementsDf` (builds the unified `Measurements` DataFrame by cloning required
  `--TEST / --STRESN / --DY` columns across all eligible domains).
- **`src/utils/views-creation-utils.ts`** — `createTableView(studyId, name, factory, helpUrl)`
  wraps a `StudyTableViewParams` factory into a `DG.TableView` and hides validation columns.
  Every table-style view uses this.
- **`src/types/types.ts`**, **`src/types/validation-result.ts`**, **`src/constants/*.ts`** —
  shared types and string constants (SEND column names, view names, domain categories). Import
  column names from `columns-constants.ts` rather than writing string literals.
- **`dockerfiles/`** — Python Celery worker that downloads and runs the
  [cdisc-rules-engine](https://github.com/cdisc-org/cdisc-rules-engine) binary (fetched on first
  run from GitHub releases, cached inside the container) against SENDIG 3.1. Called from TS via
  `grok.functions.call('PreclinicalCase:run_core_validate', {...})`.

### Views (two patterns)

- **ViewBase subclasses** (`StudySummaryView`, `ConfigurationView`) — have a `load()` method for
  lazy init, own their DOM. Set `helpUrl` in constructor.
- **Table views** (Validation, Matrix, Measurements, Microscopic Findings, Observation Timelines)
  — factory functions returning `StudyTableViewParams`, wrapped by `createTableView()` in
  `views-creation-utils.ts`. `helpUrl` is passed through `createTableView`.

Each view has a corresponding help file in `views_help/` wired via `helpUrl`.

Adding a new view:

1. Create the view function/class in `src/views/`.
2. Add the view name constant to `constants/view-names-constants.ts`.
3. Add it to `SUPPORTED_VIEWS` and `VIEW_CREATE_FUNC` in `package.ts`.
4. Add a `views_help/<name>.md` help file and wire `helpUrl`.
5. The tree browser picks it up automatically from `SUPPORTED_VIEWS`.

### Validation lifecycle

`PreclinicalStudy.validate()` runs in the background via the Docker container. Views must handle
the case where validation hasn't completed yet — check `study.validated` and subscribe to
`study.validationCompleted`.

## Glossary — SEND concept → code

| Concept                        | Code                                                | Notes                                                                                     |
|--------------------------------|-----------------------------------------------------|-------------------------------------------------------------------------------------------|
| Domain (`dm`, `lb`, `mi`, ...) | field on `PreclinicalDomains`                       | CDISC-standardized codes; full list in `constants/domains-constants.ts`.                  |
| Supplemental domain (`supp*`)  | `PreclinicalDomains.supp`                           | Always a list; matched by `df.name.startsWith('supp')`.                                   |
| Subject                        | column `USUBJID` (`SUBJECT_ID`)                     | Primary join key between domains and DM.                                                  |
| Study day                      | `--DY` column (`LBDY`, `BWDY`, `MIDY`, ...)         | Per-domain integer; `VISITDY` is the cross-domain fallback used by `createVisitDayStrCol`. |
| Test name / result             | `--TEST` / `--STRESN` (e.g. `LBTEST` / `LBSTRESN`)  | Test name + numeric standardized result; same pattern across domains.                     |
| Validation issue               | `ValidationResult.Issue_Details[i]` (`IssueDetail`) | Row is 1-indexed in the JSON — `addValidationColumnsToDomain` subtracts 1.                |
| Per-column error columns       | `<VAR>_hasErrors`, `<VAR>_errors` (JSON string)     | Added by `addValidationColumnsToDomain`; hidden by `hideValidationColumns`.               |

## Conventions specific to this package

- **`.d42` cache can go stale.** After first load, all domain DataFrames are serialized to a single
  `.d42` binary file via `grok.dapi.files.writeBinaryDataFrames()`. Subsequent loads read from this
  cache instead of re-parsing individual `.xpt` files. If you change data-processing logic,
  previously-cached `.d42` files are stale and must be regenerated.
- **DM join tagging — don't remove.** When domains are joined with DM, columns originating from
  DM are tagged with `COLUMN_FROM_DM_TAG`. Validation checks this tag to avoid flagging DM-sourced
  columns against other domains' rules.
- **BG domain is excluded from unified measurements.** Body weight gain (BG) is excluded from
  `createAllMeasurementsDf` because it uses period-based semantics (BGENDY) rather than
  point-in-time study days. Handle it separately if you need it.
- **Baseline: `BLFL='Y'` → else earliest study day.** Two-tier logic per SENDIG v3.1.1 — `BLFL`
  flag is authoritative where present, earliest study day is the fallback when it's absent. Don't
  simplify to a single rule.
- **File storage path.** Study data lives at `System:AppData/Preclinicalcase/SEND/<StudyName>/`.
  The lowercase `Preclinicalcase` (not `PreclinicalCase`) matches `package.json` — do not
  "correct" it. Per-study cached files: `<StudyName>.d42` (DataFrames), `validation_results.json`
  (CDISC CORE output), `study_config.json` (parsed study config).
- **SEND column naming is strict.** Use the standardized names (`USUBJID`, `LBTEST`, `LBSTRESN`,
  `--DY`, etc.). Import from `columns-constants.ts` rather than writing string literals.
