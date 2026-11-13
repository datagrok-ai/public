# CLAUDE.md

## Overview

**NCA** is a Datagrok package for Non-Compartmental Analysis of preclinical pharmacokinetic data. It takes a
concentration-time table to a validated set of PK parameters (Cmax, AUC, λz, t½, CL, Vz, Vss, MRT, %F, TK
ratios, dose proportionality) and lets the scientist review, stratify and save the result in studies.

The app is in **Apps > Compute > Non-Compartmental Analysis**.

## Product principles

1. **The math lives in `@datagrok-libraries/sci-comp`, not here.** `nca.computeNca`, λz best-fit, BLQ
   handling, AUC methods and derived parameters are validated against PKNCA in sci-comp. This package calls
   them on DataFrames and presents the results. Never re-derive or vendor NCA math in this package; if a
   parameter or primitive is missing, add it to `libraries/sci-comp` first.
2. **The audience is the DMPK scientist**: discovery screens (rank compounds by CL, t½, Vss, %F) and
   IND-enabling / TK studies (serial and sparse sampling, defensible λz).
3. **Interactivity is the value**: click-to-include λz points with live refit, linked selection between
   profiles and parameter tables, problem-first ordering of the results.
4. **Show honest uncertainty.** Small N, fragile λz fits (low adj. R², few points, short span) and high
   %AUCextrap are flagged, never hidden behind a single number. Warnings from
   `ComputeResult.provenance.warnings` reach the user.

## Architecture

```
Datagrok layer: views, viewers, dialogs, events       src/views, src/tour
Controller: recompute, λz overrides, review state     src/controller
Adapter: DataFrame <-> sci-comp, validation, storage  src/adapter
Headless pipeline functions                           src/pipeline
Computation core                                      @datagrok-libraries/sci-comp (nca)
```

- **adapter/**: column contract and auto-detection, validators, unit conversion, profile index and inputs
  (Float64Array at the core boundary), results writer, aggregation, F%, TK, dose proportionality, sparse
  composite profiles, acceptance verdicts, study grid, storage (`ems-nca-repository.ts`).
- **controller/**: `NcaController`, which owns recompute, manual λz point overrides and review decisions.
- **pipeline/**: `ncaDetectMapping`, `ncaCompute`, `ncaSummarize`, `ncaSaveToStudy`, `ncaSeedDemoStudies`:
  headless functions for Flow and scripts. They must not import `views/` (enforced by `.eslintrc.json`).
- **views/**: Hub, analysis view (mapping pane, results grid, faceted concentration-time plot, λz editor,
  review pane), study views, Browse tree.

## Storage

Studies and analyses are stored in the platform domain schema `nca`, defined in
`databases/nca/schema.json`. `grok api` (run by `grok build`) generates the typed client
`src/generated/db.ts` from it. Never edit the generated file; change the schema and rebuild. Table snapshots
are kept in a `file` column, and every write of an analysis and its child rows is one versioned transaction.

## Build and test

```bash
grok build --typecheck                      # bundle + type-check (from this folder)
pnpm run lint
grok test --host local                      # all tests; requires a running Datagrok server
grok test --host local --category "Adapter: column validators"
grok publish local                          # deploy a debug build
```

The reference suite (`src/tests/reference-*-tests.ts`) checks the parameters against PKNCA results for
theoph, indometh, IV infusion, the rat routes, the multi-compound screen, TK and the NONMEM-format datasets in
`files/`; `sparse-nca-tests.ts` covers the sparse composite profiles. A parameter that is not
covered by a reference assertion is not considered validated. `files/REGEN.md` describes where each fixture
comes from.

Runtime dependency: PowerGrid (the `Tags` renderer draws the λz acceptance-reason pills).

## Conventions

- A grid gets **one** `onCellPrepare` and **one** `onCellTooltip` handler. A second registration replaces
  the first and silently removes its tooltips or styling.
- Verdict colours (Accepted / Flagged / Not estimated) come only from `views/status-palette.ts`; the
  results grid, facet tiles and editor header must read the same constants.
- The colour-banded results columns are derived from `buildColourRules`; do not restate the set elsewhere.
- Percent-derived parameter labels use the `%` prefix (`%AUCextrap`, `%F`, `%CV`); code identifiers keep
  word prefixes (`pctExtrap`, `F_PCT`).
- UI strings are English. Use platform viewers and inputs; no custom canvas or SVG renderers.
