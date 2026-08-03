# CLAUDE.md — Biologics package

Demo-only Datagrok package that deploys a self-contained PostgreSQL schema for biologics
data (antibody chains, peptides, ADCs, linkers, drugs, assays, dose-response curves) and
wires the schema into Datagrok's DB Explorer so identifiers like `GROKMOL-000001` become
clickable entity references everywhere in the UI.

**Not for production.** No real samples or assays — every row is generated from bundled
files at `createDemoBiologicsData` time.

## Layout

```
packages/Biologics/
  databases/biologics/        # PostgreSQL schema, migrations applied in lexicographic order
    0000_init.sql             # sequences, drugs, linkers, adc, assay_types, target_organisms,
                              # assay_results, purification_batches, expression_batches
    0001_init.sql             # adds adc_id FK to assay_results
    0002_peptides.sql         # peptides table (HELM) + peptide_id FK
    0003_regions_targets.sql  # sequence_regions, targets, sequence_properties, sequence_liabilities,
                              # expanded assay_types panel (fit_model, category), sequence_id FK
    0004_seq_chains.sql       # splits sequences.sequence into heavy_chain + light_chain;
                              # adds chain column to regions/properties/liabilities/expression_batches/
                              # purification_batches; adds assay_curves
  queries/
    demo.sql                  # 6 named queries (see below)
    antibody_profile.sql      # comprehensive antibody profile by organism + target
    antibody_profile.layout   # ViewLayout with filters/grid/form/pivot/pie chart
  files/
    antibodies.csv            # 493 rows: Antibody_ID, AntibodyHC, AntibodyLC, Antigen_ID, Antigen, Y
    random_curves.csv         # 36 4PL dose-response curves as JSON (single 'curve' column)
  src/
    package.ts                # entry point: info, createDemoBiologicsData, populateAdcGlyphs, autostartBiologics
    config.ts                 # biologicsConfig — drives all DB Explorer behavior
    helms.ts                  # 91 head-to-tail cyclized HELM V2.0 strings used to seed peptides
    randsmiles.ts             # 1000 small-molecule SMILES used to seed drugs
    glyphs/glyphs.ts          # auto-generated; exports glyphPool (3 base64 PNGs)
    glyphs/glyphConverter.js  # one-shot Node script: re-runs to regenerate glyphs.ts when PNGs change
    package-test.ts           # standard test harness (no actual tests in this package)
  webpack.config.js           # standard datagrok externals + exceljs, html2canvas
```

There is no `src/tests/` folder. PNGs never reach webpack — they are pre-inlined as base64
in `glyphs.ts` at codegen time.

## Schema overview (post-0004)

| Table | Identifier prefix | Notes |
|---|---|---|
| `sequences` | `GROKSEQ-` | `heavy_chain`, `light_chain` (the old single `sequence` column was dropped in 0004) |
| `drugs` | `GROKMOL-` | small-molecule SMILES |
| `linkers` | `GROKLINKER-` | type SMALL or PROTEIN; one of smiles/sequence is non-null |
| `adc` | `GROKADC-` | (antibody_id, linker_id, drug_id) + glyph (base64 PNG) |
| `peptides` | `GROKPEP-` | HELM notation |
| `targets` | `GROKTGT-` | molecular targets (HER2, PD-1, …) |
| `target_organisms` | `GROKORG-` | E. coli, CHO, … |
| `assay_types` | — | `category` + `fit_model` drive curve eligibility |
| `assay_results` | — | FKs to assay_types, target_organism, target, and one of {adc_id, peptide_id, sequence_id} |
| `assay_curves` | `GROKCRV-` | curve JSON (4PL dose-response), FK to assay_results, ON DELETE CASCADE |
| `sequence_regions` | — | per-chain (`chain` ∈ HC/LC); positions are within that chain |
| `sequence_properties` | — | one row per chain (HC and LC each get MW, pI, charge, etc.) |
| `sequence_liabilities` | — | per-chain; `position` is chain-relative |
| `purification_batches` | `GROKPUR-` | per-chain (HC and LC purified separately) |
| `expression_batches` | `GROKEXP-` | per-chain |

Identifier columns are `GENERATED ALWAYS AS … STORED UNIQUE`, so never insert them — let
PostgreSQL compute them from the SERIAL `id`.

## Demo data flow (`createDemoBiologicsData`)

Order matters; FKs cascade. The function:

1. Wipes existing data in dependency order (assay_curves → assay_results → sequence_*
   → adc → expression/purification batches → linkers/drugs/sequences/peptides). Explicit
   delete on `assay_curves` is redundant with ON DELETE CASCADE but kept for clarity.
2. Loads [files/antibodies.csv](files/antibodies.csv) and [files/random_curves.csv](files/random_curves.csv) via
   `_package.files.readCsv(...)` (returns a `DG.DataFrame`).
3. Inserts `helms` (91 strings) into `peptides`.
4. Picks up to 400 antibody rows; inserts `(sequence_type='PROTEIN', heavy_chain, light_chain, name)`.
5. Generates per-chain regions (VH/CH1/Hinge/CH2/CH3 + CDR-H1..3 with `chain='HC'`,
   VL/CL + CDR-L1..3 with `chain='LC'`), positions clamped to actual chain length.
6. Generates per-chain properties (one row per chain) and liabilities (scans both chains).
7. Inserts 200 drugs from `randsmiles.smi` and 10 linkers (5 PROTEIN, 5 SMALL).
8. Creates ADCs from random (antibody, drug, linker) triples.
9. Creates per-sequence (legacy), per-peptide, and comprehensive (every ADC × every new-panel
   assay type) `assay_results`. The shared `insertAssayResultsBatch` helper uses
   `RETURNING id, assay_id` and queues a curve from the pool for any result whose assay
   type has a non-null `fit_model`.
10. Bulk-inserts queued curves into `assay_curves` (chunk size 50 — curve JSON is large).
11. Backfills `assay_results.sequence_id` from `adc.antibody_id` for ADC-linked rows.
12. Calls `populateAdcGlyphs()` to assign random base64 PNGs to ADCs.

When extending demo data: stick to `RETURNING id` for any FK source, escape strings with the
local `escape()` (single-quote → `''` and strips `@` to avoid `@param` interpretation), and
batch inserts (chunk size 200 typical, 50 for curve JSON).

## DB Explorer integration (`biologicsConfig` in src/config.ts)

`autostartBiologics` calls `DBExplorer.initFromConfig(biologicsConfig)`. This single config
object drives a lot of UI behavior — edit it carefully:

- **`entryPoints`** — each key is a Datagrok semantic type (e.g. `DG_BIOLOGICS_DRUG_ID`).
  Each entry registers a `SemValueObjectHandler` (powers context panel / tooltip / card view)
  AND a `DG.SemanticValue.registerRegExpDetector` so any matching string anywhere in the UI
  (grid cells, scripts, search) auto-resolves to that semantic type and becomes a clickable
  reference. **Removing an entry** un-registers both — identifiers stop being clickable.
- **`joinOptions`** — each entry adds a `LEFT JOIN <tableName> ON <columnName> = <onColumn>`
  with the listed `select` aliases to every query whose result has that FK column. Entity
  cards depend on these; e.g. an ADC card pulls `compound_structure` from `drugs` and
  `antibody_heavy_chain`/`antibody_light_chain` from `sequences` via this mechanism.
- **`headerNames`** — `table → column`. When a card has multiple rows from a related table
  (e.g. several `sequence_liabilities` for one sequence), the accordion uses this column's
  value as the pane label instead of "Row 1 / Row 2". Currently set for `linkers`,
  `sequence_liabilities`, `sequence_regions`, `sequence_properties`,
  `expression_batches`, `purification_batches`.
- **`uniqueColumns`** — `table → column`. Cells matching this column render as a blue
  clickable link that opens the entity's full context panel. Always `'identifier'` here.
- **`customSelectedColumns`** — `table → string[]`. Restricts and orders columns shown in
  that table's card. Columns not listed are hidden even if the query returned them.
- **`addCustomRenderer(predicate, renderer)`** (called after `initFromConfig`) — predicate
  takes `(tableName, colName, value)` and short-circuits the default rendering when it
  returns true. First match wins. `autostartBiologics` registers four:
  SMILES → 2D molecule, HELM (lowercase starts with `peptide`) → HELM widget,
  glyph/image/png with base64 prefix `iVBORw0KGgo` → canvas image. The Curves package's
  own `fit` semantic-type detector picks up `curve` JSON columns automatically.

When adding a new table to the schema, decide whether it deserves an `entryPoint` (only if
it has a stable user-facing identifier prefix like `GROKxxx-######`), `uniqueColumns`,
`customSelectedColumns`, and any join options to enrich its card with denormalized fields.

## Queries

Named queries become functions discoverable across the platform; the `meta.searchPattern`
makes them findable from natural-language search.

- `assaysByOrganism(organism)` — assay results joined with ADC, linker, drug, sequence,
  including the per-result curve.
- `ADCsWithCapsazeActivityHigherThan(minActivity)` — caspase-activity filter.
- `ADCsWithIC50HLThan(value, valueTarget)` — IC50 directional filter (`higher`/`lower`).
- `adcsLinkedToDrug(drugID)` — reverse lookup from drug identifier.
- `getBiologicsPeptideHelmByIdentifier(peptideIdentifier)` — returns HELM string for
  `GROKPEP-######`. Output `semType: Macromolecule, units: helm`.
- `getBiologicsCurveByAssayResultIdentifier(curveIdentifier)` — returns curve JSON for
  `GROKCRV-######`.
- `antibodyProfileByOrganismAndTarget(organism, target)` — comprehensive view paired with
  [antibody_profile.layout](queries/antibody_profile.layout); the layout configures filters,
  pivot table, grid, and form views, and tags `heavy_chain`/`light_chain` as Macromolecule
  fasta sequences and `curve` with `cell.renderer: fit`.

If you change a query's projected columns, update the layout's column metadata — orphaned
formulas (e.g. `Length(${sequence})` pre-0004) silently break.

## Build & deploy

Standard Datagrok build (run from this package folder):

```
npm run build              # grok api && grok check && webpack
npm run debug-biologics    # webpack && grok publish
npm run release-biologics  # webpack && grok publish --release
```

Migrations are applied automatically by the platform when the package is published — they
run in lexicographic order, forward-only. Never edit a published migration; add a new one.

To regenerate `src/glyphs/glyphs.ts` after changing PNGs:

```
node src/glyphs/glyphConverter.js
```

## Common gotchas

- `sequences.sequence` no longer exists (dropped in 0004). Use `heavy_chain` / `light_chain`.
  Anything still referencing the old column will silently return null.
- The `escape()` helper in `package.ts` strips `@` because Datagrok's query layer treats
  `@name` as a parameter placeholder. Don't bypass it.
- Identifier columns are GENERATED — passing a value into them on INSERT raises a Postgres
  error.
- `assay_results` has three mutually-meaningful FKs (`adc_id`, `peptide_id`, `sequence_id`).
  ADC-linked results also get `sequence_id` backfilled from the ADC's antibody. Peptide
  results have `peptide_id` only.
- The semType `DG_BIOLOGICS_*` strings are part of the public contract — other packages or
  saved layouts may reference them. Renaming a semType orphans them.
- `*.g.ts` and `*-api.g.ts` are auto-generated by `grok api` — never edit by hand.
