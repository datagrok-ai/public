# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**HitTriage** (`@datagrok/hit-triage`) is a Datagrok plugin for chemists to assess hit compound quality and manage
molecule campaigns. It bundles **three applications** in one package:

- **Hit Triage** — Upload a molecule dataset, filter/compute properties, submit results to a chosen function, share campaigns.
- **Hit Design** — Sketch molecules in a spreadsheet, calculate properties, organize into stages, share campaigns.
- **PeptiHit** — Like Hit Design but for peptides (HELM notation). Converts sequences to atomic level via the Bio package.

Category: **Cheminformatics**. Browse paths: `Chem` (Hit Triage, Hit Design), `Peptides` (PeptiHit).

## Build Commands

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run lint               # eslint src --ext .ts
npm run lint-fix
npm run link-all           # npm link datagrok-api @datagrok-libraries/compute-utils @datagrok-libraries/utils
npm run build-all          # Builds js-api → compute-utils → utils → this package
```

## Architecture

### Entry Point — `src/package.ts`

- Defines `HTPackage` extending `DG.Package` with an LRU SMILES cache and campaigns cache.
- `PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:
  - Three `@app` entries: Hit Triage, Hit Design, PeptiHit
  - Three `@treeBrowser` entries for Browse panel integration
  - Demo data source functions (`hitTriageDataSource` role)
  - Demo submit function (`hitTriageSubmitFunction` role)
  - `registerMoleculesToViD` — batch-registers all campaign molecules to V-iD dictionary
  - `Hit Design V-iD` panel — semantic value panel for `HIT_DESIGN_VID` semtype
  - `gasteigerRenderer` — custom grid cell renderer for Gasteiger PNG images
  - Package settings editor widget

### App Classes (in `src/app/`)

| File | Class | Purpose |
|---|---|---|
| `hit-app-base.ts` | `HitAppBase<T>` | Abstract base for all three apps. Handles compute functions discovery, permissions, dataframe union/join, campaign locking, CSV download |
| `hit-triage-app.ts` | `HitTriageApp` | Hit Triage flow: InfoView → PickView (filter) → SubmitView. Uses `DG.MultiView` |
| `hit-design-app.ts` | `HitDesignApp<T>` | Hit Design flow: InfoView → DesignView (grid + sketching) → SubmitView. Handles molecule joining, single-cell calculations, V-iD registration, tile viewer, auto-save |
| `pepti-hit-app.ts` | `PeptiHitApp` | Extends `HitDesignApp` for HELM/peptide molecules. Converts HELM→atomic level via `Bio:toAtomicLevel` |
| `base-view.ts` | `HitBaseView<T,A>` | Base view class extending `DG.ViewBase` with campaign deletion and template access |

### Types — `src/app/types.ts`

Key types:
- `AppName` = `'Hit Triage' | 'Hit Design' | 'PeptiHit'`
- `HitTriageCampaign` — campaign state: name, template, status, filters, ingest config, save path, permissions, layout, column types
- `HitDesignCampaign` — like HitTriageCampaign but without filters/ingest; has tile viewer sketch state
- `HitTriageTemplate` / `HitDesignTemplate` — template definitions with compute config, submit config, campaign fields, stages, layout
- `TriagePermissions` — `{edit: string[], view: string[]}` storing group IDs

### Constants — `src/app/consts.ts`

Important constants and their roles:
- `HitTriageComputeFunctionTag` / `HitDesignerFunctionTag` — tags for discoverable compute functions
- `HitTriageDataSourceTag` / `HitTriageSubmitTag` — tags for data source and submit functions
- `ViDColName = 'V-iD'`, `ViDSemType = 'HIT_DESIGN_VID'` — virtual ID column
- `TileCategoriesColName = 'Stage'` — column used for tile viewer lanes
- `HitDesignMolColName = 'Molecule'`, `PeptiHitHelmColName = 'Helm'`
- `i18n` — UI label strings
- `CampaignGrouping` enum — None, Template, Status, Author, Last Modified User

### Views (per app)

**Hit Triage views** (`src/app/hit-triage-views/`):
- `info-view.ts` — Landing page: list campaigns, create new campaign/template
- `submit-view.ts` — Submit filtered results to a configured function

**Hit Design views** (`src/app/hit-design-views/`):
- `info-view.ts` — Landing page with campaign table, grouping, sorting, create/continue campaigns
- `submit-view.ts` — Submit results
- `tiles-view.ts` — Tile/Kanban viewer for stage-based molecule organization

**PeptiHit views** (`src/app/pepti-hits-views/`):
- `info-view.ts` — Extends Hit Design info view for peptide context

### Accordions (`src/app/accordeons/`)

UI accordions for creating new campaigns and templates:

| File | Purpose |
|---|---|
| `new-campaign-accordeon.ts` | Hit Triage new campaign form: file/query data source, campaign fields |
| `new-template-accordeon.ts` | Hit Triage template creator: name, key, compute config, submit function, layout, campaign fields |
| `new-hit-design-campaign-accordeon.ts` | Hit Design/PeptiHit new campaign: creates initial DataFrame with Molecule, Stage, V-iD columns |
| `new-hit-design-template-accordeon.ts` | Hit Design/PeptiHit template creator with stages editor |
| `layout-input.ts` | Reusable layout file input that parses `.layout` files into `DG.ViewLayout` |

### Dialogs (`src/app/dialogs/`)

| File | Purpose |
|---|---|
| `functions-dialog.ts` | Compute configuration dialog: descriptor tree selection + tagged compute functions/scripts/queries with their parameter editors |
| `save-campaign-dialog.ts` | Simple dialog to name a campaign before saving |
| `permissions-dialog.ts` | Edit view/edit group permissions for a campaign. Default: All Users |

### Utilities

**`src/app/utils.ts`** — General helpers:
- `loadCampaigns()` — loads campaign JSONs from `files/` storage, checks view permissions
- `modifyUrl()` — updates URL query params without page reload
- `checkEditPermissions()` / `checkViewPermissions()` — group-based permission checks
- `addBreadCrumbsToRibbons()` — navigation breadcrumbs in ribbon
- `joinQueryResults()` — joins query result dataframe back into main dataframe by molecule column
- Campaign grouping/sorting helpers for info view tables

**`src/app/utils/calculate-single-cell.ts`** — Runs compute pipeline (descriptors, functions, scripts, queries) on molecule values. Used for both batch column calculations and single-cell recalculations in Hit Design.

**`src/app/utils/molreg.ts`** — V-iD (Virtual ID) molecule registration:
- `obfuscateSmiles()` / `deobfuscateSmiles()` — XOR + base64 obfuscation using package secret key (`meta.sok`)
- `registerMol()` — registers one molecule, returns V-iD string (e.g., `V000001`)
- `registerMolsBatch()` — batch registration in groups of 50
- `registerAllCampaignMols()` — scans all Hit Design campaigns and registers unregistered molecules

**`src/packageSettingsEditor.ts`** — Custom settings editor widget. Settings:
- Default sharing groups (view/edit) for new campaigns
- Default campaign storage folder (default: `System.AppData/HitTriage`)

**`src/pngRenderers.ts`** — `GasteigerPngRenderer`: custom `DG.GridCellRenderer` for base64 PNG images in grid cells.

## Database Schema (`databases/hitdesign/`)

PostgreSQL schema `hitdesign` with two migrations:

**0000_init.sql** — Lock tables:
- `campaign_locks(app_name, campaign_id, expires_at, locked_by)` — 30-second TTL campaign locks for concurrent edit protection
- `update_logs(app_name, campaign_id, updated_at)` — tracks last save time

**0001_dict.sql** — V-iD dictionary:
- `vid_dictionary(id SERIAL, vid VARCHAR(15) GENERATED AS 'V' || LPAD(id, 6, '0'), mh_string TEXT UNIQUE)` — canonical SMILES → auto-generated V-iD
- `campaign_vids(app_name, vid, campaign_id, created_by)` — tracks which V-iDs appear in which campaigns

## SQL Queries (`queries/`)

**vid.sql** (connection: `HitTriage:hitdesign`):
- `addMolecule` — upsert single molecule, return V-iD
- `addMolecules` — batch upsert via `UNNEST`, preserves input order
- `getMoleculeByVid` — lookup canonical SMILES by V-iD
- `getCampaignsByVid` — find all campaigns containing a V-iD

**locks.sql** (connection: `HitTriage:hitdesign`):
- `acquireCampaignLock` — acquire 30-second lock (auto-cleans expired)
- `releaseCampaignLock` — release lock and log update timestamp
- `getLastModified` — get last save timestamp for a campaign

## Files Storage (`files/`)

Campaign and template data stored in `System.AppData/HitTriage/`:

```
files/
  Hit Triage/campaigns/{ID}/campaign.json + enriched_table.csv
  Hit Triage/templates/{Name}.json
  Hit Design/campaigns/{ID}/campaign.json + enriched_table.csv
  Hit Design/templates/{Name}.json
  PeptiHit/campaigns/{ID}/campaign.json + enriched_table.csv
  PeptiHit/templates/{Name}.json
```

Each campaign has a `campaign.json` (metadata, template snapshot, permissions, layout) and `enriched_table.csv` (the molecule data).

## Extension Points

The package discovers external functions via tagged roles:

| Tag/Role | Purpose | Expected Signature |
|---|---|---|
| `HitTriageFunction` | Compute function for property calculation | `(df: DataFrame, molCol: Column, ...args) → void` (mutates df) |
| `HitDesignerFunction` | Same but specific to Hit Design context | Same as above |
| `hitTriageDataSource` | Data source for Hit Triage campaigns | `(...args) → DataFrame` (with molecule column) |
| `hitTriageSubmitFunction` | Submit handler called on campaign submission | `(df: DataFrame, molecules: string) → void` |

## Package Settings (in `package.json`)

- `properties.view` / `properties.edit` — default user group IDs for new campaign sharing
- `properties.defaultCampaignFolder` — storage path (default: `System.AppData/HitTriage`)
- `meta.sok` — secret key for SMILES obfuscation in V-iD dictionary

## Key Dependencies

- `@datagrok-libraries/compute-utils` — compute pipeline utilities
- `@datagrok-libraries/utils` — `u2.appHeader`, `ItemsGrid` UI components
- `uuid` — campaign ID generation
- `typeahead-standalone` — autocomplete UI

## Quick Lookups

| Looking for... | Check first |
|---|---|
| App registration, tree browsers, demo functions | `src/package.ts` |
| Campaign types, template types | `src/app/types.ts` |
| UI labels, column names, tag constants | `src/app/consts.ts` |
| Hit Triage campaign flow | `src/app/hit-triage-app.ts` |
| Hit Design campaign flow (+ PeptiHit base) | `src/app/hit-design-app.ts` |
| Peptide-specific logic | `src/app/pepti-hit-app.ts` |
| Campaign load/save, permissions, URL management | `src/app/utils.ts` |
| Molecule registration (V-iD) | `src/app/utils/molreg.ts` |
| Property calculation pipeline | `src/app/utils/calculate-single-cell.ts` |
| Compute function picker dialog | `src/app/dialogs/functions-dialog.ts` |
| Campaign creation UI | `src/app/accordeons/` |
| Database schema | `databases/hitdesign/` |
| SQL queries (V-iD, locks) | `queries/` |
| Package settings editor | `src/packageSettingsEditor.ts` |
