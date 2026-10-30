# CLAUDE.md

## Overview

**HitTriage** (`@datagrok/hit-triage`) is a Datagrok plugin for chemists to assess hit compound quality and manage
molecule campaigns. It bundles **three applications** in one package:

- **Hit Triage** — Upload a molecule dataset, filter/compute properties, submit results to a chosen function, share campaigns.
- **Hit Design** — Sketch molecules in a spreadsheet, calculate properties, organize into stages, share campaigns.
- **PeptiHit** — Like Hit Design but for peptides (HELM notation). Converts sequences to atomic level via the Bio package.

Category: **Cheminformatics**. Browse paths: `Chem` (Hit Triage, Hit Design), `Peptides` (PeptiHit).

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
- `submit-view.ts` — `HitDesignSubmitView`: the Submit dialog (u2: read-only summary rows, Name and File Path edited in
  place, Status suggestions, Submit function; SAVE / SUBMIT) and `submit()` / `submitWithStatus()` the AI calls
- `tiles-view.ts` — `showProgressTracker(app, host)`: the platform Tile Viewer in a u2 dialog pinned to the campaign
  view (`host.show`, modal; the session's `openTracker` keeps one open), one lane per stage (then one per other stage a
  row holds); a move applies at once, and every close (CLOSE, ✕, Escape, the view closing) keeps the moves and saves
  them — a close-only dialog (`closeOnly`: CLOSE, no CANCEL); Modify Stages (u2 `RecordsInput`, through
  `host.show` too). The card layout kept with the campaign
  (`tilesViewerFormSketch`) is pointed at the campaign's own table (`state.table`) before the viewer reads it

### The campaign view (`src/app/shared/campaign-session.ts`, u2)

A campaign opens as a platform `TableView`; `HitDesignApp.openView` adds it and opens its session with
`openCampaignSession(app, view, options)`: a plain `CampaignSession` (`view`, `host`, `needsRefresh`, the save
feedback, `openTracker`) built under the view's host (`ViewHost.of(view)`), so everything the app adds to the view ends
when the view closes. Its dialogs (merge, compute, designer, a dropped file, the tracker, Submit) are shown through the
host (`host.show`): hidden while another view is current, closed with the view.
- Ribbon through the host (`host.ribbon`): the package's group (refresh with the "new changes" hint, add row, merge,
  compute, designer, permissions, apply template layout, PROGRESS TRACKER, SUBMIT...) pinned to the ribbon's end and
  folding into "…" when narrow; the crumb "Hit Design › <campaign>" is the shell's (`breadcrumbs`: `tv.parentView` +
  `showParentCrumb`), a press on it also closes the view and rebuilds the landing; Save off (`hideSave`).
- The grid's edits through `host.signals.grid` (a layout replaces the grid; an effect re-subscribes): a Stage edit
  saves, an entry edit registers the V-iD and runs the single-cell calculation (`registerFirst`: Hit Design registers
  first, PeptiHit converts the HELM first).
- A Stage cell is a dropdown of the template's stages: the column's `.choices` (`DG.TAGS.CHOICES`), set on open
  (`_prepareColumns`) and by `setStages`, never saved (left out of `columnTags`; the CSV and the layout hold no
  tags). A row in another stage keeps it, in its dropdown too, until another is picked. The template editor's and
  Modify Stages' stage names stay free text.
- The check for other sessions' saves every 3 s while the view is current (a `host.task`: it starts once the last
  save is read, also when that read failed, and never for a session closed meanwhile); save feedback (the grid's "Saving
  Campaign...", the status bar while waiting for another session's lock, "View only — changes are not saved").
- A layout's viewers bound to columns the table lacks (a new campaign's computed ones): the columns are added empty
  and marked (`PlaceholderTag`, `utils/merge-table.ts`); a viewer's columns are read as the platform reads them
  (`…ColumnName(s)`, the Filters viewer's filters); a merge gives one the incoming type, a calculation its own type
  (a number placeholder keeps integers in place, so a filter on it stays), a save leaves it out until it holds values.
- "Calculate additional properties": the compute editor (`computeEditor`, `shared/compute-editor.ts`) in a u2 dialog
  (`HitDesignApp.editCompute`); what is new or changed is calculated over the table (statistics `calculateColumns`:
  the functions get a copy of the table's listed columns, the canonical SMILES under the molecule column's own name;
  what they add or change fills the columns of that name, and the table fires its values changed).
- Hidden `~` columns (Chem's `~Molecule.canonicalSmiles` cache) are never the molecule column
  (`HitAppBase.moleculeColumn(s)`), are dropped when a campaign's table opens and are never written to its CSV or
  campaign.json; Merge with table neither offers nor picks one of the incoming table, and no merge copies one into
  the campaign (`getSourceMappings` reads the incoming columns through u2/dg `Platform.listedColumns`, as the AI
  functions' column lists and Hit Triage's save do).

Hit Design and PeptiHit differ only by `CampaignSessionOptions` (`HitDesignApp.sessionOptions`, overridden by
`PeptiHitApp`) and `prepareTable` (PeptiHit's HELM tags).

Saving (`HitDesignApp.saveCampaign`): one save at a time, a save of a campaign asked for meanwhile folded into the next
one of that campaign; a save writes the campaign open when it was asked for, also when another one opened since (the
app keeps what a pending save needs of the campaign it switched away from; `setTemplate` reads everything first and
takes the new campaign at once). Refused before the lock without edit rights; under the lock a save over another
session's newer save merges that session's table (at the `savePath` its campaign.json names) in first
(`mergeServerTable` in `utils/merge-table.ts`: new rows and columns added, a cell this session left as it read it
takes the other session's value, a molecule several rows hold keeps ours); campaign.json keeps keys this version does
not know; a failure lets the lock go (without logging a save when nothing was written).

**PeptiHit views** (`src/app/pepti-hits-views/`):
- `info-view.ts` — Extends Hit Design info view for peptide context

### Accordions (`src/app/accordeons/`)

UI accordions for creating new campaigns and templates:

| File | Purpose |
|---|---|
| `new-campaign-accordeon.ts` | Hit Triage new campaign form: file/query data source, campaign fields |
| `new-template-accordeon.ts` | Hit Triage template creator: name, key, compute config, submit function, layout, campaign fields |
| `new-hit-design-campaign-accordeon.ts` | Hit Design/PeptiHit new campaign: the initial one-row DataFrame (Molecule, Stage, V-iD; Helm first for PeptiHit); Molecule's units are molblock (the sketcher writes the molecule as drawn, campaign.json keeps the units in `columnTags`); an existing campaign keeps its own |
| `new-pep-triage-template-accordeon.ts` | PepTriage template creator: Hit Triage's plus the sequence and molecule columns |
| `layout-input.ts` | The platform file input a template's layout is picked with; `readLayout` (a `.layout` file's view state and its share path), `layoutFile` (a template's layout as a file the input shows) |

### Shared landing modules (`src/app/shared/`, u2)

The four landings are one composition, configured per app by `HitLandingView` (`base-view.ts`; `init` builds the home
under `Scope.root('hit-home')`: the landing is a `DG.MultiView` child, which no view host ends). Each piece is a u2
function component (`component(…)`, its handle an interface of the same name), except the catalog, a builder run under
the home's owner; their former statics are module functions:

| File | Purpose |
|---|---|
| `hit-home.ts` | `hitHome(options)` → `HitHome` (`catalog`, `picker`): the app header, then a u2 `Page` (in a plain block, so the view scrolls the whole landing) of two sections: "Continue Campaign" with the table and its icons, "New Campaign" with the picker (720px, `PICKER_PX`) |
| `campaign-catalog.ts` | `campaignCatalog(options)` → `CampaignCatalog` (`table`, `chrome()`, `reload(fresh)`): the u2 `CatalogTable` of the campaigns — columns, groupings, row actions, the last-save order, edit rights looked up once per group; `catalogColumns`, `catalogGroupings`, the HD*/HT* localStorage keys (`catalogPersistence`), Date-field values in every stored shape (`fieldValue`, `fieldText`, `dayText`) |
| `template-picker.ts` | `templatePicker(options)` → `TemplatePicker` (`names`): the Template choice with New / Clone / Delete, the form of the picked template (stale picks dropped), the "no templates" state |
| `new-campaign-form.ts` | `newCampaignForm(options)` → `NewCampaignForm` (`form`; Hit Design, PeptiHit): Campaign Name + the template's fields, START disabled with the reason |
| `permissions-editor.ts` | `PermissionsEditor`: Edit Permissions as a u2 `AccessGrid` over the stored `{view, edit}` group ids; `defaultPermissions` |
| `template-editor.ts` | `openTemplateEditor` + `templateEditor(options)` → `TemplateEditor` (`dirty`; Hit Design, PeptiHit): the "New Template" view (`ViewHost.create`; parent crumb, which closes it after a Discard? once something was typed; Templates help), a u2 `Page` (820px, its own label columns): Details / Additional fields / Stages (u2 `RecordsInput`) / Compute / Submit over the page's pinned footer, CREATE disabled with the reason; a clone's free name and key; writes `<app>/templates/<name>.json` in the keys' historical order |
| `compute-editor.ts` | `computeEditor(options)` → `ComputeEditor` (`ready`, `edited`, `read()`): a template's compute as a u2 `MasterDetail` (functions with checkboxes; descriptors as a searchable tree with per-group counts; a function's `funcForm` and re-run flag, off with the reason while unchecked — every detail built at once, so `funcForm` writes its computed defaults as the compute dialog does; "Reorder…" over `HTFunctionOrderingLS`); `read()` answers what the compute dialog writes (`templateCompute`). `descriptors` is optional (no Descriptors row; the given ones kept). No HitTriage import but its CSS; for Chem MPO it still lacks the dialog's single choice (`singleSelect`) |

Hit Triage and PepTriage keep their own new-campaign form (`accordeons/new-campaign-accordeon.ts`) and template editor
(`accordeons/new-template-accordeon.ts`, `new-pep-triage-template-accordeon.ts`) under the picker.

### Dialogs (`src/app/dialogs/`)

| File | Purpose |
|---|---|
| `functions-dialog.ts` | Compute configuration dialog (Hit Triage campaigns, the Hit Triage / PepTriage template editors): descriptor tree selection + tagged compute functions/scripts/queries with their parameter editors; `rerunOnOpenOption`, `descriptorsTree` |
| `save-campaign-dialog.ts` | Simple dialog to name a campaign before saving |
| `merge-table-dialog.ts` | Merge with table (u2): a file from the computer or a share (typed, or picked in the platform's "Select a file" through u2/dg `fileInput`'s "Add file from Files") merged by canonical SMILES or V-iD; a share file is kept in `mergeConfig.filePath` and, with Auto-merge, merged again on every open (`autoMergeOnCampaignOpen`, no save); the column choices are the incoming columns the platform lists (u2/dg `Platform.listedColumns`: no hidden `~` or nameless one) and the defaults are guessed among them (`guessMoleculeColumn`, `guessVidColumn`); each option off says why; OK off with the reason in the footer until the merge has what it needs |

### Utilities

**`src/app/utils.ts`** — General helpers:
- `loadCampaigns()` — loads campaign JSONs from `files/` storage, checks view permissions
- `modifyUrl()` — updates URL query params without page reload
- `checkEditPermissions()` / `checkViewPermissions()` — group-based permission checks
- `addBreadCrumbsToRibbons()` — navigation breadcrumbs in ribbon
- `joinQueryResults()` — joins query result dataframe back into main dataframe by molecule column

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
- `getLastModified` — get last save timestamp for a campaign; `getLastModifiedAll` — every campaign's of an app (the landing's order)
- `dropCampaignLock` — drop the caller's lock (`lockedBy` = `HitAppBase.lockHolder`, user + app instance) without logging a save
- `deleteCampaignRecords` — a deleted campaign's update_logs, campaign_locks and campaign_vids rows

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
- `@datagrok-libraries/u2` — the landings (`src/app/shared/`): `AppHeader`, `Page`, `CatalogTable`, `Form`,
  `MasterDetail`, `AccessGrid`, `confirmDialog`; the package's CSS styles only its own classes (u2's look through
  `--u2-*` properties set on them)
- `@datagrok-libraries/utils` — `ItemsGrid` UI components
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
