# RevvitySignalsLink

Datagrok plugin that integrates with **Revvity Signals** — a cloud platform for managing scientific data (compounds, batches, assays, registrations) used in drug discovery labs.
Users browse libraries and entity types (assets, batches, ...), build and run searches through
a visual query builder, visualize results with molecule structures, save/load named searches,
and see info panels for compound IDs.

## Architecture

| File | What it holds |
|------|---------------|
| `src/package.ts` | Entry point. All platform-registered functions (`//name:` comments): `init`, `revvitySignalsLinkApp`, `revvitySignalsLinkAppTreeBrowser` (browse tree), `searchEntities`, `searchEntitiesWithStructures`, `getUsers`, `getLibraries`, `getTags`, `searchTerms`, `getTermsForField`, `getStructureById`, `entityTreeWidget` (revvity-id panel), `revvityLabelWidget` (revvity-label panel). |
| `src/revvity-api.ts` | Low-level REST client. Single `request()` with 429-retry (10 attempts, 1s delay), 403 passthrough for `/users`. Exports typed wrappers: `queryUsers`, `queryLibraries`, `search`, `queryEntityById`, `queryMaterialById`, `queryTerms`, `queryTags`, `queryStructureById`. |
| `src/signals-search-query.ts` | Revvity's `$match/$and/$or/$gt/$chemsearch/...` query types + **`convertComplexConditionToSignalsSearchQuery(ComplexCondition)`** — bridges the `@datagrok-libraries/utils` QueryBuilder `ComplexCondition` to the Revvity search body. |
| `src/search-utils.ts` | QueryBuilder UI. `runSearchQuery` builds the final query (prepends `materialsCondition` + `assetTypeEid` + `type` filters, defaults sort to `modifiedAt desc`). `initializeFilters` renders the filters panel with save/load icons. `runSearch` executes and handles "Load all" for results > `MAX_RETURN_ROWS`. |
| `src/view-utils.ts` | Browse-tree navigation. `createViewFromPreDefinedQuery` opens a table preview, `openRevvityNode` walks/expands nodes, `handleInitialURL` parses the app path, `updateView` swaps in a new DataFrame + reapplies layout. `openedView` is the module-cached currently-open view. |
| `src/compounds.ts` | Material helpers. `materialsCondition` (always merged into every search via AND), `addMoleculeStructures` (batched 100-at-a-time with cancellable `TaskBarProgressIndicator`), `getConditionForLibAndType`. |
| `src/libraries.ts` | Module-cached library list. `getRevvityLibraries()` (memoized), `resetRevvityLibraries()`, `getLibrariesWithEntityTypes` (calls `/materials/libraries` + `searchTerms` per library). Also renders the homepage stats table. |
| `src/utils.ts` | `transformData` (Revvity JSON `data[].attributes` → DataFrame, flattens tags into columns, looks up users, infers type from `getPropertiesForLibAndEntityType`), `createRevvityWidgetByCorporateId` + `createRevvityWidgetByEntityId` (info-panel renderers), `reorderColumns` (`FIRST_COL_NAMES` / `LAST_COL_NAMES`), `calculatePropForColumn` (add computed column from panel "+"), `createWidgetByRevvityLabel` (revvity-label input picker). |
| `src/properties.ts` | `DG.Property` defaults shown in every QueryBuilder (`createdAt`, `modifiedAt`, `createdBy`, `editedBy`, `Structure`). `REVVITY_USER` semtype + `RevvityUserConditionEditor` (typeahead by first/last name). `REVVITY_FIELD_TO_PROP_TYPE_MAPPING` (`double`/`text`/`date`/`boolean` → `DG.TYPE.*`). `NOT_IN_TAGS` — fields that do NOT get `"in": "tags"` added in queries. |
| `src/users.ts` | Module-cached Revvity users with 403 fallback (`getUsersAllowed = false`; users then extracted from `response.included` in `searchEntities`). |
| `src/layout.ts` | Per-`libName\|entityType` grid layout (`look` + per-column tags) in `grok.userSettings` (`LAYOUT_STORAGE`). Debounce-saved on `onPropertyValueChanged` / `onMetadataChanged` in `updateView`. |
| `src/detectors.ts` | `convertIdentifierFormatToRegexp` turns Revvity numbering formats like `DGS-{#######}-{###}` into regexes, registered as `revvity-id` in `init`. |
| `src/credentials-utils.ts` | `apiKey` + `apiUrl` from `_package.getCredentials()` (parameter names in `constants.ts`). Both module-cached after first read. |
| `src/constants.ts` | All strings: column names (`MOL_COL_NAME`, `HIDDEN_ID_COL_NAME = '~id'`), storage keys (`STORAGE_NAME`, `LAYOUT_STORAGE`), column-ordering lists, exclusion lists, semtype names, `MAX_RETURN_ROWS = 100`, `REVVITY_SEARCH_RES_TOTAL_COUNT` df-tag name. |

Generated (never edit): `src/package.g.ts`, `src/package-api.ts`.

## Glossary

- **Library** — a Revvity materials library ("Compounds"). ID looks like `assetType:686ecf60e3c7095c954bd94f`.
- **Entity type / compound type** — `asset`, `batch`, `sample`, etc. UI pluralizes via `compoundTypeAndViewNameMapping`: `asset↔assets`, `batch↔batches`. Use `getCompoundTypeByViewName` / `getViewNameByCompoundType`.
- **Tags** — Revvity's per-entity key/value metadata. Most query operators need `"in": "tags"` added to the field clause; the exceptions are in `NOT_IN_TAGS` (system columns like `createdAt`, `type`, `isMaterial`).
- **`revvity-id` semtype** — corporate compound IDs that match a library numbering regex. Triggers `entityTreeWidget` info panel.
- **`revvity-label` semtype** — values from `REVVVITY_LABEL_FIELDS` (currently `materials.project`). Triggers `revvityLabelWidget` — a choice input scoped to the current library/entityType.
- **`REVVITY_USER` semtype** — applied to `createdBy`/`editedBy` QueryBuilder properties; its custom editor typeaheads by name.
- **Saved search** — a `ComplexCondition` JSON. Stored in `grok.userSettings` under `SAVED_SEARCH_STORAGE`, keyed `libName|entityType`, each value is `{searchName: JSON-stringified-condition}`.
- **`funcs` namespace** — `import { funcs } from './package-api'`. Call platform-registered functions this way (not direct TS imports) so `meta.cache: all` and client-cache kicks in. Most `getX` in `package.ts` have `meta.cache: all` + `invalidateOn: 0 0 * * *`.
- **`openedView`** — the package's module-level current table view. `createViewFromPreDefinedQuery` closes it and replaces with a new `addTablePreview`.
- **`currentQueryBuilderConfig`** — `{libId, libName, entityType, qb?}`. Set by `initializeFilters`, read by `createWidgetByRevvityLabel` to know which library to search against.

## Conventions

- **All REST traffic goes through `grok.dapi.fetchProxy`** — never raw `fetch`. Every external call funnels through the single `request()` in `revvity-api.ts`. Don't add new call sites — add a new wrapper function there.
- **Always call cached package functions via `funcs.*`** (from `package-api.ts`). Direct TS imports bypass `meta.cache`. Examples in code: `funcs.getLibraries()`, `funcs.searchEntities(...)`, `funcs.getTermsForField(...)`.
- **Module-level caches must have resets.** `libraries` (reset via `resetRevvityLibraries()`), `users` (no public reset — cleared by reloading), `apiKey`/`apiUrl`, `filterProperties`. The "Refresh" button in the app clears the `getLibraries` client cache + calls `resetRevvityLibraries()`.
- **Prefixed-with-tilde = hidden.** `HIDDEN_ID_COL_NAME = '~id'` makes the column invisible in the grid. `REVVITY_SEARCH_RES_TOTAL_COUNT = '~RevvitySearchResTotalCount'` is a DataFrame tag, not a column, but follows the same convention.
- **Path format.** The app uses `apps/Revvitysignalslink/<lib>/<viewName>[/search/<urlEncodedJSON>]` or `.../saved searches/<lib>/<viewName>/<searchName>`. `createPath` + `handleInitialURL` are the only places that construct / parse it.
- **Query shape.** `runSearchQuery` wraps user conditions as `<sectionCondition> AND assetTypeEid=<libId> AND type=<compoundType> AND <user conditions>`. Today `<sectionCondition>` is always `materialsCondition` (`isMaterial=true AND type!=assetType`) because this plugin only covers the **materials section** of Signals. Revvity Signals also has experiments, notebooks, journals, etc. — supporting those means branching on section, not extending `materialsCondition`. Don't add peer top-level conditions; extend `materialsCondition` or append to the inner `cond.conditions`.
- **QueryBuilder fields come from two places:** `getDefaultProperties()` (hardcoded in `properties.ts`) + tags fetched via `funcs.getTags(entityType, libId)`, merged in `getPropertiesForLibAndEntityType` (cached per `libId|entityType` in `filterProperties`).
