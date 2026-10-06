# MolTrack — confirmed defects

Adversarially verified against the working tree (2026-09-09). Severity is user-facing. IDs MT-xx.
MolTrack is 0.0.2 (pre-release) — several highs are cheap to fix now, before the surface freezes.

## High

- **MT-01 — The "+" column icon attaches values to the wrong rows.**
  `src/widgets/moltrack-property-panel.ts:21-26`: the fetched batch data is mapped to grid rows
  positionally (`batchData[i]`), but `getBatchInfoBySynonym` (`utils/utils.ts:27-43`) returns rows in
  backend order with no key echo and generally FEWER rows than the column (non-matching/null ids
  dropped) — wrong values on rows, `undefined` for the tail, all swallowed by the bare catch
  (:27-29). Also creates an OBJECT-typed column for scalar values. Reachable from every Grok-ID
  context panel (`package.ts:310`).

- **MT-02 — Package init runs the setup queries in arbitrary order; failure corrupts permanently.**
  `src/package.ts:26-34`: when the DB is uninitialized, init lists ALL queries of the moltrack
  connection (including plain SELECTs) and calls each with no ordering guarantee — but
  `insertProperties` depends on `insertAdminUser` and `insertSemanticTypeSynonym`
  (`queries/setup.sql:23-50`); wrong order yields rows with NULL `created_by`/`semantic_type_id`,
  and `on conflict do nothing` makes the damage permanent while `checkDBInitialized` (row-existence
  only) reports success forever after. No error handling around the loop; the `catch (e) {}` around
  `Chem:initChemAutostart` is silent.

- **MT-03 — The app home downloads entire entity tables to display five counts.**
  `src/utils/view-utils.ts:150-174`: `getStatisticsWidget` calls `MolTrack:retrieveEntity` for every
  scope — `GET /v1/{scope}/` with no paging (`moltrack-docker-service.ts:105-113`), full JSON →
  flatten → DataFrame — just to read `rowCount`, on every app-home render. The dominant scalability
  problem in the package. (The `flatten: true` argument passed at :157 is undeclared and ignored.)

- **MT-04 — Schema Save reclassifies every property as 'measured'.**
  `src/views/schema-view.ts:102-135`: `editor.properties` is the FULL fetched schema, not just new
  rows; all of it is re-posted with `property_class: 'measured'` (:122), and `buildPropertyOptions`
  doesn't carry the original class through — one Save reclassifies the seeded DECLARED properties
  (`corporate_compound_id`/`corporate_batch_id`, `queries/setup.sql:44,50`). Only the backend's
  conflict behavior stands between this and schema corruption; the client is unambiguously wrong.

## Medium

- **MT-05 [resolved] — Loaded saved-search/history aggregations are silently not applied.**
  `src/views/search.ts:275`: `loadSearchQuery` reassigns the `aggregations` PARAMETER; the Search/
  Save/validation closures keep the original array (:377, :401, :411, :424) — after loading a saved
  search the UI shows the loaded aggregation rows while the executed query carries the previous
  (often empty) list; stale rows also linger in the captured array.

- **MT-06 [resolved] — Filters on `0` / `false` are sent as null.**
  `search.ts:762`: `value: cond.value ? ... : null` — a numeric `= 0` filter or a boolean "is false"
  condition (directly selectable in the query builder) silently becomes `value: null`.

- **MT-07 [resolved] — Saved-search names with underscores become permanently unloadable.**
  `search.ts:76-77`: `key.split('_')[1]` truncates the name (`high_activity` → `high`); the truncated
  key reads back `''` → `JSON.parse('')` throws → "Failed to load saved search query", and the
  duplicate-name validator is defeated. (Scopes with underscores are excluded from the UI, so that
  half is latent.)

- **MT-08 [resolved] — Registration CSV is built with naive join(',').**
  `registration-entity-base.ts:109-113`: user-entered property values with commas shift every
  following column; newlines break the row — silent payload corruption on registration.

- **MT-09 [resolved] — Opening a registration view hijacks the session's sketcher.**
  `registration-entity-base.ts:52`: `DG.chem.currentSketcherType = 'Ketcher'` overwrites the js-api
  global that every later Sketcher reads and that Chem initializes from the user's saved preference —
  session-wide override just by opening the view (resets on reload; also set after the local
  sketcher was constructed, so even the local effect is order-dependent).
  Resolution: Ketcher is required for registration; the assignment now precedes the sketcher
  construction. The session-wide override is accepted behaviour.

- **MT-10 [resolved] — A search round-trip per sketch edit, uncancelled.**
  `registration-entity-base.ts:42-51`: `Sketcher.onChanged` is undebounced and fires
  `getCorporateCompoundIdByExactStructure` (container search) per edit; out-of-order responses can
  leave `compoundExists`/the Register button matching a stale structure.

- **MT-11 — Bulk-registration join keyed on the first CSV column.**
  `registration-tab.ts:211-222`: INNER join on `columns.names()[0]` with no uniqueness/existence
  validation — a duplicate-valued first column multiplies rows and `createSummary` reports counts
  over the multiplied set.

- **MT-12 [resolved] — Search error path wipes the previous results.**
  `search.ts:567-569`: any throw (including transient backend failures) replaces `tv.dataFrame` with
  an empty frame — results lost, only a toast remains.

- **MT-13 [resolved] — Raw `fetch()` against GitHub for schemas/demo data.**
  `utils/fetch-utils.ts:6,14` (pinned raw.githubusercontent.com URL, `constants.ts:70-71`) — the
  repo-wide hard rule requires `grok.dapi.fetchProxy` for external URLs; the rest of the package
  does it right. Reached via `initDB` (explicit call, not package init) — which also re-registers
  GitHub demo data into the registration DB for whoever calls it.

- **MT-14 — Five independent `openedView` globals fragment view lifecycle.**
  `registration-tab.ts:15`, `registration-entity-base.ts:15`, `registration-view-base.ts:9`,
  `schema-view.ts:11`, plus `openedSearchView` (`search.ts:65`): each `show()` closes only its own
  module's view — navigating Compound → Bulk → Schema accumulates stale views; `EntityBaseView`
  shadows the base `show()` so path updates read a null tracker on base-shown views.

- **MT-15 [resolved] — Grid-layout subscriptions accumulate per search run.**
  `search.ts:940-957`: `grid.onPropertyValueChanged` is re-subscribed on every `updateView` (the
  grid survives dataframe swaps) with no disposal — N full tag-serialization + localStorage writes
  per property change after N searches.

## Low

- **MT-16 [resolved]** — `handleSearchURL` saved-search-not-found branch calls `createSearchView` twice
  (`search.ts:1043-1045`) — the guard makes the second a no-op, but the first is fire-and-forget and
  the function returns `grok.shell.tv` instead of the created view.
- **MT-17 [resolved]** — Package-name casing drift in cross-calls (`'Moltrack:registerBulk'`,
  `'MolTrack:RegisterBulk'` in the generated api, `'MolTrack:search'`) — resolution is
  case-insensitive, but the inconsistency invites grep misses.
- **MT-18 [resolved]** — `search.ts:548` stray `; ;`; `createSearchPanel`'s `ui.onSizeChanged` sub never
  disposed (:325-327, bounded by view count); `waitForGrid` setInterval theoretically uncleared;
  hardcoded `'green'/'red'/'orange'` status-icon colors and `'#7990A5'` history color vs design
  tokens; hand-rolled inline SVG funnel icon.
  Resolution: stray `; ;` removed; size-changed sub kept in a single replaced handle; `waitForGrid` kept as is
  (theoretical only); colours moved to `moltrack-*` classes on design tokens. The SVG
  funnel is kept on purpose (the arrow marks the external filter panel; FA `filter` has no arrow).
- **MT-19** — `getCorporateCompoundIdByExactStructure` uses `IS SIMILAR, threshold 1` as "exact
  match" (tautomer/canonicalization semantics untested) and logs via `console.error`.
- **MT-20** — **Zero test coverage**: 17 registered functions, no tests beyond the harness entry
  point (`package-test.ts` imports nothing); no Playwright. For a registration system (data-entry,
  DB-mutating) this is the single biggest process gap.
