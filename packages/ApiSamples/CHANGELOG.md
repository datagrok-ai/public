# API Samples changelog

## v.next

* GROK-20298: Added dapi/domains/facade-form.js — the domain form plain: `issues.form({values})` mounted in a view with a Save button, everything else from the registry
* GROK-20298: Moved the domain samples into their own `dapi/domains/` folder and dropped the `domains-` prefix (dapi/domains-typed-client.js → dapi/domains/typed-client.js, dapi/domains.js → dapi/domains/crud.js, and so on)
* GROK-20298: Added the rest of the `domains` facade catalog — dapi/domains/facade-app.js (the whole CRUD app in one body line), dapi/domains/facade-list-view.js (a query-scoped list page + the reference picker), dapi/domains/facade-grid.js (an editable grid over a filtered subset, plus a grid's machine surface) and dapi/domains/facade-composed.js (a master-detail page with one dirty state and one gate, then a domain grid live next to a bar chart)
* GROK-20298: Added the `domains` facade samples — dapi/domains/facade-form-dialog.js (an insert/edit form in a dialog, 2 lines), dapi/domains/facade-form-view.js (a routed form page through `grok.shell.preview`, 2 lines), dapi/domains/facade-customization.js (replace one input, add one validator) and dapi/domains/facade-machine-surface.js (what an AI agent sees: `DG.Widget.getAll()` → `getWidgetStatus()` → `props` → `getFunctions()`)
* GROK-20298: Added dapi/domains/app-view.js — the domain app framework in two views: a one-expression DomainAppView list page and the DomainEntityAppView page of one row (rowFrom wrapping a row the query already returned)
* GROK-20298: Added dapi/domains/query-state.js — DomainQuery as the one serializable view state: deep-link URL parameters (one key per filter element) round-tripped back, the silent REST spec, the recorded run, and the query behind a live Domain View
* GROK-20298: Added dapi/domains/embedded-view.js and dapi/domains/row-editor.js — DG.DomainView created programmatically and embedded with a permanent filter, and the platform row dialogs as one-line openers (create/edit/picker/conflict/audit/grants); dapi/domains/async-form.js now opens the REAL lookup picker instead of routing to /domains
* GROK-20298: Added dapi/domains/handler.js — DomainObjectHandler: a per-table handler that overrides only renderCard, plus the reflective members every domain table gets for free (properties, capabilities, row permissions, detail tabs, capability-gated actions, editor)
* GROK-20298: Added dapi/domains/registry-reflection.js — registry reflection (rowProperties/tableInfo/resolveNames) and permission-driven table capabilities
* GROK-20298: Added dapi/domains/render-grid.js — customizing every grid over a domain table's rows (the built-in Domain View grid included) via ObjectHandler.renderGrid
* GROK-20298: Added dapi/domains/async-form.js — the "New Grit issue" async form: FK lookup via ui.typeAhead's callback source, server-provided choices, debounced server-side uniqueness validation, and per-field mapping of server rejections
* GROK-20605: Rewrote the domains samples to the typed patterns — fluent builder + bound-condition helpers (apostrophe teaching), getByKey/first/count/exists, aggregateDf, typed tuple transactions with `opIndex`, typed conflict errors + `updateWithRetry`, schema grants with the no-fan-out truth; NEW dapi/domains/typed-client.js (upsert, save/updateWithRetry, typed errors, bounded deleteWhere cleanup); every remaining per-row delete loop is gone
* GROK-20601: Replaced dapi/domains-dapi2.js (the generated dapi2 domains client is removed at typed-surface parity) with dapi/domains/schema.js — schema lifecycle handle: manifest export + whole-schema audit
* GROK-20600: dapi/domains/batch.js cleanup now uses one `deleteWhere` filtered bulk delete instead of the query + per-row delete loop
* GROK-20591: Added dapi/domains/filters.js — batched facets (filter-panel counts) and shareable saved filter presets on a domain table
* GROK-20572: Added data-access/db/query_builder_domain.js — TableQueryBuilder over a domain schema's virtual connection
* GROK-20316: Added dapi/domains-dapi2.js — querying domain-table rows via the generated grok.dapi2 REST client
* GROK-20315: Added dapi/domains/batch.js, domains/transaction.js, domains/aggregate.js, domains/dataframe.js, domains/idempotency.js — batch upsert, multi-entity transactions, aggregation, queryDf → grid, idempotent retries + optimistic concurrency
* GROK-20307: Added dapi/domains/crud.js — domain-table row CRUD via grok.dapi.domains

## 1.2.1 (2-24-11-24)

* New samples, version bump

## 1.0.12 (2024-10-17)

* Add samples before for scripts
* Sample fixes 

## 1.0.11 (2024-07-09)

* Add samples for HelmInput

## 1.0.7 (2023-07-27)

* Dependency: datagarok-api >= 1.16.0*

### Features

* All examples now have auto tests 
