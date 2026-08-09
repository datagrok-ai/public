# Grit changelog

## v.next

* GROK-20298: Rewrote the app on domain-ui/the domain UI facade — the hand-rolled `GritApp` view (grid, detail pane, dialogs, timeline) is deleted; the app is now 'Issues' (bug icon) in the browse tree: the main view is the platform Domain View over `grit.issue` with the filter panel open, with 'Create Issue' (`domains.formView`) and 'Projects' (`domains.table('grit.project').app()`) as child tree nodes carrying `Issues / <name>` breadcrumbs and their own routes (`/apps/Grit/create-issue`, `/apps/Grit/projects` — published on open, restored on cold deep links via the app's `path` input)
* GROK-20298: The app acquires the typed schema handle ONCE (`gritUiDb()` from the `grok api --ui` codegen, cached per session) and builds every child view from it — `db.issues.form()` / `db.projects.app()` with columns and insert payloads compile-checked; the data clients are the plural `gritDb.projects` / `gritDb.issues` / ... after the codegen rename
* GROK-20298: `GritIssueHandler` now extends `DG.DomainObjectHandler` — keeps only the genuine customizations (badges on cards/tooltips/markup, search-handle resolution, a `renderGrid` column tweak); Entity View, context panel, and commands fall through to the platform defaults
* GROK-20604: Migrated to the fluent typed surface — builder queries with bound conditions (string-built filters gone), `expand('details:comment')` types the comment array (hand-written ExpandedIssue deleted), `updateWithRetry` replaces the read-then-write patch, `nextNumber` uses `aggregate max(number)`, handle resolution uses business-key `getByKey` lookups, audit rendering is typed `DomainAuditEntry`
* GROK-20602: BREAKING regen — `grok api` codegen v2: datetime fields are dayjs (comment timeline sorts by `valueOf()`), `status`/`priority` are literal unions checked against the schema, expand keys and transactions are typed, db.ts clients are lazy
* GROK-20591: Added a declared default filter panel to `issue` (status, assignee, priority, and the project name via FK path) — the Domain View filter panel pre-opens with these
* GROK-20409: Marked `issue.title` as the primary display-name column (`isName`) — issue titles now caption cards, tooltips, and entity views, and become the friendly names of promoted issues
* GROK-20357: Resolved an issue handle (`<PROJECTKEY>-<number>`, e.g. `GRITEST-1`) typed into global search to its issue — registers a `<KEY>-\d+` `grit.issue` detector per project and teaches the ObjectHandler to resolve the handle (project key + number), render a search card, and open the issue's Entity View on click; the handler claims only the `<KEY>-<number>` shape so the `grit.issue:<uuid>` colon form falls through to the generic domain-handle resolver
* GROK-20336: Moved the static issue-badge presentation into a `grit-badge` CSS class (status/priority color stays inline)
* GROK-20336: Added a `renderProperties` override to the `grit.issue` ObjectHandler so the context panel shows badges, issue fields, and the timeline under whole-handler override
* GROK-20336: Registered a `grit.issue` ObjectHandler (status/priority badges, issue timeline) that overrides the generic domain-row rendering in the Domain View, context panel, and Entity View
* GROK-20308: Introduced Grit (GRok Issue Tracker) — reference app for entity-mapped domain schemas: projects/issues/comments/labels manifest, minimal UI over `grok.dapi.domains`, issue CRUD tests
* GROK-20319: Moved the app and tests to the `grok api`-generated typed clients (`src/generated/db.ts`); the issue detail view now loads comments in the same query via `expand: ['details:comment']` instead of a second per-issue query

## 0.0.1 (2026-07-04)
