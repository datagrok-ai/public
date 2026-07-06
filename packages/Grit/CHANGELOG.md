# Grit changelog

## v.next

* GROK-20357: Resolved an issue handle (`<PROJECTKEY>-<number>`, e.g. `GRITEST-1`) typed into global search to its issue — registers a `<KEY>-\d+` `grit.issue` detector per project and teaches the ObjectHandler to resolve the handle (project key + number), render a search card, and open the issue's Entity View on click; the handler claims only the `<KEY>-<number>` shape so the `grit.issue:<uuid>` colon form falls through to the generic domain-handle resolver
* GROK-20336: Moved the static issue-badge presentation into a `grit-badge` CSS class (status/priority color stays inline)
* GROK-20336: Added a `renderProperties` override to the `grit.issue` ObjectHandler so the context panel shows badges, issue fields, and the timeline under whole-handler override
* GROK-20336: Registered a `grit.issue` ObjectHandler (status/priority badges, issue timeline) that overrides the generic domain-row rendering in the Domain View, context panel, and Entity View
* GROK-20308: Introduced Grit (GRok Issue Tracker) — reference app for entity-mapped domain schemas: projects/issues/comments/labels manifest, minimal UI over `grok.dapi.domains`, issue CRUD tests
* GROK-20319: Moved the app and tests to the `grok api`-generated typed clients (`src/generated/db.ts`); the issue detail view now loads comments in the same query via `expand: ['details:comment']` instead of a second per-issue query

## 0.0.1 (2026-07-04)
