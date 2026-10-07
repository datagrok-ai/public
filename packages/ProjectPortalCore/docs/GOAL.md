# Goal

A project portal where every project team member is seamlessly connected to their project's data and team.

1. Start from your molecule or project, not from raw platform objects.
2. Find any molecule or project from one landing page, with recent and recommended actions.
3. See project status and the next deliverables in one place.
4. Get governed project data in one click, with no spreadsheets and no assembling files by hand.
5. Know what data is ready and what gap is blocking model or document use.
6. Find the right model for your molecule, unit operation, stage and data, explained by the scientific question it answers.
7. See which model runs have already been done.
8. Generate documents and evidence tables from linked data and runs, so slides are the last resort.
9. See which deliverable documents exist, which are missing, and what evidence backs them.
10. Keep discussions, questions and decisions with the project, not in email or chat.
11. Link every decision to its evidence, so the reasoning can be found and reused.
12. Trace any output back to its source data, model version, document and decision.
13. Share live data and analysis instead of saved links and emailed screenshots.
14. Use it on any device, in real time, as a guided simple mode or a technical advanced mode.
15. See only the data, models, documents and decisions your source-system permissions allow.

## Who delivers each point

The core package delivers the capability. A customer package supplies only data and configuration: its registry
and data sources, data rules, model catalogue and cards, gate deliverables, document templates, and the feed of
its own entitlements into groups. The platform delivers what a package cannot.

| # | Core | Customer package | Platform |
|---|---|---|---|
| 1 | Board, program page, workspace builder | Vocabulary, workspace tree | — |
| 2 | Search, recent programs, My programs | — | — |
| 3 | Stage track, gate criteria; a date and an owner per deliverable (needs a core schema field: `criterion_status.owner` and `due_on` exist but nothing sets or shows them, and deliverables outside gate criteria have no table) | Stage model, deliverables per gate | — |
| 4 | Program dataset (F-09) | Data source | — |
| 5 | Readiness and data conditions (F-20) | Data rules | — |
| 6 | Catalogue, model cards, recommendation (F-12, F-13), scoped to the portfolio's own catalogue, by unit operation, with the question shown | Models and cards | — |
| 7 | Saved model runs (F-14) | Installed models | Row grants for non-author writes (GROK-21041) |
| 8 | Evidence tables (F-40); document generation (F-50, new) | Document templates | Row grants for non-author writes (GROK-21041) |
| 9 | Documents and document cards (F-27) | Deliverable list | Row grants for non-author writes (GROK-21041) |
| 10 | Discussion and Ask SME (F-51, new); decisions (F-25) | Where Ask SME is routed | Row grants for non-author writes (GROK-21041); a notification API |
| 11 | Decisions with evidence (F-25) | — | Row grants for non-author writes (GROK-21041) |
| 12 | Provenance (F-28) | — | Soft references to function calls (GROK-21042) |
| 13 | Stable program addresses; a share-link button | — | Sharing, live views; scheduling a package function (new-data notices wait for someone to sync) |
| 14 | Simple and advanced modes (F-08); tablet-first surface (F-45, not built) | Persona mapping | Real-time push; scheduling a package function |
| 15 | Program roles and grants (F-16); keeping each portfolio's records to its own people (requests, sync history, compounds, follows, the catalogue) | Entitlement feed (F-46) | Sharing a single row (GROK-21041) |

GROK-21041 is not only point 15's: until a row can be granted, only a program's author and admins write its gate
log, evidence, documents, runs and observations, so points 7 to 11 work for writers only as admin. Point 12's run
reference waits on GROK-21042. Points 13 and 14 need the platform to run a package function on a schedule and to
push to the browser. Keeping one portfolio's records and catalogue from another's users is core work.
