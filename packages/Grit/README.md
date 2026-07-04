# Grit

Grit (GRok Issue Tracker) is the reference app for **entity-mapped domain schemas**
(`databases/grit/schema.json`): plugin-owned relational tables inside the Datagrok database
with platform row/column security, Datlas-managed CRUD, and an in-transaction audit trail.

It demonstrates the relational-app story:

* `project` — `table`-mode security, business key on `key`
* `issue` — `row` mode with lazy promotion (an issue becomes shareable as an entity on first
  share), per-project numbering (`GRIT-123`, assigned by the app), `user` columns
  (reporter, assignee), `status`/`priority` choices dictionaries
* `comment`, `issue_label` — `master` mode delegating security to the issue; cascade on
  issue delete
* `label` — `table`-mode dictionary

The app (`Apps | Grit`) is a minimal client of the generic `grok.dapi.domains` JS API:
issue grid per project, detail pane with status/priority editing, comments, and the audit
timeline surfaced as issue history.

See `core/docs/features/db-schemas/ARCHITECTURE.md` (§9.4) in the Datagrok monorepo for the design.
