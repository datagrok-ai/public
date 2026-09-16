# Stockroom changelog

## v.next

* GROK-20753: Locations reads as a page: the two panes carry "Locations" / "Containers" headers, the status line says "6 containers in Lab 101, including sublocations", Enter (or a double-click) on a container opens its entity page, and an "All locations" row is the way back to every container once a node has been picked
* GROK-20753: The Stockroom app follows the server (`live: true`) — another session's substance shows up within a poll instead of on the next reload
* GROK-20753: Added the **Locations** app — the location tree beside the containers stored under the selected node: `domains.tree(locations)` drives a `domains.dataTable` over a live container source through `location_id under "<id>"`, and nothing selected is the whole table
* GROK-20753: Declared `location` a hierarchy (`"hierarchy": true`, schema 1.1.0) — the self-referencing `parent_id` is now the registry's parent column, so `pathTo(id)` answers a row's ancestors and the `under` filter term matches a whole subtree
* GROK-20753: Seeded the Pilot plant branch (`0005_pilot_plant_shelves.sql`) — two shelves in the warehouse with a container each, so both roots of the tree carry a subtree
* GROK-20753: Added `Stockroom: trash` — a location deleted, found in the trash through `get(id, {deleted: 'only'})`, restored, and its `insert,delete,undelete` audit trail; `Stockroom: schema` covers the hierarchy flag, `pathTo` and `under` on both the tree table and the ref column into it
* GROK-20753: Fixed the unstyled domain chrome — the breadcrumb, the child tabs, the History section, the filter box, the grid and the value editors now come dressed: the package imports the one `u2/src/dg/domain/styles.js` instead of a hand-kept list of sheets that had drifted
* GROK-20753: Introduced Stockroom — a chemical stockroom on the GHS classification and the zero-code reference app for entity-mapped domain schemas: twelve tables declared in `databases/stockroom/schema.json` (constraints, ref filters, searchable columns, custom permissions, a self-referencing location tree, N:N hazards, a file column, auto-numbered labels), the UNECE GHS vocabulary (29 hazard classes, 80 H-statements, 97 P-statements) plus a small demo stockroom as seed scripts, a three-line `package.ts` over `domains.table(...).app()`, a `dg-ui/1` spec for the configuration tier, and a schema smoke test
