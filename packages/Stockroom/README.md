# Stockroom

Stockroom is a chemical stockroom — reagents, containers, locations, orders, disposals — built on the
GHS classification (UNECE's hazard classes, pictograms, H-statements and P-statements). It is the
**zero-code reference app** for [entity-mapped domain schemas](../../help/develop/how-to/db/domain-schemas.md):
everything is declared in `databases/stockroom/schema.json`, the data is seeded by SQL scripts next to
it, and the app is three lines:

```ts
export async function stockroomApp(): Promise<DG.ViewBase> {
  return (await domains.table('stockroom.substance')).app();
}
```

Grit (`../Grit`) is the code-tier counterpart — actions, validators, renderers, an app subclass.
`src/app.spec.json` shows the middle tier: the same table as a `dg-ui/1` spec (a `u2-domain-source`,
search, filters, a card list, the form, the child tabs, the history and the Save/Discard buttons)
that the designer edits and the app framework renders without a line of TypeScript.

## Tables and what each demonstrates

| Table | Declared | Demonstrates |
|---|---|---|
| `substance` | `name` (isName) and `cas` searchable and unique, `smiles` with `semType: Molecule`, `notes` as a textarea, `hazards` N:N relation, relation-path filters | searchable columns, molecule rendering, declared filter panel |
| `hazard_class`, `h_statement`, `p_statement` | `code` business keys, `h_statement.hazard_class_id` ref and `pictogram` choices, `p_statement.kind` choices | lookups in natural order, the real GHS taxonomy |
| `substance_hazard` | `master` mode delegating to `substance_id`, business key over both refs, no audit | the junction behind the `hazards` relation |
| `container` | auto-numbered `label` (from 1000), `owner` user column, `received`/`expires` dates, constraints `expires >= received` and `quantity <= initial_quantity` with messages, `location_id` picker filtered by `site = $site` | child-collection tabs, user column, autoNumber, CHECK constraints the form pre-validates, dependent picker |
| `location` | `kind` choices, `parent_id` self-reference, `site` | self-referencing lookup (a plain ref — see follow-ups) |
| `vendor`, `purchase_order`, `order_line` | auto-numbered order `number`, `status` with a default, `permissions: {approve}`, `order_line` in `master` mode under its order | classic master–detail; a draft parent referenced by its children saves in one transaction |
| `sds_document` | `file` column, `revision`, `language` choices, `master` mode under the substance | file columns, per-row history |
| `disposal` | `container_id` ref, `reason` textarea, `approved_by` user column, `permissions: {approve}` | custom permission gating an action, audit trail |

## Seed data

`databases/stockroom/000N_*.sql` run once, in name order, after the schema deploys (recorded in
`package_db_ups`; a **debug** publish skips them — publish with `--release` to land the data). Every
script is idempotent: ids are `md5('<table>:<natural key>')::uuid` and every insert is
`ON CONFLICT (id) DO NOTHING`.

| Script | Rows |
|---|---|
| `0001_hazard_classes.sql` | 29 GHS hazard classes (2.1–2.17, 3.1–3.10, 4.1–4.2; UNECE GHS Rev.10 names) |
| `0002_h_statements.sql` | 80 hazard statements — the current single-code set with its hazard class and pictogram; statements withdrawn by Rev.8 (H200–H203, H205) and the EU CLP sub-variants (H350i, H360D/F…, H361d/f…) are out |
| `0003_p_statements.sql` | 97 precautionary statements by kind (general, prevention, response, storage, disposal); combined statements are label compositions and are out |
| `0004_demo_data.sql` | 8 locations on two sites, 1 vendor, 20 common lab reagents with CAS, SMILES and 71 hazard links, 7 containers, 1 draft purchase order with 3 lines; the auto-number counters are moved past the seeded labels |

## Follow-ups (need platform work, not faked here)

* `location` as a real hierarchy: the phase-3 `hierarchy` key (tree navigation, ancestor filters) does
  not exist yet — `parent_id` is a plain self-reference.
* Soft-delete restore UI: rows are soft-deleted by the engine but there is no restore action in the app.
* Row visibility by site through grants (a user sees only their site's containers): needs row-level
  grants keyed on a column, phase 3.
* `disposal.approve` and `purchase_order.approve` are declared and grantable, but the zero-code app
  has no "Approve" action bound to them — that is the code tier (Grit's `actions.add`).
* Molecule rendering of `smiles` relies on the Chem package being installed for the `Molecule`
  semantic type.
