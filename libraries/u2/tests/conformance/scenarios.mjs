/* What every domain backend must answer the same way (holistic review 4, A4). One module, two
   harnesses: `tests/backend-conformance.test.js` runs it headless over `MemoryTable`, and U2Demo's
   `U2: domain conformance` runs it on the stand over the platform's `DomainTable`. A backend that
   drifts from the other is a red test here, not a surprise in the browser.

   The fixture both harnesses build: table `folder` — `code` (string, required, the business key),
   `name` (string, the name column), `parent_id` (ref → `folder`), `hierarchy: true`, soft delete
   on. Plain ESM importing NOTHING, so a package with no u2 link can run it too. */

/** @type {import('./scenarios.d.ts').ConformanceScenario[]} */
export const scenarios = [
  {
    name: 'ancestors: a deleted seed still answers its chain',
    requires: ['hierarchy', 'softDelete'],
    seed: [
      {key: 'root', code: 'anc-root', name: 'Site'},
      {key: 'mid', code: 'anc-mid', name: 'Room', parent_id: '$ref:root'},
      {key: 'leaf', code: 'anc-leaf', name: 'Shelf', parent_id: '$ref:mid'},
    ],
    async run(table, ids, t) {
      const chain = [{id: ids.root, name: 'Site'}, {id: ids.mid, name: 'Room'}];
      t.deepEqual(await table.ancestors(ids.leaf), chain, 'root first, and without the row itself');
      await table.transaction([{op: 'delete', table: table.address, id: ids.leaf}]);
      t.deepEqual(await table.ancestors(ids.leaf), chain,
        'the seed level takes a deleted row: a row opened from the trash still has its breadcrumb');
      await table.transaction([{op: 'delete', table: table.address, id: ids.mid}]);
      t.deepEqual(await table.ancestors(ids.leaf), [], 'while the walk truncates at a deleted ancestor');
    },
  },
  {
    name: 'under: the subtree of a hierarchy, the root itself included',
    requires: ['hierarchy'],
    seed: [
      {key: 'root', code: 'sub-root', name: 'Site'},
      {key: 'mid', code: 'sub-mid', name: 'Room', parent_id: '$ref:root'},
      {key: 'leaf', code: 'sub-leaf', name: 'Shelf', parent_id: '$ref:mid'},
      {key: 'other', code: 'sub-other', name: 'Other site'},
    ],
    async run(table, ids, t) {
      const codes = async (property) =>
        (await table.query({filter: [{property, operator: 'under', value: ids.root}], limit: 100}))
          .map((row) => row.code).sort();
      t.deepEqual(await codes('id'), ['sub-leaf', 'sub-mid', 'sub-root'],
        'off the id column: the whole subtree, the seed row included');
      t.deepEqual(await codes('parent_id'), ['sub-leaf', 'sub-mid'],
        'off a ref column: the rows pointing into that subtree');
    },
  },
  {
    name: 'updateWhere: a limit of 0 is clamped to one row and says the filter matched more',
    requires: ['updateWhere'],
    seed: [
      {key: 'a', code: 'cap-a', name: 'First'},
      {key: 'b', code: 'cap-b', name: 'Second'},
    ],
    async run(table, ids, t) {
      const filter = 'code in ("cap-a", "cap-b")';
      const capped = await table.updateWhere(filter, {name: 'Capped'}, {limit: 0});
      t.equal(capped.updated, 1, 'a limit of 0 is clamped to one row, never an infinite drain');
      t.equal(capped.hasMore, true, 'and the rest is reported, not silently dropped');
      const rest = await table.updateWhere(filter, {name: 'Capped'});
      t.equal(rest.updated, 2, 'the default cap takes everything the filter matches');
      t.equal(rest.hasMore, false);
    },
  },
];
