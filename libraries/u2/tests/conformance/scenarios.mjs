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
    name: 'an unsupported op is refused by name: the member and the flag agree',
    requires: [],
    seed: [],
    async run(table, ids, t) {
      const support = table.support;
      t.ok(support !== undefined, 'the backend declares what it can do');
      for (const [member, flag] of [['restore', 'restore'], ['ancestors', 'ancestors'],
        ['updateWhere', 'writes'], ['batch', 'writes'], ['probe', 'probe'], ['audit', 'audit']]) {
        t.equal(table[member] !== undefined, support[flag] === true,
          `${member} is installed exactly when support.${flag} says so`);
      }
      t.ok(support.systemColumns.includes('id'), 'the projection always carries the id');
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
  {
    name: 'restore: a landed delete is soft, and restore brings the row back',
    requires: ['softDelete', 'restore'],
    seed: [
      {key: 'gone', code: 'res-gone', name: 'Gone'},
      {key: 'kept', code: 'res-kept', name: 'Kept'},
    ],
    async run(table, ids, t) {
      await table.transaction([{op: 'delete', table: table.address, id: ids.gone}]);
      const live = await table.query({limit: 100});
      t.deepEqual(live.map((row) => row.code), ['res-kept'], 'a deleted row leaves a live read');
      t.equal(live[0]['~is_deleted'], undefined, 'which does not project the column at all');
      const trash = await table.query({deleted: 'only', limit: 100});
      t.deepEqual(trash.map((row) => row.code), ['res-gone']);
      t.equal(trash[0]['~is_deleted'], true, 'and a trash read does');
      await table.restore(ids.gone);
      t.deepEqual((await table.query({limit: 100})).map((row) => row.code).sort(),
        ['res-gone', 'res-kept'], 'restore puts it back');
      await t.rejects(() => table.restore(ids.gone), 'not-found',
        'restoring a live row is not found, not a no-op');
    },
  },
  {
    name: 'restore: a deleted parent vetoes a child\'s restore, and the veto names the column',
    requires: ['softDelete', 'restore', 'hierarchy'],
    seed: [
      {key: 'parent', code: 'veto-parent', name: 'Site'},
      {key: 'child', code: 'veto-child', name: 'Room', parent_id: '$ref:parent'},
    ],
    async run(table, ids, t) {
      await table.transaction([{op: 'delete', table: table.address, id: ids.child}]);
      await table.transaction([{op: 'delete', table: table.address, id: ids.parent}]);
      let code = null;
      let message = '';
      try {
        await table.restore(ids.child);
      } catch (e) {
        code = e.code;
        message = `${e.message}`;
      }
      t.equal(code, 'restrict', 'a child under a deleted parent is refused, not restored');
      t.ok(message.includes('parent_id'), 'and the refusal names the ref column');
      await table.restore(ids.parent);
      await table.restore(ids.child);
      t.deepEqual((await table.query({limit: 100})).map((row) => row.code).sort(),
        ['veto-child', 'veto-parent'], 'the parent first, then the child');
    },
  },
  {
    name: 'updateWhere: the filter is required, and an unwritable column is refused before any write',
    requires: ['updateWhere'],
    seed: [
      {key: 'a', code: 'req-a', name: 'First'},
      {key: 'b', code: 'req-b', name: 'Second'},
    ],
    async run(table, ids, t) {
      await t.rejects(() => table.updateWhere('', {name: 'Nope'}), 'validation');
      await t.rejects(() => table.updateWhere('code = "req-a"', {}), 'validation');
      await t.rejects(() => table.updateWhere('code = "req-a"', {id: ids.b}), 'validation');
      const rows = await table.query({sort: 'code', limit: 100});
      t.deepEqual(rows.map((row) => row.name), ['First', 'Second'], 'a refusal writes nothing');
      t.deepEqual(rows.map((row) => row.version), [1, 1]);
    },
  },
  {
    name: 'updateWhere: deleted rows are never touched',
    requires: ['updateWhere', 'softDelete'],
    seed: [
      {key: 'live', code: 'upd-live', name: 'Live'},
      {key: 'gone', code: 'upd-gone', name: 'Gone'},
    ],
    async run(table, ids, t) {
      await table.transaction([{op: 'delete', table: table.address, id: ids.gone}]);
      const done = await table.updateWhere('code in ("upd-live", "upd-gone")', {name: 'Touched'});
      t.equal(done.updated, 1, 'the trashed row is not a live row');
      const [gone] = await table.query({deleted: 'only', limit: 100});
      t.equal(gone.name, 'Gone', 'and was not written');
    },
  },
  {
    name: 'batch: duplicates are skipped and reported; errorOnDuplicate makes them errors; upsert merges',
    requires: ['batch'],
    seed: [{key: 'a', code: 'bat-a', name: 'Existing'}],
    async run(table, ids, t) {
      const report = await table.batch([{code: 'bat-a', name: 'Again'}, {code: 'bat-b', name: 'New'},
        {code: 'bat-b', name: 'Twice'}], {allOrNothing: false});
      t.equal(report.inserted, 1);
      t.equal(report.updated, 0);
      t.equal(report.skipped, 2, 'the live clash and the second occurrence inside the batch');
      t.equal(report.errorCount, 0);
      t.deepEqual([...report.rows].sort((x, y) => x.index - y.index).map((row) => row.status),
        ['duplicate', 'inserted', 'duplicate']);
      t.equal([...report.rows].find((row) => row.index === 0).existingId, ids.a,
        'a live clash names the row it clashed with');

      const strict = await table.batch([{code: 'bat-a', name: 'Again'}],
        {allOrNothing: false, errorOnDuplicate: true});
      t.equal(strict.errorCount, 1);
      t.equal(strict.skipped, 0);
      t.equal(strict.rows[0].errors[0].code, 'unique');

      const merged = await table.batch([{code: 'bat-a', name: 'Renamed'}], {mode: 'upsert'});
      t.equal(merged.updated, 1);
      t.equal(merged.inserted, 0);
      const [a] = await table.query({filter: 'code = "bat-a"', limit: 10});
      t.equal(a.name, 'Renamed', 'upsert merged into the row the key names');

      // the payload is checked before any row is: a column that is not there, one the caller may
      // not write, and an upsert with no business key to merge by
      await t.rejects(() => table.batch([{nope: 1}]));
      await t.rejects(() => table.batch([{id: ids.a, name: 'x'}]));
      await t.rejects(() => table.batch([{name: 'No key'}], {mode: 'upsert'}));
    },
  },
  {
    name: 'batch: allOrNothing writes nothing and reports `error`; false applies the good rows',
    requires: ['batch'],
    seed: [],
    async run(table, ids, t) {
      const payload = [{code: 'aon-a', name: 'Fine'}, {name: 'No code'}];
      const aborted = await table.batch(payload);
      t.equal(aborted.error !== undefined && aborted.error !== null, true, 'the abort is reported');
      t.equal(aborted.errorCount, 1);
      t.deepEqual(aborted.rows.map((row) => [row.index, row.status]), [[1, 'error']]);
      t.equal((await table.query({limit: 100})).length, 0, 'and nothing was written');

      const partial = await table.batch(payload, {allOrNothing: false});
      t.equal(partial.inserted, 1);
      t.equal(partial.errorCount, 1);
      t.deepEqual((await table.query({limit: 100})).map((row) => row.code), ['aon-a']);
    },
  },
  {
    name: 'validate: the dry run\'s verdicts equal the commit\'s, and it writes nothing',
    requires: ['validate'],
    seed: [{key: 'a', code: 'val-a', name: 'Existing'}],
    async run(table, ids, t) {
      const payload = [{code: 'val-a', name: 'Dup'}, {code: 'val-b', name: 'New'}, {name: 'No code'}];
      const before = (await table.query({limit: 100})).length;
      const dry = await table.validate(payload, {allOrNothing: false});
      t.equal(dry.validateOnly, true);
      t.equal(dry.rowCount, 3);
      t.equal((await table.query({limit: 100})).length, before, 'a dry run writes nothing');
      t.deepEqual([...dry.rows].sort((x, y) => x.index - y.index).map((row) => row.predicted),
        ['skip', 'insert', 'error']);
      t.equal(dry.willInsert, 1);
      t.equal(dry.willUpdate, 0);
      t.equal(dry.willSkip, 1);
      t.equal(dry.errorCount, 1);

      const report = await table.batch(payload, {allOrNothing: false});
      const landed = {inserted: 'insert', updated: 'update', duplicate: 'skip', error: 'error'};
      t.deepEqual([...report.rows].sort((x, y) => x.index - y.index).map((row) => landed[row.status]),
        [...dry.rows].sort((x, y) => x.index - y.index).map((row) => row.predicted),
        'the prediction is the commit, row for row');
    },
  },
  {
    name: 'captions: a ref column carries the target\'s display name; an unknown one is one refusal',
    requires: [],
    seed: [
      {key: 'parent', code: 'cap-parent', name: 'Site'},
      {key: 'child', code: 'cap-child', name: 'Room', parent_id: '$ref:parent'},
    ],
    async run(table, ids, t) {
      const rows = await table.query({captions: ['parent_id'], sort: 'code', limit: 100});
      t.deepEqual(rows.map((row) => `${row.code}=${row['~caption_parent_id']}`),
        ['cap-child=Site', 'cap-parent=null'],
        'the target\'s display name, and null where there is no target');
      await t.rejects(() => table.query({captions: ['name']}), 'filter');
      await t.rejects(() => table.query({captions: ['nope']}), 'filter');
    },
  },
  {
    name: 'probe: an unscoped read answers a token that moves once per write transaction',
    requires: ['probe'],
    seed: [{key: 'a', code: 'prb-a', name: 'First'}],
    async run(table, ids, t) {
      const first = await table.probe();
      t.equal(first.count, -1, 'nothing is counted: the whole live table is the change token');
      await table.transaction([{op: 'update', table: table.address, id: ids.a, values: {name: 'Second'}}]);
      const second = await table.probe();
      t.equal(second.count, -1);
      t.equal(second.last !== first.last, true, 'and the token moved');
      const scoped = await table.probe({filter: 'code = "prb-a"'});
      t.equal(scoped.count, 1, 'a scoped read still counts rows');
    },
  },
];
