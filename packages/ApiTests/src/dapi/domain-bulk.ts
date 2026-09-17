import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// Bulk edit is one call: `updateWhere(filter, values)` patches every row the filter
// matches AND the caller may edit, oldest first, in ONE transaction, capped at 1000
// rows (`options.limit` lowers it, `hasMore` says to loop). It rides the same per-row
// engine as update — writability, immutability, validation, audit — so a value the
// engine refuses refuses the WHOLE call and nothing is written. Fixture: apitests.item;
// every test cleans its own sku prefix up.
category('Dapi: domain bulk', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const like = (property: string, prefix: string): any =>
    ({property, operator: 'like', value: `${prefix}%`});

  async function cleanup(prefix: string): Promise<void> {
    try {
      await items().deleteWhere(like('sku', prefix));
    } catch (e) {
      console.error(`bulk fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  test('updateWhere over id in (…): the named rows only, one version bump and one audit entry each', async () => {
    const prefix = `bulk-sel-${stamp()}`;
    const rows = await items().insert([0, 1, 2].map((i) =>
      ({sku: `${prefix}-${i}`, name: 'Bulk', quantity: 1})));
    try {
      const picked = [rows[0].id, rows[1].id];
      const before = (await items().get(picked[0])).version;
      const report = await items().updateWhere(
        `id in (${picked.map((id) => `"${id}"`).join(', ')})`, {quantity: 7});
      expect(report.updated, 2, `updateWhere reported ${JSON.stringify(report)}`);
      expect(report.hasMore, false, 'two rows under the cap must not report hasMore');

      const after = await items().query({filter: like('sku', prefix), sort: 'sku'});
      expect(after.map((r: any) => r.quantity).join(','), '7,7,1',
        'the update did not land exactly on the selected rows');
      expect(after[0].version, before + 1, 'the bulk update did not bump the version once');
      expect(after[2].version, before, 'a row outside the selection was written');

      const audit = await items().audit(picked[0]);
      expect(audit.filter((a) => a.op === 'update').length, 1,
        `one update entry per row expected: ${audit.map((a) => a.op).join(', ')}`);
    } finally {
      await cleanup(prefix);
    }
  });

  test('limit caps the call and hasMore drives the loop', async () => {
    const prefix = `bulk-cap-${stamp()}`;
    await items().insert([0, 1, 2].map((i) => ({sku: `${prefix}-${i}`, name: 'Capped', quantity: 1})));
    const pending = (): any => [like('sku', prefix), 'and', {property: 'quantity', operator: '=', value: 1}];
    try {
      const first = await items().updateWhere(pending(), {quantity: 9}, {limit: 2});
      expect(first.updated, 2, `the limit did not cap the call: ${JSON.stringify(first)}`);
      expect(first.hasMore, true, 'a capped call with more matching rows must report hasMore');

      const second = await items().updateWhere(pending(), {quantity: 9}, {limit: 2});
      expect(second.updated, 1, `the loop did not drain the rest: ${JSON.stringify(second)}`);
      expect(second.hasMore, false, 'the drained filter still reports hasMore');
      expect(await items().count([like('sku', prefix), 'and',
        {property: 'quantity', operator: '=', value: 9}] as any), 3, 'the loop left a row behind');
    } finally {
      await cleanup(prefix);
    }
  });

  test('a refused value refuses the whole call: immutable, unknown, service and relation columns', async () => {
    const prefix = `bulk-refuse-${stamp()}`;
    // Two inserts, oldest first: `updateWhere` walks the selection in `created_on, id`
    // order, so the row that does NOT refuse is attempted BEFORE the one that does —
    // a server that wrote before refusing would leave its 'origin' behind.
    const [plain] = await items().insert({sku: `${prefix}-0`, name: 'Plain', quantity: 1});
    const [seeded] = await items().insert({sku: `${prefix}-1`, name: 'Seeded', quantity: 1, origin: 'seed'});
    const filter = like('sku', prefix);
    try {
      expect((await items().query({filter, sort: 'created_on'}))[0].id, plain.id,
        'the fixture rows do not carry the order updateWhere attempts them in');
      const version = (await items().get(seeded.id)).version;
      const immutable = await thrown(() => items().updateWhere(filter, {origin: 'rewritten'}));
      expect(immutable instanceof DG.DomainValidationError, true,
        `a write-once column must refuse the call: ${immutable?.constructor?.name}: ${immutable}`);
      // The refusal is the whole transaction: the row whose 'origin' was still null
      // (and therefore writable) did not get the value either.
      expect((await items().get(plain.id)).origin, null, 'a refused updateWhere wrote a row');

      for (const values of [{created_on: '2020-01-01'}, {'~is_deleted': false}, {tags: []}]) {
        const refusal = await thrown(() => items().updateWhere(filter, values as any));
        expect(refusal instanceof DG.DomainValidationError, true,
          `${Object.keys(values)[0]} must be refused: ${refusal?.constructor?.name}: ${refusal}`);
      }

      const empty = await thrown(() => items().updateWhere('', {quantity: 3}));
      expect(empty instanceof DG.DomainValidationError, true,
        `an empty filter must be refused: ${empty?.constructor?.name}: ${empty}`);

      const untouched = await items().query({filter, sort: 'sku'});
      expect(untouched.map((r: any) => r.quantity).join(','), '1,1', 'a refused call wrote a value');
      expect((await items().get(seeded.id)).version, version, 'a refused call bumped a version');
    } finally {
      await cleanup(prefix);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
