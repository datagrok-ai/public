import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// Soft delete is reversible: `deleted: 'exclude' | 'include' | 'only'` decides which
// rows a read sees (anything but 'exclude' projects ~is_deleted), and restore is the
// Delete right spent backwards — audited as 'undelete'. A row whose reference points
// at a still-deleted parent cannot come back alone. Fixture: apitests.item and its
// cascade child apitests.item_event; every test cleans its prefix up.
category('Dapi: domain trash', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const events = () => grok.dapi.domains.table('apitests.item_event');
  const stamp = () => `${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const like = (property: string, prefix: string): any =>
    ({property, operator: 'like', value: `${prefix}%`});

  async function cleanup(prefix: string): Promise<void> {
    try {
      for (const id of (await events().query({filter: like('kind', prefix), deleted: 'only'})).map((r) => r.id))
        await events().restore(id);
      for (const id of (await items().query({filter: like('sku', prefix), deleted: 'only'})).map((r) => r.id))
        await items().restore(id);
      await events().deleteWhere(like('kind', prefix));
      await items().deleteWhere(like('sku', prefix));
    } catch (e) {
      console.error(`trash fixture ${prefix} not cleaned up: ${e}`);
    }
  }

  test('delete → deleted: only → restore: version +2, an undelete audit entry, back in a default query', async () => {
    const prefix = `dt-round-${stamp()}`;
    const [row] = await items().insert({sku: `${prefix}-0`, name: 'Trash round trip'});
    const filter = like('sku', prefix);
    try {
      const version = (await items().get(row.id)).version;
      await items().delete(row.id);
      expect((await items().query({filter})).length, 0, 'a deleted row is still in a default query');
      const [deleted] = await items().query({filter, deleted: 'only'});
      expect(deleted?.id, row.id, 'the deleted row is not in the trash');
      expect((deleted as any)['~is_deleted'], true, `~is_deleted is not set: ${JSON.stringify(deleted)}`);
      expect((await items().query({filter, deleted: 'include'})).length, 1, 'include does not see the deleted row');
      expect(await items().count(filter, {deleted: 'only'}), 1, 'the trash count does not agree with its rows');
      expect(await items().count(filter), 0, 'a default count still counts the deleted row');

      const restored = await items().restore(row.id);
      expect(restored.restored, true, `restore did not report success: ${JSON.stringify(restored)}`);
      expect(restored.version, version + 2, 'the delete and the restore did not bump the version once each');
      expect((await items().get(row.id)).version, restored.version, 'restore reported another version than it wrote');
      const audit = await items().audit(row.id);
      expect(audit.some((a) => a.op === 'undelete'), true,
        `no undelete entry in the trail: ${audit.map((a) => a.op).join(', ')}`);
      expect((await items().query({filter})).length, 1, 'the restored row is not in a default query');
      expect((await items().query({filter, deleted: 'only'})).length, 0, 'the restored row is still in the trash');

      const again = await thrown(() => items().restore(row.id));
      expect(again?.code, 'not-found', `restoring a live row must be a not-found: ${again}`);
    } finally {
      await cleanup(prefix);
    }
  });

  test('a restore under a deleted parent is refused naming the column', async () => {
    const prefix = `dt-parent-${stamp()}`;
    const [item] = await items().insert({sku: `${prefix}-0`, name: 'Parent'});
    const [event] = await events().insert({item_id: item.id, kind: `${prefix}-ev`, amount: 1});
    try {
      await items().delete(item.id);
      const [cascaded] = await events().query({filter: like('kind', prefix), deleted: 'only'});
      expect(cascaded?.id, event.id, 'the cascade did not delete the child row');

      const refusal = await thrown(() => events().restore(event.id));
      expect(refusal?.code, 'restrict', `a child under a deleted parent must be refused: ${refusal}`);
      expect(refusal?.body?.column, 'item_id', `the refusal does not name the column: ${JSON.stringify(refusal?.body)}`);
      expect(`${refusal}`.includes('item_id'), true, `the message does not name the column: ${refusal}`);

      await items().restore(item.id);
      expect((await events().restore(event.id)).restored, true, 'the child did not come back under a live parent');
      expect((await events().query({filter: like('kind', prefix)})).length, 1, 'the restored child is not queryable');
    } finally {
      await cleanup(prefix);
    }
  });

  test('toCsv() of a deleted: include frame carries no service column', async () => {
    const prefix = `dt-csv-${stamp()}`;
    await items().insert({sku: `${prefix}-0`, name: 'Kept'});
    const [gone] = await items().insert({sku: `${prefix}-1`, name: 'Gone'});
    try {
      await items().delete(gone.id);
      const df = await items().queryDf({filter: like('sku', prefix), deleted: 'include', sort: 'sku'});
      expect(df.rowCount, 2, 'include did not return both rows');
      expect(df.col('~is_deleted') != null, true, `~is_deleted is not in the frame: ${df.columns.names().join(', ')}`);
      for (const name of DG.DOMAIN_SERVICE_COLUMNS.filter((c) => df.col(c) != null)) {
        expect(df.col(name)!.meta.includeInCsvExport, false, `${name} is not excluded from csv export`);
        expect(df.col(name)!.meta.includeInBinaryExport, false, `${name} is not excluded from binary export`);
      }
      expect(df.get('~is_deleted', 0), false, 'the live row reads as deleted');
      expect(df.get('~is_deleted', 1), true, 'the deleted row reads as live');
      expect(df.toCsv().includes('~'), false, `a service column leaked into toCsv(): ${df.toCsv()}`);
    } finally {
      await cleanup(prefix);
    }
  });
}, {owner: 'askalkin@datagrok.ai'});
