import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {domains} from '@datagrok-libraries/u2/src/dg/index.js';

// The schema as deployed plus the GHS vocabulary the seed scripts land: opened through the
// same handle the app opens, counted through the table clients. Nothing is written.
category('Stockroom: schema', () => {
  const count = (table: string) => grok.dapi.domains.table(`stockroom.${table}`).count();

  test('substance handle: name column, searchable columns, hazards relation', async () => {
    const substances = await domains.table('stockroom.substance');
    expect(substances.info.nameColumn, 'name');
    expect(substances.info.searchableColumns.includes('cas'), true, 'cas is searchable');
    expect(substances.properties.some((p) => p.name === 'smiles' && p.semType === 'Molecule'), true);
  });

  test('container: constraints and the site-bound location filter', async () => {
    const containers = await domains.table('stockroom.container');
    expect(containers.info.constraints.map((c) => c.name).sort().join(','), 'dates,quantity');
    expect(containers.info.refFilters.location_id, 'site = $site');
  });

  test('seeded GHS vocabulary', async () => {
    expect(await count('hazard_class'), 29);
    expect(await count('h_statement'), 80);
    expect(await count('p_statement'), 97);
  });

  // The demo rows are asserted by their own keys, not by the table totals: on a stand people
  // have used, substances and containers they made by hand sit beside the seeded ones.
  test('seeded stockroom', async () => {
    const seededCas = ['67-64-1', '64-17-5', '67-56-1', '67-63-0', '108-88-3', '110-54-3', '75-09-2',
      '67-66-3', '109-99-9', '141-78-6', '75-05-8', '67-68-5', '64-19-7', '7647-01-0', '7664-93-9',
      '1310-73-2', '7647-14-5', '7722-84-1', '26628-22-8', '1336-21-6'];
    const substances = grok.dapi.domains.table('stockroom.substance');
    expect(await substances.count(`cas in (${seededCas.map((cas) => `"${cas}"`).join(', ')})`), 20);
    expect(await count('location'), 10);
    expect(await count('vendor'), 1);
    expect(await grok.dapi.domains.table('stockroom.container').count('site = "Pilot plant"'), 2);
    expect(await substances.count('cas = "67-64-1"', {search: 'aceto'}), 1);
  });

  test('location is a hierarchy: ancestors and the subtree filter', async () => {
    const locations = await domains.table('stockroom.location');
    expect(locations.info.hierarchy, true, 'location must declare "hierarchy": true');
    expect(locations.info.parentColumn, 'parent_id');

    const client = grok.dapi.domains.table('stockroom.location');
    const shelf = await client.first({filter: 'site = "Main campus" and name = "Flammables cabinet"'});
    const path = await client.pathTo(shelf!.id);
    expect(path.map((p) => p.name).join(' > '), 'Main campus > Building A > Lab 101',
      'pathTo is root-first and excludes the row itself');
    expect(await client.count(`id under "${path[0].id}"`), 6, 'the Main campus subtree');

    // the same term through a ref column INTO the tree — what a tree node click filters by
    const pilot = await client.first({filter: 'site = "Pilot plant" and kind = "site"'});
    expect(await client.count(`id under "${pilot!.id}"`), 4);
    expect(await grok.dapi.domains.table('stockroom.container')
      .count(`location_id under "${pilot!.id}"`), 2, 'both Pilot plant containers, two levels down');
  });
});

// Soft delete on a location of its own: deleted, found in the trash, restored, and gone from
// the live rows in between. The probe row is created and removed by the test.
category('Stockroom: trash', () => {
  test('a deleted location is addressable in the trash and restorable', async () => {
    const client = grok.dapi.domains.table('stockroom.location');
    const name = `Trash probe ${Date.now() % 1e10}`;
    const [probe] = await client.insert({name, kind: 'cabinet', site: 'Main campus'});
    try {
      await client.delete(probe.id);
      expect((await client.query({filter: `id = "${probe.id}"`})).length, 0, 'gone from the live rows');

      const trashed = await client.get(probe.id, {deleted: 'only'});
      expect(trashed['~is_deleted'], true, 'a trash row carries ~is_deleted');
      expect(await client.count(`id = "${probe.id}"`, {deleted: 'only'}), 1, 'the trash count agrees');

      const restored = await client.restore(probe.id);
      expect(restored.restored, true);
      expect((await client.query({filter: `id = "${probe.id}"`})).length, 1, 'back among the live rows');

      const audit = await client.audit(probe.id);
      expect(audit.map((a) => a.op).join(','), 'insert,delete,undelete');
    } finally {
      await client.delete(probe.id);
    }
  });
});
