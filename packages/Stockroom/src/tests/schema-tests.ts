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

  test('seeded stockroom', async () => {
    expect(await count('substance'), 20);
    expect(await count('location'), 8);
    expect(await count('vendor'), 1);
    expect(await count('container'), 7);
    expect(await grok.dapi.domains.table('stockroom.substance').count('cas = "67-64-1"', {search: 'aceto'}), 1);
  });
});
