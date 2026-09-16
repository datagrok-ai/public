import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {thrown} from './domain-lifecycle';

// A table that declares `hierarchy` has exactly one ref column pointing at itself, and
// two surfaces over it: `pathTo(id)` — the row's ancestors, root first and the row itself
// excluded (the breadcrumb in front of it) — and the `under` filter term, which matches
// everything in a node's subtree (through the hierarchy table's own id, or through any
// ref column targeting it). Both walk the
// View predicate at every level, so an invisible ancestor truncates rather than leaks.
// The fixture is a throwaway user-managed schema (root › mid › leaf plus two rows
// hanging off it); the category skips cleanly without the CreateDomainSchema privilege.
category('Dapi: domain hierarchy', () => {
  const name = `zzh${`${Date.now()}`.slice(-8)}`;
  const folders = () => grok.dapi.domains.table(`${name}.folder`);
  const things = () => grok.dapi.domains.table(`${name}.thing`);
  const id: {[code: string]: string} = {};
  let skip: string | null = null;

  before(async () => {
    try {
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Hierarchy probe'});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'forbidden' || e.status === 403)) {
        skip = 'no CreateDomainSchema privilege';
        return;
      }
      throw e;
    }
    await grok.dapi.domains.schema(name).apply({tables: {
      folder: {
        hierarchy: true,
        businessKey: ['code'],
        columns: {
          code: {type: 'string', required: true, unique: true},
          name: {type: 'string', isName: true},
          parent_id: {type: 'ref', ref: 'folder'},
        },
      },
      thing: {
        columns: {
          name: {type: 'string', required: true, isName: true},
          folder_id: {type: 'ref', ref: 'folder'},
        },
      },
    }});
    for (const [code, parent] of [['root', null], ['mid', 'root'], ['leaf', 'mid']] as [string, string][]) {
      const [row] = await folders().insert({code, name: code, parent_id: parent == null ? null : id[parent]});
      id[code] = row.id;
    }
    await things().insert([{name: 'in-mid', folder_id: id['mid']}, {name: 'in-leaf', folder_id: id['leaf']}]);
  });

  after(async () => {
    if (skip == null)
      await grok.dapi.domains.schema(name).delete();
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  test('tableInfo says which table is a tree, and on which column', async () => {
    if (skipped())
      return;
    const folder = await grok.dapi.domains.registry.tableInfo(`${name}.folder`);
    expect(folder.hierarchy, true, `folder must be a hierarchy: ${JSON.stringify(folder)}`);
    expect(folder.parentColumn, 'parent_id', `the self-ref column is not named: ${JSON.stringify(folder)}`);
    const thing = await grok.dapi.domains.registry.tableInfo(`${name}.thing`);
    expect(thing.hierarchy === true, false, 'a table with no self-ref must not be a hierarchy');
    expect(thing.parentColumn, null, 'a non-hierarchy table named a parent column');
  });

  test('pathTo: the ancestors of a 3-level leaf, root first and the row itself excluded', async () => {
    if (skipped())
      return;
    const path = await folders().pathTo(id['leaf']);
    expect(path.map((p) => p.name).join(' / '), 'root / mid', `unexpected path: ${JSON.stringify(path)}`);
    expect(path[0].id, id['root'], 'the path does not start at the root');
    expect((await folders().pathTo(id['mid'])).map((p) => p.name).join(' / '), 'root',
      'the middle node has the root as its only ancestor');
    expect((await folders().pathTo(id['root'])).length, 0, 'a root row has no ancestors');
  });

  test('under: the whole subtree of a node, through id and through a child ref', async () => {
    if (skipped())
      return;
    const codes = async (client: any, filter: string): Promise<string> =>
      (await client.query({filter, sort: 'code'})).map((r: any) => r.code).join(',');
    expect(await codes(folders(), `id under "${id['root']}"`), 'leaf,mid,root',
      'under must return the node and everything below it');
    expect(await codes(folders(), `id under "${id['mid']}"`), 'leaf,mid', 'under returned the wrong subtree');
    expect(await codes(folders(), `id under "${id['leaf']}"`), 'leaf', 'a leaf is its own subtree');

    const names = async (filter: string): Promise<string> =>
      (await things().query({filter, sort: 'name'})).map((r: any) => r.name).join(',');
    expect(await names(`folder_id under "${id['mid']}"`), 'in-leaf,in-mid',
      'a ref column into the tree must match the whole subtree');
    expect(await names(`folder_id under "${id['leaf']}"`), 'in-leaf', 'the deepest subtree matched too much');
  });

  test('a table with no hierarchy refuses both surfaces', async () => {
    if (skipped())
      return;
    const [thing] = await things().query({limit: 1});
    const path = await thrown(() => things().pathTo(thing.id));
    expect(path instanceof DG.DomainUnsupportedError && path.op === 'ancestors', true,
      `pathTo on a flat table must be refused as unsupported: ${path?.constructor?.name}: ${path}`);
    // No oracle: `under` on a column that is not a tree reads as an unknown column.
    const filter = await thrown(() => things().query({filter: `name under "${thing.id}"`}));
    expect(filter instanceof DG.DomainFilterError, true,
      `under on a non-ref column must be refused: ${filter?.constructor?.name}: ${filter}`);
  });
}, {owner: 'askalkin@datagrok.ai'});
