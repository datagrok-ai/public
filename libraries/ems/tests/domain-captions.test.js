/* A table named by a number — an autonumbered key with a display format — reads in that format
   wherever a row is named: the list's caption, the picker's items, the name a ref resolves to and
   the caption a query projects for a ref (which the server answers bare); the picker finds a row
   by the number a typed id carries, and by its other searchable columns where the text has none. */

import {emsBackends} from '../src/core/domain-backend.js';
import {DomainSource} from '../src/core/domain-source.js';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope} from '@datagrok-libraries/u2';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains, DomainTable} = await import('../src/dg/facade.js');
const {DomainPick} = await import('../src/dg/pick.js');
const {Access, Rows} = await import('@datagrok-libraries/u2');

const SCHEMA = {
  name: 'reg',
  tables: {
    compound: {
      businessKey: ['regno'], friendlyName: 'Compounds',
      columns: {regno: {type: 'int', format: 'DG-000000', isName: true}, name: {type: 'string', searchable: true},
        smiles: {type: 'string'}},
    },
    batch: {
      friendlyName: 'Batches',
      columns: {name: {type: 'string'}, compound_id: {type: 'ref', ref: 'compound'},
        batch_regno: {type: 'int', format: 'DGB-000000'}, potency: {type: 'float'}},
    },
    lot: {columns: {number: {type: 'int', format: 'L-000000', isName: true}}},
    label: {
      businessKey: ['number'], friendlyName: 'Labels',
      columns: {number: {type: 'int'}},
    },
  },
};
const ROWS = {compound: [{id: 'c1', regno: 7, name: 'Aspirin', smiles: 'CCO'},
  {id: 'c2', regno: 123, name: 'Caffeine', smiles: 'CCN'}],
batch: [{id: 'b1', name: 'B1', compound_id: 'c1', batch_regno: 45, potency: 1.2e-9}],
lot: [{id: 't1', number: 5}],
label: [{id: 'l1', number: 1000}, {id: 'l2', number: 2000}]};

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    const dg = globalThis.DG;
    try {
      globalThis.DG = {...dg, format: (x, format) => `${format.replace(/0+$/, '')}${String(x).padStart(6, '0')}`};
      emsBackends.domain = new MemoryDomainBackend(SCHEMA, {rows: structuredClone(ROWS)});
      await body();
    } finally {
      globalThis.DG = dg;
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

scoped('the caption of a row named by a formatted number is the formatted number', async () => {
  const table = await domains.table('reg.compound');
  const source = table.source({});
  await flush();
  assert.equal(table.renderer.caption(source.rows.byKey('c2')), 'DG-000123');
  source.dispose();
});

scoped('the picker searches a numeric name by its digits and names items in the format', async () => {
  assert.deepEqual((await DomainPick.search('reg.compound', 'DG-000123')).map((i) => [i.id, i.name]),
    [['c2', 'DG-000123']]);
  assert.deepEqual((await DomainPick.search('reg.compound', '7')).map((i) => i.name), ['DG-000007']);
  assert.equal((await DomainPick.resolve('reg.compound', 'c1')).name, 'DG-000007');
  assert.equal(DomainPick.digits('DG-000123'), '000123');
  assert.equal(DomainPick.digits('DG-'), '');
});

scoped('captionText leaves a column without a format as its text', async () => {
  const table = await domains.table('reg.compound');
  assert.equal(DomainTable.captionText(table.properties, 'smiles', 'CCO'), 'CCO');
});

scoped('a number without a display format is named and searched as before', async () => {
  const table = await domains.table('reg.label');
  assert.equal(DomainSource.formattedNumber(table.properties, 'number') !== undefined, false);
  assert.deepEqual((await DomainPick.search('reg.label', '')).map((i) => i.name), ['1000', '2000']);
});

scoped('a projected caption of a numbered target reads in the target\'s format — on a form and in a data table',
  async () => {
    const table = await domains.table('reg.batch');
    const source = table.source({pageSize: 10});
    await flush();
    const row = source.rows.byKey('b1');
    const [answered] = await emsBackends.domain.tableSync('reg.batch').query({captions: ['compound_id']});
    assert.equal(answered[Rows.caption('compound_id')], '7', 'the query answers the name bare');
    assert.equal(row[Rows.caption('compound_id')], 'DG-000007', 'the source formats it where the frame lands');
    source.currentRow.value = row;
    const form = domains.form(source, {access: Access.from({can: {view: true}, fields: {compound_id: 'readonly'}})});
    const shown = () =>
      form.root.querySelector('[data-u2-name="compound_id"] [data-u2-part="readonly-value"]').textContent;
    await flush();
    assert.equal(shown(), 'DG-000007');
    form.dispose();

    const dt = domains.dataTable(source, {columns: ['compound_id', 'batch_regno', 'potency']});
    document.body.append(dt.root);
    await flush();
    dt.table.value.root.clientHeight = 400;
    await flush();
    const cells = [...dt.root.querySelector('.u2-data-table-row[data-index="0"]').children].map((c) => c.textContent);
    assert.deepEqual(cells, ['DG-000007', 'DGB-000045', '1.2e-9'],
      'the target\'s format, an int column\'s own, and an unformatted float never rounded to 0');
    dt.dispose();
    source.dispose();
  });

scoped('a form opened before the source has formatted its captions still reads the target\'s format', async () => {
  const formats = DomainSource.nameFormat;
  let release;
  DomainSource.nameFormat = (target) => new Promise((resolve) => release = () => resolve(formats(target)));
  try {
    const table = await domains.table('reg.batch');
    const source = table.source({pageSize: 10});
    await flush();
    source.currentRow.value = source.rows.byKey('b1');
    const form = domains.form(source, {access: Access.from({can: {view: true}, fields: {compound_id: 'readonly'}})});
    await flush();
    assert.equal(form.root.querySelector('[data-u2-name="compound_id"] [data-u2-part="readonly-value"]').textContent,
      'DG-000007');
    release();
    form.dispose();
    source.dispose();
  } finally {
    DomainSource.nameFormat = formats;
  }
});

scoped('nameFormat: a numbered target\'s format for a bare name; none where the name is not a formatted number',
  async () => {
    const format = await DomainPick.nameFormat('reg.compound');
    assert.equal(format('7'), 'DG-000007');
    assert.equal(format('B-7'), 'B-7', 'a name that is not bare digits stands');
    assert.equal(await DomainPick.nameFormat('reg.label'), null);
    assert.equal(await DomainPick.nameFormat('reg.nonesuch'), null, 'a table that cannot be read formats nothing');
  });

scoped('text without digits searches a numbered table\'s other searchable columns; with none there, nothing', async () => {
  assert.deepEqual((await DomainPick.search('reg.compound', 'aspirin')).map((i) => [i.id, i.name]), [['c1', 'DG-000007']]);
  assert.deepEqual(await DomainPick.search('reg.lot', 'abc'), [], 'never the first rows unfiltered');
  assert.deepEqual((await DomainPick.search('reg.lot', 'L-5')).map((i) => i.id), ['t1']);
  assert.deepEqual(DomainPick.filter('regno', 'DG-000045', undefined, true),
    [{property: 'regno', operator: '=', value: 45}]);
  assert.equal(DomainPick.filter('regno', 'aspirin', undefined, true), undefined);
  for (const propertyType of ['int', 'bigint', 'double', 'float', 'qnum'])
    assert.equal(DomainSource.formattedNumber([{name: 'x', propertyType, format: '0.00'}], 'x') !== undefined, true, propertyType);
  assert.equal(DomainSource.formattedNumber([{name: 'x', propertyType: 'string', format: 'x'}], 'x') !== undefined, false);
});

scoped('lookUp: names to ids by the name column — a formatted number by its digits — a hundred at a time', async () => {
  const ids = await DomainPick.lookUp('reg.compound', ['DG-000007', '123', 'aspirin', 'DG-999999']);
  assert.deepEqual([...ids], [['DG-000007', 'c1'], ['123', 'c2']], 'a name no row holds is absent');
  const compound = emsBackends.domain.tableSync('reg.compound');
  const query = compound.query.bind(compound);
  const asked = [];
  compound.query = (spec) => {
    asked.push(spec.filter.value.length);
    return query(spec);
  };
  const many = await DomainPick.lookUp('reg.compound', Array.from({length: 250}, (_, i) => `DG-${i + 1}`));
  assert.deepEqual(asked, [100, 100, 50]);
  assert.deepEqual([...many], [['DG-7', 'c1'], ['DG-123', 'c2']]);
  assert.deepEqual([...await DomainPick.lookUp('reg.label', ['1000', 'x', '20'])], [['1000', 'l1']],
    'an unformatted number by its value');
});
