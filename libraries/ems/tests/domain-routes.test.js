/* `domains.route` (WO 3-7) over the memory backend: what a `/domains/<schema>/<table>` address
   resolves to — the u2 app as a view, the row a segment names (a business key, a composite one,
   an id), the query and the search forwarded, an app already open taking the address instead, a
   not-found view for an unknown table, and null wherever the platform keeps the address. The fixture
   here is a table nothing refers to: the route builds the app with its panes, and the children
   pane's platform grid has no stub to build over. */

import {emsBackends} from '../src/core/domain-backend.js';
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/facade.js');
const {DomainApp} = await import('../src/dg/app.js');
const {DomainAddress} = await import('../src/dg/address.js');
const grok = await import('datagrok-api/grok');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
    assert.equal(DomainApp.live.size, 0, 'no app left behind');
  });
}

const SCHEMA = {
  name: 'demo',
  tables: {
    widget: {
      businessKey: ['code'], friendlyName: 'Widgets',
      columns: {
        code: {type: 'string', required: true},
        name: {type: 'string', isName: true, searchable: true},
      },
    },
  },
};

const UUID = '11111111-2222-3333-4444-555555555555';
const widgets = () => new MemoryDomainBackend(SCHEMA, {rows: {widget: [
  {id: UUID, code: 'A-1', name: 'Alpha'}, {id: 'w2', code: 'B2', name: 'Beta'}]}});

scoped('the table route opens the app as a view, under the /domains base', async () => {
  emsBackends.domain = widgets();
  const view = await domains.route('/domains/demo/widget');
  document.body.append(view.root);
  const a = DomainApp.of(view);
  assert.equal(view.name, 'Widgets');
  assert.equal(a.base, '/domains/demo/widget', 'the route is the base, not a rebased app route');
  assert.equal(a.entityPath, true);
  assert.equal(a.page.value, 'list');
  await flush();
  assert.equal(a.listSource.rows.items.value.length, 2);
  assert.equal(view.path, '/domains/demo/widget');
  a.dispose();
});

scoped('the platform route mounts a live app — one table, two addresses, one behaviour', async () => {
  emsBackends.domain = widgets();
  const view = await domains.route('/domains/demo/widget');
  document.body.append(view.root);
  const a = DomainApp.of(view);
  await flush();
  assert.equal(a.listSource.live.value, true, 'the canonical address polls like a package app does');
  a.dispose();
});

scoped('a row segment opens its entity page: a business key with a dash in it, and an id', async () => {
  emsBackends.domain = widgets();
  const view = await domains.route('/domains/demo/widget/A-1');
  document.body.append(view.root);
  const a = DomainApp.of(view);
  assert.equal(a.page.value, 'entity');
  await flush();
  assert.equal(a.form.value.input('name').value.value, 'Alpha', 'a single-column key keeps its dashes');
  assert.equal(a.path.value, '/domains/demo/widget/A-1');
  assert.equal(await domains.route(`/domains/demo/widget/${UUID}`) === view, true,
    'the open app answers for the id too, as itself');
  await flush();
  assert.equal(a.form.value.input('name').value.value, 'Alpha');
  assert.equal(a.path.value, '/domains/demo/widget/A-1', 'the path settles on the key the row carries');
  a.dispose();
});

scoped('a composite key travels as one segment, split by the key\'s arity', async () => {
  emsBackends.domain = backend();
  const view = await domains.route('/domains/grit/issue/p1-2');
  document.body.append(view.root);
  const a = DomainApp.of(view);
  await flush();
  assert.equal(a.form.value.input('title').value.value, 'Ibuprofen', 'the int component is read as an int');
  assert.equal(a.path.value, '/domains/grit/issue/p1-2');
  a.dispose();
});

scoped('the query and the search travel with the address', async () => {
  emsBackends.domain = widgets();
  const view = await domains.route(`/domains/demo/widget?q=${encodeURIComponent('code = "B2"')}&search=bet`);
  document.body.append(view.root);
  const a = DomainApp.of(view);
  await flush();
  assert.equal(a.listSource.query.value, 'code = "B2"');
  assert.equal(a.listSource.search.value, 'bet');
  assert.equal(a.listSource.rows.items.value.length, 1);
  assert.equal(a.path.value, '/domains/demo/widget?q=code%20%3D%20%22B2%22&search=bet');
  a.dispose();
});

scoped('find-or-activate: an app already open takes the address and is answered as itself', async () => {
  emsBackends.domain = widgets();
  const table = await domains.table('demo.widget');
  const open = table.app({path: '/domains/demo/widget'});
  document.body.append(open.root);
  const a = DomainApp.of(open);
  await flush();
  // THE view, never null: the platform reads null as "no u2 route" and opens the Dart domain
  // view over the app that is already showing the address (Back off `?trash=1` did exactly that)
  const pushes = [];
  const saved = globalThis.history;
  globalThis.history = {pushState: (_state, _title, url) => pushes.push(url)};
  try {
    assert.equal(await domains.route('/domains/demo/widget/B2') === open, true, 'the open app answers as itself');
    await flush();
    assert.equal(a.page.value, 'entity');
    assert.equal(a.form.value.input('name').value.value, 'Beta');
    assert.equal(grok.shell.v, open, 'and its view comes to the front');
    const back = await domains.route(`/domains/demo/widget?q=${encodeURIComponent('code = "B2"')}`);
    assert.equal(back === open, true);
    await flush();
    assert.equal(a.page.value, 'list', 'the whole address, not just the row');
    assert.equal(a.listSource.query.value, 'code = "B2"');
    // the address came FROM the history: restoring it must not push another entry onto it
    assert.deepEqual(pushes, [], 'a route into an open app costs no history entry');
  } finally {
    globalThis.history = saved;
  }
  grok.shell.v = undefined;
  a.dispose();
});

scoped('null where the platform keeps the address: anything that is not a table route',
  async () => {
    emsBackends.domain = widgets();
    assert.equal(await domains.route('/domains'), null);
    assert.equal(await domains.route('/domains/demo'), null, 'the schema gallery stays Dart');
    assert.equal(await domains.route('/domains/demo/widget/A-1/extra'), null);
    assert.equal(await domains.route('/apps/Grit'), null);
  });

scoped('a table the registry does not know opens a not-found view; any other failed read stays the platform\'s',
  async () => {
    const refusing = (e) => ({table: () => Promise.reject(e)});
    for (const [be, address, name] of [[widgets(), '/domains/demo/nope/A-1?x=1', 'demo.nope'],
      [refusing({code: 'validation', message: 'Unknown domain table "nope.widget"'}), '/domains/nope/widget', 'nope.widget'],
      [refusing({code: 'forbidden'}), '/domains/demo/widget', null],
      [refusing(new Error('Failed to fetch')), '/domains/demo/widget', null]]) {
      emsBackends.domain = be;
      const view = await domains.route(address);
      assert.equal(view?.root.textContent ?? null, name && `No table matches "${name}". It may have been deleted, ` +
        'renamed, or you may not have access to it.', address);
      if (view === null)
        continue;
      assert.equal(view.path, address, 'the address stays as typed');
      assert.deepEqual(view.ribbonPanels, [], 'no New, no import');
      globalThis.grok_Widget_Kill(view.root);
    }
  });

test('DomainAddress: the platform\'s row route, and what an address carries under an app\'s bases', () => {
  assert.deepEqual(DomainAddress.ROUTE.exec('/domains/demo/widget').slice(1),
    ['demo', 'widget', undefined]);
  assert.deepEqual(DomainAddress.ROUTE.exec('/domains/demo/widget/A-1').slice(1), ['demo', 'widget', '/A-1']);
  assert.equal(DomainAddress.ROUTE.exec('/domains/demo'), null, 'the schema gallery is not a table route');
  assert.equal(DomainAddress.ROUTE.exec('/domains/demo/widget/A-1/extra'), null);
  assert.equal(DomainAddress.entityPath('/domains/demo/widget'), true);
  assert.equal(DomainAddress.entityPath('/apps/Grit'), false, 'an app addresses rows by ?entity=');

  // both routes at once: the one the shell mounted the view at, and the one it was built with
  const bases = ['/apps/Grit', '/domains/grit/issue'];
  assert.equal(DomainAddress.restOf('/apps/Grit', bases), '', 'the base itself');
  assert.equal(DomainAddress.restOf('/APPS/grit/p1-1', bases), '/p1-1', 'the case of an address is not its own');
  assert.equal(DomainAddress.restOf('/domains/grit/issue/p1-1', bases), '/p1-1', 'the fallback route too');
  assert.equal(DomainAddress.restOf('/apps/GritLabs', bases), null, 'a segment boundary, not a prefix');
  assert.equal(DomainAddress.restOf('/apps/Other/x', bases), null);
  assert.equal(DomainAddress.segmentOf('/apps/Grit/p1-1', bases), 'p1-1');
  assert.equal(DomainAddress.segmentOf('/apps/Grit/p1%2D1/extra', bases), 'p1-1', 'the first segment names the row');
  assert.equal(DomainAddress.segmentOf('/apps/Grit', bases), null, 'a base names no row');
  assert.equal(DomainAddress.segmentOf('/apps/Other/x', bases), null);
});

test('DomainAddress: a key equal to a route word, a dot segment, "new" or led by ~ round-trips a segment escaped', () => {
  // RESERVED is checked against a literal copy of the Dart DomainRowKey.reserved list (the Datlas suite
  // domain_row_key_test pins that list to the router's route table)
  assert.deepEqual(DomainAddress.RESERVED, ['access', 'aggregate', 'audit', 'batch', 'delete', 'facets', 'filters',
    'query', 'update', 'version', 'watch', 'new', '.', '..']);
  const bases = ['/domains/grit/issue'];
  for (const key of [...DomainAddress.RESERVED, '~x', '~~', '~', 'a/b', 'a%2Fb', '1,a%2Cb', 'p.q', 'GRIT', UUID]) {
    const segment = DomainAddress.segment(key);
    assert.equal(DomainAddress.unescape(DomainAddress.segmentOf(`/domains/grit/issue/${segment}`, bases)), key,
      `${key} → ${segment} → back`);
  }
  assert.equal(DomainAddress.segment('watch'), '~watch');
  assert.equal(DomainAddress.segment('~watch'), '~~watch');
  assert.equal(DomainAddress.segment('..'), '~..', 'a browser folds a bare dot segment away');
  assert.equal(DomainAddress.segment('a%2Fb'), 'a%252Fb', 'the canonical id is percent-encoded once more');
  assert.equal(DomainAddress.segment(UUID), UUID);
  assert.equal(DomainAddress.segment('Version'), 'Version', 'exact and case-sensitive');
  assert.equal(DomainAddress.segmentOf('/domains/grit/issue/~new', bases), '~new', 'segmentOf keeps the escape: the draft sentinel is read off the raw segment');
  assert.equal(DomainAddress.unescape('new'), 'new');
});

test('DomainAddress: how a row is spelled in a path, and read back out of one', () => {
  assert.equal(DomainAddress.keyOf({id: 'x', code: 'A-1'}, ['code']), 'A-1',
    'a single-column key keeps its dashes');
  assert.equal(DomainAddress.keyOf({id: 'x', project_id: 'p1', number: 2}, ['project_id', 'number']), 'p1-2');
  assert.equal(DomainAddress.keyOf({id: 'x', project_id: 'p-1', number: 2}, ['project_id', 'number']), 'x',
    'the split would be ambiguous');
  assert.equal(DomainAddress.keyOf({id: 'x', code: null}, ['code']), 'x', 'a null component is not addressable');
  assert.equal(DomainAddress.keyOf({id: 'x'}, []), 'x', 'no business key, no key address');

  const schema = {properties: [{name: 'project_id', type: 'string'}, {name: 'number', type: 'int'}]};
  const query = DomainAddress.keyQuery('p1-2', ['project_id', 'number'], schema);
  assert.deepEqual(query.nodes.map((n) => [n.property, n.operator, n.value]),
    [['project_id', '=', 'p1'], ['number', '=', 2]], 'the int component is read as an int');
  assert.equal(DomainAddress.keyQuery(UUID, ['code'], {properties: []}), null, 'a uuid IS the id');
  assert.equal(DomainAddress.keyQuery('A-1', [], {properties: []}), null);
  assert.equal(DomainAddress.keyQuery('p1-1-2', ['project_id', 'number'], {properties: []}), null,
    'an arity that does not match is read as an id');
});
