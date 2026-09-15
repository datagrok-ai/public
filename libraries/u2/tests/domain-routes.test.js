/* `domains.route` (WO 3-7) over the memory backend: what a `/domains/<schema>/<table>` address
   resolves to — the u2 app as a view, the row a segment names (a business key, a composite one,
   an id), the query and the search forwarded, an app already open taking the address instead, and
   null wherever the platform keeps the address (a schema route, an unknown table). The fixture
   here is a table nothing refers to: the route builds the app with its panes, and the children
   pane's platform grid has no stub to build over. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {MemoryDomainBackend} from '../src/sources/memory-domain.js';
import {notify} from '../src/components/display/notify.js';
import {backend} from './domain-fixtures.mjs';

register('./dg-stub.mjs', import.meta.url);
const {domains} = await import('../src/dg/domain/index.js');
const {DomainApp} = await import('../src/dg/domain/app.js');
const grok = await import('datagrok-api/grok');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
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
  backends.domain = widgets();
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

scoped('a row segment opens its entity page: a business key with a dash in it, and an id', async () => {
  backends.domain = widgets();
  const view = await domains.route('/domains/demo/widget/A-1');
  document.body.append(view.root);
  const a = DomainApp.of(view);
  assert.equal(a.page.value, 'entity');
  await flush();
  assert.equal(a.form.value.input('name').value.value, 'Alpha', 'a single-column key keeps its dashes');
  assert.equal(a.path.value, '/domains/demo/widget/A-1');
  assert.equal(await domains.route(`/domains/demo/widget/${UUID}`), null, 'the open app answers for the id too');
  await flush();
  assert.equal(a.form.value.input('name').value.value, 'Alpha');
  assert.equal(a.path.value, '/domains/demo/widget/A-1', 'the path settles on the key the row carries');
  a.dispose();
});

scoped('a composite key travels as one segment, split by the key\'s arity', async () => {
  backends.domain = backend();
  const view = await domains.route('/domains/grit/issue/p1-2');
  document.body.append(view.root);
  const a = DomainApp.of(view);
  await flush();
  assert.equal(a.form.value.input('title').value.value, 'Ibuprofen', 'the int component is read as an int');
  assert.equal(a.path.value, '/domains/grit/issue/p1-2');
  a.dispose();
});

scoped('the query and the search travel with the address', async () => {
  backends.domain = widgets();
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

scoped('find-or-activate: an app already open over the table takes the address and nothing is docked', async () => {
  backends.domain = widgets();
  const table = await domains.table('demo.widget');
  const open = table.app({path: '/domains/demo/widget'});
  document.body.append(open.root);
  const a = DomainApp.of(open);
  await flush();
  assert.equal(await domains.route('/domains/demo/widget/B2'), null, 'the open app answers for it');
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.form.value.input('name').value.value, 'Beta');
  assert.equal(grok.shell.v, open, 'and its view comes to the front');
  assert.equal(await domains.route(`/domains/demo/widget?q=${encodeURIComponent('code = "B2"')}`), null);
  await flush();
  assert.equal(a.page.value, 'list', 'the whole address, not just the row');
  assert.equal(a.listSource.query.value, 'code = "B2"');
  grok.shell.v = undefined;
  a.dispose();
});

scoped('null where the platform keeps the address: a schema route, a table the registry does not know',
  async () => {
    backends.domain = widgets();
    assert.equal(await domains.route('/domains'), null);
    assert.equal(await domains.route('/domains/demo'), null, 'the schema gallery stays Dart');
    assert.equal(await domains.route('/domains/demo/widget/A-1/extra'), null);
    assert.equal(await domains.route('/apps/Grit'), null);
    assert.equal(await domains.route('/domains/demo/nope'), null, 'an unknown table falls back');
    assert.equal(await domains.route('/domains/nope/widget'), null);
  });
