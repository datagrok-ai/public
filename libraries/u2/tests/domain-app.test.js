/* `domains.app` (WO 2-7) over the memory backend: the list and entity pages, `path` for the list,
   a query and an entity, `open()` from a path, the gate in front of a move, breadcrumbs, the
   summary following the page, the entity source sharing the session (one transaction for an
   edit on each page), find-or-activate through `DomainTable.open`, and New as a draft page. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {SharedSession} from '../src/sources/session.js';
import {Rows} from '../src/sources/rows-like.js';
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

const BASE = '/apps/T/Issues';
const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);
const crumbs = (app) => [...app.breadcrumbs.root.querySelectorAll('.u2-breadcrumbs-item, .u2-breadcrumbs-current')]
  .filter((el) => el.style.display !== 'none' && el.textContent !== '…').map((el) => el.textContent);

async function app(options = {}) {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const a = domains.app({table, base: BASE, pageSize: 10, ...options});
  document.body.append(a.root);
  await flush();
  return {table, app: a};
}

scoped('the list page first; path follows the query; goTo the entity page and back', async () => {
  const {app: a} = await app();
  assert.equal(a.root.dataset.u2, 'domain-app');
  assert.equal(a.page.value, 'list');
  assert.equal(a.entity.value, null);
  assert.equal(a.path.value, BASE);
  assert.equal(a.listSource.rows.items.value.length, 3);
  assert.equal(a.summary.value, '3 issues');
  a.listSource.query.value = 'done = true';
  assert.equal(a.path.value, `${BASE}?q=done%20%3D%20true`);
  assert.equal(await a.goTo('entity', 'i2'), true);
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, 'i2');
  assert.equal(a.path.value, `${BASE}?entity=i2`);
  const source = a.entitySource.value;
  assert.equal(source.session, a.session, 'the entity source shares the session');
  await flush();
  assert.equal(a.form.value.input('title').value.value, 'Ibuprofen');
  assert.deepEqual(crumbs(a), ['Issues', 'Ibuprofen']);
  assert.equal(a.summary.value, 'Ibuprofen', 'the entity page says which row, not a count');
  assert.equal(a.root.querySelector('[data-u2-part="list-page"]').hidden, true);
  assert.equal(a.root.querySelector('[data-u2-part="entity-page"]').hidden, false);
  fire(a.breadcrumbs.root.querySelector('.u2-breadcrumbs-item'), 'click');
  await flush();
  assert.equal(a.page.value, 'list');
  assert.equal(a.entitySource.value, null, 'the entity source is released');
  assert.equal(a.form.value, null);
  assert.equal(a.path.value, `${BASE}?q=done%20%3D%20true`, 'the list keeps its query');
  a.dispose();
});

scoped('open(path): ?entity= restores the entity page, ?q= the list under the query, nothing the list', async () => {
  const {app: a} = await app();
  assert.equal(await a.open(`${BASE}?entity=i1`), true);
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, 'i1');
  assert.equal(await a.open(`${BASE}?q=${encodeURIComponent('title like "a"')}`), true);
  assert.equal(a.page.value, 'list');
  assert.equal(a.listSource.query.value, 'title like "a"');
  assert.equal(a.path.value, `${BASE}?q=title%20like%20%22a%22`);
  await flush();
  assert.equal(a.listSource.total.value, 2);
  await a.open(BASE);
  assert.equal(a.listSource.query.value, '');
  await a.open(`${BASE}?entity=nope`);
  await flush();
  assert.match(a.form.value.root.querySelector('[data-u2-part="empty"]').textContent, /not found/);
  a.dispose();
});

scoped('a path is authoritative: open(\'\') is the list even while the address bar holds ?entity=', async () => {
  const {app: a} = await app();
  assert.equal(await a.open('?entity=new'), true);
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, DomainApp.NEW);
  location.search = '?entity=new';
  assert.equal(await a.open(''), true);
  assert.equal(a.page.value, 'list');
  assert.equal(a.entity.value, null);
  a.dispose();
});

scoped('the gate: a move with unsaved changes asks; cancel keeps the page and the changes', async () => {
  const {app: a} = await app();
  await a.goTo('entity', 'i1');
  await flush();
  a.form.value.input('title').value.value = 'Aspirin 100';
  await flush();
  assert.equal(a.session.isDirty.value, true);
  assert.equal(a.summary.value, '1 unsaved change', 'the session\'s summary while dirty');
  const answer = a.goTo('list');
  await flush();
  assert.notEqual(document.body.querySelector('.u2-dialog'), null);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await answer, false);
  assert.equal(a.page.value, 'entity');
  assert.equal(a.session.isDirty.value, true);
  assert.equal(a.form.value.input('title').value.value, 'Aspirin 100');
  assert.equal(await a.goTo('entity', 'i1'), true, 'staying put asks nothing');
  const back = a.goTo('list');
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  assert.equal(await back, true);
  assert.equal(a.page.value, 'list');
  assert.equal(a.session.isDirty.value, false);
  a.dispose();
});

scoped('one session: an edit on the entity page and one on the list save as one transaction', async () => {
  const {app: a} = await app();
  let transactions = 0;
  const real = backends.domain.transaction.bind(backends.domain);
  backends.domain.transaction = (ops) => {
    transactions++;
    return real(ops);
  };
  await a.goTo('entity', 'i1');
  await flush();
  a.form.value.input('title').value.value = 'Aspirin 100';
  a.listSource.rows.byKey('i3').title = 'Naproxen 500';
  await flush();
  assert.equal(a.session.changeCount.value, 2);
  assert.equal(a.summary.value, '2 unsaved changes', 'the entity source and the list source are one table');
  assert.equal(await a.session.save(), true);
  assert.equal(transactions, 1);
  const rows = backends.domain.tableSync('grit.issue').rows;
  assert.deepEqual(rows.map((r) => r.title), ['Aspirin 100', 'Ibuprofen', 'Naproxen 500']);
  a.dispose();
});

scoped('find-or-activate: DomainTable.open lands on the app\'s entity page; a registered view comes to the front', async () => {
  const {table, app: a} = await app();
  assert.equal(DomainApp.baseOf('grit.issue'), BASE);
  const other = table.source();
  await flush();
  table.open(other.rows.byKey('i3'));
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, 'i3');
  const view = {name: 'Issues'};
  const unregister = DomainApp.register(a, view);
  assert.equal(DomainApp.of(view), a);
  assert.equal(DomainApp.activate(BASE, 'i2'), true);
  assert.equal(grok.shell.v, view, 'the view is made current');
  await flush();
  assert.equal(a.entity.value, 'i2');
  assert.equal(DomainApp.activate('/apps/Other', 'i2'), false);
  unregister();
  assert.equal(DomainApp.of(view), undefined);
  other.dispose();
  a.dispose();
  assert.equal(DomainApp.baseOf('grit.issue'), undefined);
});

scoped('New: the entity page over a pristine draft; once saved the page is the row\'s', async () => {
  const {app: a} = await app();
  const ribbon = a.ribbon();
  assert.equal(ribbon.length, 2);
  assert.deepEqual(ribbon[0].map((c) => c.root.dataset.u2), ['new-button', 'save-button', 'discard-button']);
  assert.deepEqual(ribbon[1].map((c) => c.root.dataset.u2), ['domain-search', 'domain-filters']);
  assert.equal(ribbon[1][0].root.hidden, false);
  fire(ribbon[0][0].root, 'click');
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, DomainApp.NEW);
  assert.equal(a.path.value, `${BASE}?entity=new`);
  assert.equal(ribbon[1][0].root.hidden, true, 'no search on the entity page');
  assert.equal(ribbon[1][1].root.hidden, true);
  const source = a.entitySource.value;
  assert.equal(source.isDraft, true);
  assert.equal(Rows.isDraft(source.currentRow.value), true);
  assert.equal(a.session.isDirty.value, false, 'pristine');
  assert.deepEqual(crumbs(a), ['Issues', 'New issue']);
  const form = a.form.value;
  form.input('project_id').value.value = 'p1';
  form.input('title').value.value = 'Fresh';
  await flush();
  assert.equal(await a.session.save(), true);
  await flush();
  const saved = backends.domain.tableSync('grit.issue').rows.find((r) => r.title === 'Fresh');
  assert.equal(a.entity.value, saved.id);
  assert.equal(a.path.value, `${BASE}?entity=${saved.id}`);
  assert.equal(a.session.isDirty.value, false);
  a.dispose();
});

scoped('New: discarding the draft is the way back to the list, entity and all', async () => {
  const {app: a} = await app();
  const ribbon = a.ribbon();
  fire(ribbon[0][0].root, 'click');
  await flush();
  a.form.value.input('title').value.value = 'Fresh';
  await flush();
  assert.equal(a.session.isDirty.value, true);
  fire(ribbon[0][2].button, 'click');
  await flush();
  assert.equal(a.page.value, 'list');
  assert.equal(a.entity.value, null);
  assert.equal(a.entitySource.value, null, 'the draft page is released');
  assert.equal(a.path.value, BASE, 'no ?entity= left in the URL');
  assert.equal(a.summary.value, '3 issues', 'the status bar is the list\'s again');
  a.dispose();
});

scoped('New: a refused save then Discard lands on the list too', async () => {
  const {app: a} = await app();
  const ribbon = a.ribbon();
  fire(ribbon[0][0].root, 'click');
  await flush();
  // no project: the form's own gate refuses the batch, and the draft stays pending
  a.form.value.input('title').value.value = 'Fresh';
  await flush();
  assert.equal(await a.session.save(), false);
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, DomainApp.NEW);
  fire(ribbon[0][2].button, 'click');
  await flush();
  assert.equal(a.page.value, 'list');
  assert.equal(a.entity.value, null);
  assert.equal(a.path.value, BASE);
  a.dispose();
});

scoped('the entity page carries the children and history panes under the form; both can be left out', async () => {
  const {app: a} = await app();
  const panes = () => [...a.panes.children].map((el) => el.dataset.u2);
  assert.deepEqual(panes(), [], 'nothing on the list page');
  await a.goTo('entity', 'i1');
  await flush();
  assert.deepEqual(panes(), ['domain-children', 'domain-history']);
  const history = a.panes.querySelector('[data-u2="domain-history"]');
  assert.equal(history.querySelector('[data-u2-part="hint"]').hidden, true, 'a saved row has a history');
  await a.goTo('list');
  assert.deepEqual(panes(), [], 'released with the entity page');
  a.dispose();

  const {app: bare} = await app({children: false, history: false});
  await bare.goTo('entity', 'i1');
  await flush();
  assert.deepEqual([...bare.panes.children], []);
  bare.dispose();
});

scoped('DomainTable.app(): the view over the app — name, path, the registries, and the close gate', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const view = table.app({path: BASE});
  document.body.append(view.root);
  const a = DomainApp.of(view);
  assert.notEqual(a, undefined);
  assert.equal(view.name, 'Issues', 'the display name of the table');
  assert.equal(DomainApp.titleOf({pluralName: 'order_lines'}), 'Order lines',
    'a plural the registry derived from the table name is made presentable');
  assert.equal(DomainApp.titleOf({pluralName: 'order_lines', friendlyName: 'Order Lines'}), 'Order Lines',
    'what the schema declares wins');
  assert.equal(a.base, BASE);
  assert.equal(DomainApp.baseOf('grit.issue'), BASE);
  assert.equal(view.ribbonPanels.length, 2);
  assert.equal(view.statusBarPanels.length, 1);
  await flush();
  assert.equal(view.path, BASE);
  assert.equal(view.acceptsPath('/apps/t/issues'), true);
  assert.equal(view.acceptsPath('/apps/Other'), false);
  const address = globalThis.location;
  globalThis.location = {search: '?entity=i2'};
  try {
    view.handlePath('/apps/T/Issues');
    await flush();
  } finally {
    globalThis.location = address;
  }
  assert.equal(a.entity.value, 'i2', 'a deep link opens the entity page');
  assert.equal(view.path, `${BASE}?entity=i2`, 'the address bar follows');
  assert.equal(view.statusBarPanels[0].textContent, 'Ibuprofen');
  assert.equal(await a.goTo('list'), true);
  await flush();
  assert.equal(view.path, BASE, 'back to the list path');
  a.listSource.query.value = 'done = true';
  await flush();
  assert.equal(view.path, `${BASE}?q=done%20%3D%20true`, 'the query follows too');
  assert.equal(await a.goTo('entity', 'i1'), true);
  await flush();
  assert.equal(view.path, `${BASE}?entity=i1`, 'every page switch is mirrored, not just the first');
  assert.equal(await a.goTo('list'), true);
  await flush();
  assert.equal(view.path, `${BASE}?q=done%20%3D%20true`, 'and back, under the list\'s query');
  await a.goTo('entity', 'i2');
  await flush();

  const removing = grok.events.onViewRemoving;
  const closeRequest = () => ({args: {view: {dart: view.dart}}, prevented: 0, preventDefault() {
    this.prevented++;
  }});
  let e = closeRequest();
  removing.fire(e);
  assert.equal(e.prevented, 0, 'a clean view closes without a word');
  a.form.value.input('title').value.value = 'Ibuprofen 400';
  await flush();
  e = closeRequest();
  removing.fire(e);
  assert.equal(e.prevented, 1, 'unsaved changes: the removal is cancelled');
  await flush();
  fire(buttonNamed('CANCEL'), 'click');
  await flush();
  assert.equal(view.dart.closed, 0, 'and the view stays');
  assert.equal(a.session.isDirty.value, true);
  e = closeRequest();
  removing.fire(e);
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(view.dart.closed, 1, 'decided: closed for real');
  assert.equal(a.session.isDirty.value, false);
  const other = {args: {view: {dart: {}}}, preventDefault: () => assert.fail('another view is not ours')};
  a.form.value.input('title').value.value = 'x';
  removing.fire(other);

  a.dispose();
  assert.equal(DomainApp.of(view), undefined);
  assert.equal(removing.count, 0, 'the gate left with the app');

  const defaults = table.app();
  assert.equal(DomainApp.of(defaults).base, '/domains/grit/issue', 'the route phase 3 serves');
  DomainApp.of(defaults).dispose();
});

scoped('DomainTable.app(): the docked view tells the app the route the shell mounted it at, ' +
  'and the address bar it was opened from', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const address = globalThis.location;
  globalThis.location = {pathname: '/apps/Grit', search: '?entity=i2'};
  let view;
  try {
    view = table.app();
    document.body.append(view.root);
    const a = DomainApp.of(view);
    assert.equal(a.base, '/domains/grit/issue', 'the default base until the view is docked');
    // the shell prepends its app route to every path the view reports (`View.path`)
    const own = view.path;
    Object.defineProperty(view, 'path', {get: () => `/apps/Grit${own}`, set: () => {}, configurable: true});
    grok.events.onViewAdded.fire(view);
    await flush();
    assert.equal(a.base, '/apps/Grit', 'the app rebased onto the shell\'s route');
    assert.equal(DomainApp.baseOf('grit.issue'), '/apps/Grit', 'find-or-activate looks it up there');
    assert.equal(view.acceptsPath('/apps/Grit'), true, 'and claims the route');
    assert.equal(a.page.value, 'entity', 'the cold deep link reached the app');
    assert.equal(a.entity.value, 'i2');
    a.dispose();
  } finally {
    globalThis.location = address;
  }
});

scoped('DomainTable.app(): the shell rewrites the view path only AFTER onViewAdded — ' +
  'the route is derived again, and the cold deep link still lands', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const address = globalThis.location;
  globalThis.location = {pathname: '/apps/Grit', search: '?entity=i2'};
  let view;
  try {
    view = table.app();
    document.body.append(view.root);
    const a = DomainApp.of(view);
    grok.events.onViewAdded.fire(view);
    await flush();
    assert.equal(a.base, '/domains/grit/issue', 'nothing to derive from yet');
    assert.equal(a.page.value, 'list', 'and no page restored');
    assert.equal(a.listSource.rows.items.value.length, 3, 'the list loaded all the same');
    const own = view.path;
    Object.defineProperty(view, 'path', {get: () => `/apps/Grit${own}`, set: () => {}, configurable: true});
    grok.events.onCurrentViewChanged.fire({args: {current: view}});
    await flush();
    assert.equal(a.base, '/apps/Grit', 'the rebase that lost the race lands here');
    assert.equal(DomainApp.baseOf('grit.issue'), '/apps/Grit');
    assert.equal(view.acceptsPath('/apps/Grit'), true);
    assert.equal(view.acceptsPath('/domains/grit/issue'), true, 'the fallback route is claimed too');
    assert.equal(view.acceptsPath('/apps/Other'), false);
    assert.equal(a.page.value, 'entity', 'and the cold deep link reached the app');
    assert.equal(a.entity.value, 'i2');
    grok.events.onCurrentViewChanged.fire({args: {current: view}});
    await flush();
    assert.equal(a.page.value, 'entity', 'replayed once');
    a.dispose();
  } finally {
    globalThis.location = address;
  }
});

scoped('the ambient session is the app\'s; New is hidden without insert', async () => {
  backends.domain = backend({access: {can: {view: true, insert: false, edit: true, delete: true, share: true},
    fields: {title: 'editable'}}});
  const table = await domains.table('grit.issue');
  const session = new SharedSession();
  const a = SharedSession.runWith(session, () => domains.app({table, base: BASE}));
  assert.equal(a.session, session);
  assert.equal(a.listSource.session, session);
  await flush();
  const ribbon = a.ribbon();
  assert.equal(ribbon[0][0].root.hidden, true);
  a.dispose();
});

scoped('app({app: Sub}): the subclass owns the ribbon, a preset writes the query with $me bound, a shortcut runs an action',
  async () => {
    backends.domain = backend();
    const table = await domains.table('grit.issue');
    const ran = [];
    table.actions.add({name: 'Escalate', run: (row) => ran.push(row.id)});
    class Sub extends DomainApp {
      shortcuts = {'Ctrl+Shift+E': 'Escalate'};
      ribbon() {
        const [main, tools] = super.ribbon();
        this.mine = this.presets(['Mine', 'reporter = $me'], ['Done', 'done = true']);
        return [main, [this.mine, ...tools]];
      }
    }
    const view = table.app({path: BASE, app: Sub});
    document.body.append(view.root);
    const a = DomainApp.of(view);
    assert.ok(a instanceof Sub);
    assert.equal(view.ribbonPanels[1].length, 3, 'the presets beside the search box and the filters');
    await flush();
    assert.equal(a.mine.root.dataset.u2, 'domain-presets');
    assert.deepEqual(a.mine.selected.value, [], 'no preset matches the empty query');
    const buttons = a.mine.root.querySelectorAll('.u2-btn-group-item');
    assert.deepEqual([...buttons].map((b) => b.textContent), ['Mine', 'Done']);
    fire(buttons[0], 'click');
    await flush();
    assert.equal(a.listSource.query.value, 'reporter = "user-1"', 'the stub shell\'s user, bound as $me');
    assert.deepEqual(a.mine.selected.value, ['0'], 'the preset in force is pressed');
    assert.equal(a.path.value, `${BASE}?q=${encodeURIComponent('reporter = "user-1"')}`);
    fire(buttons[0], 'click');
    await flush();
    assert.equal(a.listSource.query.value, '', 'a press on the preset in force clears the query');
    assert.deepEqual(a.mine.selected.value, [], 'and un-presses it');
    assert.equal(a.path.value, BASE);
    a.listSource.query.value = 'done = true';
    await flush();
    assert.deepEqual(a.mine.selected.value, ['1'], 'a query set elsewhere selects its preset');
    await a.goTo('entity', 'i2');
    await flush();
    assert.equal(a.mine.root.hidden, true, 'a list-page control');
    assert.equal(a.getWidgetStatus().shortcuts['Ctrl+Shift+E'], 'Escalate');
    fire(a.root, 'keydown', {key: 'e', ctrlKey: true, shiftKey: true});
    assert.deepEqual(ran, ['i2'], 'the entity page\'s row');
    fire(a.root, 'keydown', {key: 'e', ctrlKey: true});
    assert.deepEqual(ran, ['i2'], 'another chord is not the shortcut');
    a.dispose();
  });

scoped('a refused save says why in the status bar and balloons once, with no form on the page', async () => {
  const {table, app: a} = await app();
  const off = table.validators.add('priority',
    (value, row) => value === 'high' && !row.reporter ? 'Assign before escalating' : null);
  a.listSource.rows.byKey('i2').priority = 'high';
  await flush();
  assert.equal(a.summary.value, '1 unsaved change');
  assert.equal(await a.session.save(), false);
  await flush();
  assert.equal(a.summary.value, 'Cannot save: Ibuprofen: Priority: Assign before escalating',
    'the refusal outranks the change count, and names the row it is about');
  assert.equal(a.listSource.problemRow.value, 'i2', 'the offending row is named');
  assert.equal(a.list.list.root.querySelector('[data-u2-row="i2"]')
    .classList.contains('u2-domain-list-invalid'), true, 'and marked in the list');
  assert.equal(document.body.querySelectorAll('.u2-notify-error').length, 1, 'one balloon');
  a.session.discard();
  await flush();
  assert.equal(a.summary.value, '3 issues', 'a discard takes the refusal back');
  off();
  a.dispose();
});

scoped('a query the backend refuses is the list\'s error, in the status bar and over the rows', async () => {
  const {app: a} = await app();
  assert.equal(await a.open('?q=validator'), true);
  await flush();
  assert.equal(a.listSource.state.value, 'error');
  assert.equal(a.summary.value, 'Expected an operator', 'the backend\'s refusal, in the status bar');
  assert.match(a.root.querySelector('.u2-domain-list-error').textContent, /Expected an operator/);
  assert.equal(document.body.querySelector('.u2-notify-error'), null, 'a load failure is not a balloon');
  a.dispose();
});

scoped('a ?q= the backend refuses is not a dead end: the text stays in the box, New stands, Clear filter brings the rows back',
  async () => {
    const {app: a} = await app();
    const [[add], [, filters]] = a.ribbon();
    document.body.append(add.root, filters.root);
    await flush();
    assert.equal(await a.open('?q=validator'), true);
    await flush();
    const box = filters.input.value;
    assert.equal(a.listSource.state.value, 'error');
    assert.equal(box.text.value, 'validator', 'the offending text is in the box, to be fixed in place');
    assert.equal(box.validity.value !== null, true, 'red, with the reason on the editor');
    assert.equal(add.root.hidden, false, 'New follows the access, never a query the server refused');
    assert.match(a.summary.value, /Expected an operator/, 'the refusal of the query in force is the news');
    const clear = [...a.root.querySelectorAll('.u2-domain-list-error button')]
      .find((b) => b.textContent === 'Clear filter');
    assert.notEqual(clear, undefined, 'Retry alone would only run it again');
    fire(clear, 'click');
    await flush();
    assert.equal(a.listSource.query.value, '');
    assert.equal(box.text.value, '', 'the box is cleared with the query');
    assert.equal(a.listSource.rows.items.value.length, 3);
    a.dispose();
  });

scoped('an emptied box applies the empty query even when the refused one left no tree behind', async () => {
  const {app: a} = await app();
  const [, [, filters]] = a.ribbon();
  document.body.append(filters.root);
  await flush();
  assert.equal(await a.open('?q=validator'), true);
  await flush();
  const box = filters.input.value;
  box.text.value = '';
  assert.equal(box.commit(), true);
  await flush();
  assert.equal(a.listSource.query.value, '');
  assert.equal(a.listSource.rows.items.value.length, 3);
  a.dispose();
});

scoped('a refused filter leaves the rows as they were: shown stale, the count standing down', async () => {
  const {app: a} = await app();
  const [, [, filters]] = a.ribbon();
  document.body.append(filters.root);
  await flush();
  assert.equal(a.summary.value, '3 issues');
  const box = filters.input.value;
  box.text.value = 'priority = High';
  assert.equal(box.commit(), false, 'the schema refuses it');
  await flush();
  assert.equal(a.listSource.rows.items.value.length, 3, 'the rows still answer the previous filter');
  assert.equal(a.list.root.classList.contains('u2-domain-list-stale'), true);
  assert.equal(a.summary.value, '—', 'the count does not confirm rows the box no longer reads as');
  box.text.value = 'done = true';
  assert.equal(box.commit(), true);
  await flush();
  assert.equal(a.list.root.classList.contains('u2-domain-list-stale'), false);
  assert.equal(a.summary.value, '1 issue');
  a.dispose();
});

scoped('Ctrl+S saves from the ribbon — outside the app root, where Tab from the last field lands — and once from inside the form',
  async () => {
    backends.domain = backend();
    const table = await domains.table('grit.issue');
    const view = table.app({path: BASE});
    document.body.append(view.root);
    const a = DomainApp.of(view);
    grok.shell.v = view;
    try {
      let saves = 0;
      const save = a.session.save.bind(a.session);
      a.session.save = () => {
        saves++;
        return save();
      };
      await a.goTo('entity', 'i2');
      await flush();
      // the shell mounts the ribbon outside the app root, as appView hands it over
      const [, saveButton] = view.ribbonPanels[0];
      document.body.append(saveButton);
      a.form.value.input('title').value.value = 'From the ribbon';
      await flush();
      assert.equal(a.session.isDirty.value, true);

      fire(saveButton, 'keydown', {key: 's', ctrlKey: true});
      await flush();
      await flush();
      assert.equal(saves, 1, 'the ribbon is within the app\'s reach');
      assert.equal(a.session.isDirty.value, false);
      assert.equal(backends.domain.tableSync('grit.issue').rows.find((r) => r.id === 'i2').title,
        'From the ribbon');

      a.form.value.input('title').value.value = 'From the form';
      await flush();
      fire(a.form.value.root.querySelector('.u2-input-editor'), 'keydown', {key: 'Enter', ctrlKey: true});
      await flush();
      await flush();
      assert.equal(saves, 2, 'the form and the app never save twice for one keystroke');
      assert.equal(a.session.isDirty.value, false);
    } finally {
      grok.shell.v = undefined;
      DomainApp.of(view)?.dispose();
    }
  });

scoped('a preset shows in the filter box as it was written, not as the id it binds to', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  class Sub extends DomainApp {
    ribbon() {
      const [main, tools] = super.ribbon();
      this.mine = this.presets(['Mine', 'reporter = $me'], ['Done', 'done = true']);
      return [main, [this.mine, ...tools]];
    }
  }
  const view = table.app({path: BASE, app: Sub});
  document.body.append(view.root);
  const a = DomainApp.of(view);
  await flush();
  await flush();
  const filters = a.ribbon()[1].find((c) => c.root?.dataset?.u2 === 'domain-filters');
  const box = filters.input.value;
  fire(a.mine.root.querySelectorAll('.u2-btn-group-item')[0], 'click');
  await flush();
  assert.equal(a.listSource.query.value, 'reporter = "user-1"', 'the query runs bound');
  assert.equal(box.text.value, 'reporter = $me', 'the box says what the preset says');
  assert.equal(box.commit(), true, 'and Enter on it changes nothing');
  assert.equal(a.listSource.query.value, 'reporter = "user-1"');
  fire(a.mine.root.querySelectorAll('.u2-btn-group-item')[1], 'click');
  await flush();
  assert.equal(box.text.value, 'done = true', 'a preset without a binding reads as itself');
  a.dispose();
});

scoped('open(?q=) on the list page goes through the gate: cancel keeps the query and the changes', async () => {
  const {app: a} = await app();
  a.listSource.rows.byKey('i1').title = 'Edited';
  assert.equal(a.session.isDirty.value, true);
  const answer = a.open(`${BASE}?q=${encodeURIComponent('done = false')}`);
  await flush();
  assert.notEqual(document.body.querySelector('.u2-dialog'), null, 'a deep link is a query change like any other');
  assert.equal(a.listSource.query.value, '', 'nothing written yet');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await answer, false);
  assert.equal(a.listSource.query.value, '');
  assert.equal(a.session.isDirty.value, true);
  const again = a.open(`${BASE}?q=${encodeURIComponent('done = false')}`);
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  assert.equal(await again, true);
  assert.equal(a.listSource.query.value, 'done = false');
  await flush();
  assert.equal(await a.open(`${BASE}?q=${encodeURIComponent('done = false')}`), true,
    'the query it is already under asks nothing');
  a.dispose();
});

scoped('a path is the app\'s at a segment boundary only: a sibling table\'s route is not its own', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const view = table.app({path: '/domains/grit/issue'});
  await flush();
  assert.equal(view.acceptsPath('/domains/grit/issue'), true);
  assert.equal(view.acceptsPath('/domains/grit/issue?entity=i1'), true);
  assert.equal(view.acceptsPath('/domains/grit/issue/anything'), true);
  assert.equal(view.acceptsPath('/domains/grit/issue_label'), false, 'another table, not a deeper path');
  DomainApp.of(view).dispose();
});

scoped('closing through the write-back is blocked: the view is clean and still being written', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const view = table.app({path: '/domains/grit/issue'});
  await flush();
  const a = DomainApp.of(view);
  let release;
  const real = backends.domain.saveAll.bind(backends.domain);
  backends.domain.saveAll = async (edits) => {
    const landed = await real(edits);
    await new Promise((resolve) => release = resolve);
    return landed;
  };
  a.listSource.rows.byKey('i1').title = 'Edited';
  const saving = a.session.save();
  await flush();
  assert.equal(a.session.isDirty.value, false, 'the batch landed: the rows read clean already');
  const e = {args: {view: {dart: view.dart}}, prevented: 0, preventDefault() {
    this.prevented++;
  }};
  grok.events.onViewRemoving.fire(e);
  await flush();
  assert.equal(e.prevented, 1, 'the removal is cancelled through that window');
  assert.equal(view.dart.closed, 0);
  assert.match(document.body.querySelector('.u2-notify-warning')?.textContent ?? '', /Wait for the batch/);
  release();
  assert.equal(await saving, true);
  a.dispose();
});
