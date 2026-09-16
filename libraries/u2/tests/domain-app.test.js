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
import {allowedActions} from '../src/components/actions/actions.js';
import {Control} from '../src/core/component.js';
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
/** A ribbon control by its `data-u2`, wherever in the group it sits — the search box and the
 * filter travel inside one wrapping item, and a subclass puts its own controls beside it. */
function ribbonControl(group, u2) {
  for (const control of group) {
    const el = control.root.dataset.u2 === u2 ? control.root : control.root.querySelector(`[data-u2="${u2}"]`);
    if (el)
      return Control.forElement(el);
  }
  return undefined;
}

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
  assert.deepEqual(ribbon[0].map((c) => c.root.dataset.u2),
    ['new-button', 'save-button', 'discard-button', 'actions-menu', 'refresh-button']);
  // only the search box rides the ribbon: the shell's is one fixed 32px line with `overflow:
  // hidden`, and the query box is a row of the LIST PAGE instead
  assert.deepEqual(ribbon[1].map((c) => c.root.dataset.u2), ['domain-search']);
  assert.equal(ribbon[1][0].root.hidden, false);
  const listPage = a.root.querySelector('[data-u2-part="list-page"]');
  assert.equal(listPage.children[0] === a.filters.root, true, 'the query box is the list page first row');
  fire(ribbon[0][0].root, 'click');
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.entity.value, DomainApp.NEW);
  assert.equal(a.path.value, `${BASE}?entity=new`);
  assert.equal(ribbon[1][0].root.hidden, true, 'no search on the entity page');
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
    assert.equal(view.ribbonPanels[1].length, 2, 'the presets beside the search box');
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
    const [[add]] = a.ribbon();
    const filters = a.filters;
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
  const filters = a.filters;
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
  const filters = a.filters;
  document.body.append(filters.root);
  await flush();
  assert.equal(a.summary.value, '3 issues');
  const box = filters.input.value;
  box.text.value = 'priority = High';
  assert.equal(box.commit(), false, 'the schema refuses it');
  await flush();
  assert.equal(a.listSource.rows.items.value.length, 3, 'the rows still answer the previous filter');
  assert.equal(a.list.root.classList.contains('u2-domain-list-stale'), true);
  assert.match(a.summary.value, /Filter not applied|Unknown column/,
    'the count does not confirm rows the box no longer reads as: it says why');
  // and what it says is the box's own refusal, not a placeholder: a tooltip on a red box is not read
  assert.equal(a.filterProblem.value, box.problems.value[0].message);
  assert.equal(a.summary.value, a.filterProblem.value);
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
  const filters = a.filters;
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

scoped('trash mode: the ⋯ menu toggles ?trash=1, the rows are read-only with Restore, and back', async () => {
  const {app: a} = await app();
  const ribbon = a.ribbon();
  await flush();
  const [add, save, discard, more] = ribbon[0];
  assert.equal(more.root.dataset.u2, 'actions-menu');
  assert.equal(more.root.hidden, false, 'the delete grant is what the menu holds');
  // the ribbon is the shell's, and an anchored menu needs its button in the document
  document.body.append(more.root);
  fire(more.root, 'click');
  await flush();
  const items = [...document.querySelectorAll('[role="menuitem"] .u2-menu-label')].map((el) => el.textContent);
  assert.deepEqual(items, ['Import…', 'Bulk edit…', 'Trash']);
  assert.equal(a.breadcrumbs.root.hidden, true, 'a plain list has no breadcrumb: the view title names the table');

  a.listSource.edit.value.markDeleted('i2');
  assert.equal(await a.session.save(), true);
  await flush();
  assert.equal(await a.setTrash(true), true);
  await flush();
  assert.equal(a.trash.value, true);
  assert.equal(a.listSource.deleted.value, 'only');
  assert.equal(a.path.value, `${BASE}?trash=1`);
  assert.deepEqual(a.listSource.rows.items.value.map((r) => r.title), ['Ibuprofen']);
  assert.equal(a.summary.value, '1 deleted issue', 'the table is named, as the live count names it');
  assert.deepEqual(crumbs(a), ['Issues', 'Trash']);
  assert.equal(a.breadcrumbs.root.hidden, false, 'the LIST page reads "Issues › Trash" too');
  assert.equal(save.root.hidden, true, 'nothing to save in the trash');
  assert.equal(discard.root.hidden, true);
  assert.equal(add.root.hidden, true, 'and nothing to insert: the access is narrowed');

  const row = a.listSource.rows.byKey('i2');
  assert.deepEqual(a.list.actionsFor(row).map((x) => x.name), ['Restore'], 'the row action Delete became Restore');
  assert.equal(a.list.root.querySelector('[data-u2-row="i2"]').classList.contains('u2-domain-list-deleted'), true);
  assert.equal(a.listSource.access.value.row(row).can('edit'), false, 'a deleted row is read-only');
  a.list.actionsFor(row)[0].run();
  await flush();
  assert.deepEqual(a.listSource.rows.items.value.map((r) => r.title), [], 'the restored row left the trash');

  assert.equal(await a.setTrash(false), true);
  await flush();
  assert.deepEqual(a.listSource.rows.items.value.map((r) => r.title), ['Aspirin', 'Ibuprofen', 'Naproxen']);
  assert.equal(a.path.value, BASE);
  assert.equal(save.root.hidden, false);
  a.dispose();
});

scoped('the trash reads newest-deleted first and gives the sort back on the way out', async () => {
  const {app: a} = await app();
  a.listSource.sort.value = 'title';
  await flush();
  assert.equal(await a.setTrash(true), true);
  await flush();
  assert.equal(a.listSource.sort.value, '!updated_on', 'the trash order while it is on');
  assert.equal(await a.setTrash(false), true);
  await flush();
  assert.equal(a.listSource.sort.value, 'title', 'the order the user chose is back');
  a.dispose();
});

scoped('?trash=1 round-trips through open(); without the delete grant the menu drops Trash', async () => {
  const {app: a} = await app({query: 'done = false'});
  await flush();
  assert.equal(await a.open('?trash=1&q=done%20%3D%20false'), true);
  await flush();
  assert.equal(a.trash.value, true);
  assert.equal(a.path.value, `${BASE}?q=done%20%3D%20false&trash=1`);
  assert.equal(await a.open(''), true);
  assert.equal(a.trash.value, false, 'a path is authoritative about the mode');
  a.dispose();

  backends.domain = backend({access: {can: {view: true, insert: true, edit: true, delete: false, share: false},
    fields: {title: 'editable'}}});
  const table = await domains.table('grit.issue');
  const b = domains.app({table, base: BASE, pageSize: 10});
  document.body.append(b.root);
  const ribbon = b.ribbon();
  await flush();
  assert.deepEqual(allowedActions(b.menuActions(), {access: b.listSource.access.value}).map((x) => x.name),
    ['Import…', 'Bulk edit…'], 'permission ⇒ hidden: Trash alone is gone');
  assert.equal(ribbon[0][3].root.hidden, false, 'the ⋯ button stands while anything in it is allowed');
  b.dispose();

  backends.domain = backend({access: {can: {view: true, insert: false, edit: false, delete: false, share: false},
    fields: {title: 'editable'}}});
  const c = domains.app({table: await domains.table('grit.issue'), base: BASE, pageSize: 10});
  document.body.append(c.root);
  const bare = c.ribbon();
  await flush();
  assert.equal(bare[0][3].root.hidden, true, 'and the button goes when nothing in it is allowed');
  c.dispose();
});

scoped('Refresh is in the ribbon only while the page is stale, and reloads through the gate', async () => {
  const memory = backend();
  let last = '2026-01-01T00:00:00Z';
  memory.tableSync('grit.issue').probe = () => Promise.resolve({count: 3, last});
  backends.domain = memory;
  const table = await domains.table('grit.issue');
  const a = domains.app({table, base: BASE, pageSize: 10});
  document.body.append(a.root);
  await flush();
  const reload = a.ribbon()[0][4];
  await flush();
  assert.equal(reload.root.dataset.u2, 'refresh-button');
  assert.equal(reload.root.hidden, true, 'nothing is behind the server yet');

  // the poll's own tick, taken by hand: the first one only records what the server looks like
  let tick;
  const saved = globalThis.setInterval;
  globalThis.setInterval = (fn) => (tick = fn, 0);
  a.listSource.live.value = true;
  await flush();
  globalThis.setInterval = saved;
  a.listSource.edit.value.setValue('i1', 'title', 'Edited');
  tick();
  await flush();
  last = '2026-02-02T00:00:00Z';
  tick();
  await flush();
  assert.equal(a.stale.value, true, 'the collection moved under pending changes');
  assert.equal(reload.root.hidden, false);

  fire(reload.root.querySelector('button') ?? reload.root, 'click');
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  await flush();
  assert.equal(a.session.isDirty.value, false);
  assert.equal(a.stale.value, false, 'the reload cleared it');
  assert.equal(reload.root.hidden, true);
  a.dispose();
});

scoped('entering and leaving the trash is ONE history entry each, and Back leaves it', async () => {
  const {app: a} = await app();
  const saved = globalThis.history;
  const pushes = [];
  globalThis.history = {pushState: (_state, _title, url) => pushes.push(url), replaceState: () => {}};
  try {
    assert.equal(await a.setTrash(true), true);
    await flush();
    assert.deepEqual(pushes, [`${BASE}?trash=1`], 'one entry for the move into the trash');
    assert.equal(await a.setTrash(false), true);
    await flush();
    assert.deepEqual(pushes, [`${BASE}?trash=1`, BASE], 'and one back out');
    // Back: the entry before the trash is the bare list, which `open` restores without pushing
    assert.equal(await a.setTrash(true), true);
    await flush();
    assert.equal(await a.open(BASE), true);
    await flush();
    assert.equal(a.trash.value, false, 'Back leaves the trash');
    assert.equal(pushes.length, 3, 'and the restore pushed nothing');
  } finally {
    globalThis.history = saved;
    a.dispose();
  }
});

scoped('entering the trash from the entity page is still one entry', async () => {
  const {app: a} = await app();
  assert.equal(await a.goTo('entity', 'i2'), true);
  await flush();
  const saved = globalThis.history;
  const pushes = [];
  globalThis.history = {pushState: (_state, _title, url) => pushes.push(url), replaceState: () => {}};
  try {
    assert.equal(await a.setTrash(true), true);
    await flush();
    assert.equal(a.page.value, 'list');
    assert.deepEqual(pushes, [`${BASE}?trash=1`], 'the page and the mode settle together');
  } finally {
    globalThis.history = saved;
    a.dispose();
  }
});

const UUID = '11111111-2222-3333-4444-555555555555';

scoped('under /domains a row is a path segment: the business key when unambiguous, the id when not', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.project');
  const a = domains.app({table, base: '/domains/grit/project', children: false});
  document.body.append(a.root);
  await flush();
  assert.equal(a.entityPath, true, 'the platform route addresses rows by a segment');
  assert.equal(await a.open('/domains/grit/project/GRIT'), true);
  await flush();
  assert.equal(a.page.value, 'entity');
  assert.equal(a.form.value.input('name').value.value, 'Grit', 'the business key resolved the row');
  assert.equal(a.path.value, '/domains/grit/project/GRIT');
  assert.equal(await a.goTo('list'), true);
  assert.equal(a.path.value, '/domains/grit/project', 'and back, with nothing left in the URL');
  assert.equal(await a.open('/domains/grit/project/nope'), true);
  await flush();
  assert.match(a.form.value.root.querySelector('[data-u2-part="empty"]').textContent, /not found/,
    'a key that matches nothing shows the not-found, as the Dart view does');
  assert.equal(a.keyOf({id: 'x', key: 'A-B'}), 'A-B', 'a single-column key keeps its dashes');
  assert.equal(a.keyOf({id: 'x', key: null}), 'x', 'a null component is not addressable');
  a.dispose();
});

scoped('a composite key is joined by "-" and read back by arity; a dash in a component falls back to the id',
  async () => {
    backends.domain = backend();
    const table = await domains.table('grit.issue');
    const a = domains.app({table, base: '/domains/grit/issue'});
    document.body.append(a.root);
    await flush();
    assert.equal(a.keyOf(a.listSource.rows.byKey('i1')), 'p1-1');
    assert.equal(a.keyOf({id: 'i9', project_id: 'p-1', number: 2}), 'i9', 'the split would be ambiguous');
    assert.equal(a.keyOf({id: 'i9', project_id: 'p1', number: null}), 'i9');
    assert.equal(await a.open('/domains/grit/issue/p1-1'), true);
    await flush();
    assert.equal(a.form.value.input('title').value.value, 'Aspirin', 'the int component is read as an int');
    assert.equal(a.path.value, '/domains/grit/issue/p1-1');
    assert.equal(await a.open('/domains/grit/issue/p1-1-2'), true, 'an arity that does not match is read as an id');
    await flush();
    assert.match(a.form.value.root.querySelector('[data-u2-part="empty"]').textContent, /not found/);
    a.dispose();
  });

scoped('a uuid segment is read as an id, and the path settles on the key the row carries', async () => {
  backends.domain = backend({rows: {project: [{id: UUID, key: 'GRIT', name: 'Grit'}], issue: []}});
  const table = await domains.table('grit.project');
  const a = domains.app({table, base: '/domains/grit/project', children: false});
  document.body.append(a.root);
  await flush();
  assert.equal(await a.open(`/domains/grit/project/${UUID}`), true);
  assert.equal(a.path.value, `/domains/grit/project/${UUID}`, 'the address as given, until the row is in');
  await flush();
  assert.equal(a.form.value.input('name').value.value, 'Grit');
  assert.equal(a.path.value, '/domains/grit/project/GRIT');
  a.dispose();
});

scoped('an app at /apps keeps emitting ?entity= (R2); a /domains link still reaches it after a rebase', async () => {
  const {app: a} = await app();
  assert.equal(a.entityPath, false);
  assert.equal(await a.goTo('entity', 'i2'), true);
  assert.equal(a.path.value, `${BASE}?entity=i2`);
  a.dispose();

  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const b = domains.app({table, base: '/domains/grit/issue'});
  document.body.append(b.root);
  await flush();
  b.rebase('/apps/Grit');
  assert.equal(b.entityPath, false, 'mounted at an app route, it is back to ?entity=');
  assert.equal(await b.open('/domains/grit/issue/p1-1'), true, 'the fallback route still reaches it');
  await flush();
  assert.equal(b.form.value.input('title').value.value, 'Aspirin');
  assert.equal(b.path.value, '/apps/Grit?entity=i1', 'and the URL it leaves behind is the id one');
  b.dispose();
});

scoped('?search= round-trips: read on open, written into the path beside the query (U8)', async () => {
  const {app: a} = await app();
  await flush();
  assert.equal(await a.open(`${BASE}?q=${encodeURIComponent('done = false')}&search=ibu`), true);
  await flush();
  assert.equal(a.listSource.search.value, 'ibu');
  assert.equal(a.path.value, `${BASE}?q=done%20%3D%20false&search=ibu`);
  assert.equal(a.listSource.rows.items.value.length, 1, 'the search narrowed the list');
  a.listSource.search.value = '';
  assert.equal(a.path.value, `${BASE}?q=done%20%3D%20false`);
  assert.equal(await a.open(`${BASE}?search=asp`), true);
  await flush();
  assert.equal(a.listSource.query.value, '', 'a path is authoritative about both');
  assert.equal(a.listSource.search.value, 'asp');
  assert.equal(a.path.value, `${BASE}?search=asp`);
  a.dispose();
});

scoped('history: one entry per user move, none while restoring, and the push lands before the shell mirrors',
  async () => {
    backends.domain = backend();
    const table = await domains.table('grit.issue');
    const saved = globalThis.history;
    const log = [];
    globalThis.history = {pushState: (_state, _title, url) => log.push(['push', url]),
      replaceState: (_state, _title, url) => log.push(['replace', url])};
    const pushes = () => log.filter(([kind]) => kind === 'push').map(([, url]) => url);
    let a;
    try {
      const view = table.app({path: BASE});
      document.body.append(view.root);
      a = DomainApp.of(view);
      Object.defineProperty(view, 'path', {configurable: true, get: () => view.dart.path,
        set: (x) => {
          log.push(['path', x]);
          view.dart.path = x;
        }});
      await flush();
      assert.equal(await a.goTo('entity', 'i2'), true);
      assert.deepEqual(log, [['push', `${BASE}?entity=i2`], ['path', `${BASE}?entity=i2`]],
        'the entry is pushed before the shell replaces it with the same URL');
      assert.equal(await a.goTo('list'), true);
      a.listSource.query.value = 'done = true';
      await flush();
      a.listSource.search.value = 'ibu';
      await flush();
      // a page move is an entry each; the query OPENS one for the text and the search, still in
      // the same burst, rewrites it — typing is not navigating (five keystrokes were five Backs)
      assert.deepEqual(pushes(), [`${BASE}?entity=i2`, BASE, `${BASE}?q=done%20%3D%20true`],
        'two page moves and one entry for the text');
      assert.deepEqual(log[log.length - 1], ['path', `${BASE}?q=done%20%3D%20true&search=ibu`]);
      assert.equal(log.some(([kind, url]) => kind === 'replace' &&
        url === `${BASE}?q=done%20%3D%20true&search=ibu`), true, 'the search rewrote that entry');
      assert.equal(await a.open(`${BASE}?entity=i1`), true);
      await flush();
      assert.equal(await a.open(BASE), true);
      await flush();
      assert.equal(pushes().length, 3, 'a restore from the address bar pushes nothing');
      const ribbon = a.ribbon();
      fire(ribbon[0][0].root, 'click');
      await flush();
      assert.equal(pushes().length, 4, 'New is a move');
      a.form.value.input('project_id').value.value = 'p1';
      a.form.value.input('title').value.value = 'Fresh';
      await flush();
      assert.equal(await a.session.save(), true);
      await flush();
      assert.equal(a.entity.value, backends.domain.tableSync('grit.issue').rows.find((r) => r.title === 'Fresh').id);
      assert.equal(pushes().length, 4, 'the draft that became a row is the same page');
    } finally {
      globalThis.history = saved;
      a?.dispose();
    }
  });

scoped('a cold deep link replays every parameter the app owns: ?trash=1 and ?search= too', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const address = globalThis.location;
  globalThis.location = {pathname: '/domains/grit/issue', search: '?trash=1&search=asp'};
  let a;
  try {
    const view = table.app({path: '/domains/grit/issue'});
    document.body.append(view.root);
    a = DomainApp.of(view);
    grok.events.onViewAdded.fire(view);
    await flush();
    assert.equal(a.trash.value, true, 'the mode the URL carried');
    assert.equal(a.listSource.search.value, 'asp');
    assert.equal(a.path.value, '/domains/grit/issue?trash=1&search=asp');
  } finally {
    globalThis.location = address;
    a?.dispose();
  }
});

scoped('Back restores the app from the address bar, and pushes nothing while doing it', async () => {
  backends.domain = backend();
  const table = await domains.table('grit.issue');
  const savedHistory = globalThis.history;
  const address = globalThis.location;
  const pushes = [];
  globalThis.history = {pushState: (_state, _title, url) => pushes.push(url), replaceState: () => {}};
  let a;
  try {
    const view = table.app({path: BASE});
    document.body.append(view.root);
    a = DomainApp.of(view);
    grok.shell.v = view;
    await flush();
    assert.equal(await a.goTo('entity', 'i2'), true);
    await flush();
    assert.equal(pushes.length, 1);
    globalThis.location = {pathname: BASE, search: ''};
    window.dispatchEvent(new Event('popstate'));
    await flush();
    assert.equal(a.page.value, 'list', 'the entry the URL carries is restored');
    assert.equal(pushes.length, 1, 'a restore pushes nothing');
    globalThis.location = {pathname: '/apps/Somebody/Else', search: '?entity=i1'};
    window.dispatchEvent(new Event('popstate'));
    await flush();
    assert.equal(a.page.value, 'list', 'another view\'s address is not the app\'s to read');
    grok.shell.v = undefined;
    globalThis.location = {pathname: BASE, search: '?entity=i1'};
    window.dispatchEvent(new Event('popstate'));
    await flush();
    assert.equal(a.page.value, 'list', 'and neither is any address while the app is not shown');
  } finally {
    globalThis.history = savedHistory;
    globalThis.location = address;
    grok.shell.v = undefined;
    a?.dispose();
  }
});

scoped('a selection is not a pending change: the unload gate stays down over a clean list', async () => {
  const {app: a} = await app();
  a.guardUnload();
  a.list.list.root.clientHeight = 400;
  fire(a.list.list.root, 'scroll');
  await flush();
  // the guard prevents the event; a prevented dispatch answers false
  const armed = () => fire(globalThis, 'beforeunload', {cancelable: true}) === false;
  assert.equal(armed(), false, 'nothing pending, nothing to hold the browser for');

  fire(a.list.root.querySelector('[data-u2-row="i1"]'), 'click');
  fire(a.list.root.querySelector('[data-u2-row="i3"]'), 'click', {ctrlKey: true});
  await flush();
  assert.deepEqual(a.listSource.selection.value.map((r) => r.id), ['i1', 'i3']);
  assert.equal(a.session.isDirty.value, false, 'picking rows writes nothing');
  assert.equal(armed(), false, 'so leaving a clean list asks nothing');

  a.listSource.rows.byKey('i1').title = 'Edited';
  await flush();
  assert.equal(armed(), true, 'an actual change does hold it');
  a.dispose();
});
