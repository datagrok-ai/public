/* Item actions: the hover block renders the icon-bearing subset, the context menu the full list,
   and VirtualList wires right-click to it. The timestamp factory rides along — same recipe. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {timestamp} from '../src/core/elements.js';
import {rowActions, actionsMenu} from '../src/components/actions/actions.js';
import {VirtualList} from '../src/components/collections/list.js';
import {Menu} from '../src/components/navigation/menu.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {dfBindings} from '../src/sources/df-bindings.js';
import {DataFrame, FilterGroup, Property, WidgetDescriptor, platform} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {ADD_FILTER, registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');
const {Ribbon} = await import('../src/dg/designer/ribbon.js');
const {SAMPLES} = await import('../src/dg/designer/samples.js');
const {shell} = await import('datagrok-api/grok');

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function sampleActions(log) {
  return [
    {name: 'Copy', icon: 'copy', run: () => log.push('copy')},
    {name: 'Open ticket', icon: 'external-link-alt', run: () => log.push('open')},
    {name: 'Copy id', run: () => log.push('id')},
  ];
}

scoped('rowActions: icon buttons for icon-bearing actions only, tooltip carries the name', () => {
  const log = [];
  const block = rowActions(sampleActions(log));
  assert.equal(block.dataset.u2, 'row-actions');
  const buttons = [...block.querySelectorAll('[data-u2="icon-button"]')];
  assert.equal(buttons.length, 2, 'the icon-less action appears only in menus');
  assert.deepEqual(buttons.map((b) => b.getAttribute('aria-label')), ['Copy', 'Open ticket']);
  fire(buttons[0], 'click');
  assert.deepEqual(log, ['copy']);
});

scoped('actionsMenu: every action becomes an item, picking one runs it and closes', () => {
  const log = [];
  const menu = actionsMenu(sampleActions(log));
  menu.show({x: 10, y: 10});
  const items = [...document.querySelectorAll('[role="menuitem"]')];
  assert.deepEqual(items.map((el) => el.querySelector('.u2-menu-label').textContent),
    ['Copy', 'Open ticket', 'Copy id']);
  fire(items[2], 'click');
  assert.deepEqual(log, ['id']);
  assert.equal(menu.isOpen.value, false);
});

scoped('VirtualList contextActions: right-click selects the row and opens the full menu', () => {
  const log = [];
  const list = new VirtualList({
    itemHeight: 20,
    render: (item) => {
      const el = document.createElement('span');
      el.textContent = item;
      return el;
    },
    contextActions: (item) => [{name: `Act on ${item}`, run: () => log.push(item)}],
  });
  document.body.append(list.root);
  list.setItems(['a', 'b', 'c']);

  const row = list.root.querySelector('[data-index="1"]');
  fire(row.firstChild, 'contextmenu');
  assert.equal(list.selectedIndex.value, 1, 'right-click selects');
  const item = document.querySelector('[role="menuitem"]');
  assert.equal(item.textContent.trim(), 'Act on b');
  fire(item, 'click');
  assert.deepEqual(log, ['b']);
  list.dispose();
});

scoped('timestamp: short visual, full title, year only when not current, empty when invalid', () => {
  const now = new Date();
  const thisYear = new Date(now.getFullYear(), 7, 11, 14, 30);
  const el = timestamp(thisYear);
  assert.equal(el.className, 'u2-timestamp');
  assert.equal(el.textContent,
    thisYear.toLocaleDateString(undefined, {month: 'short', day: 'numeric'}));
  assert.equal(el.title, thisYear.toLocaleString(undefined,
    {month: 'short', day: 'numeric', year: 'numeric', hour: '2-digit', minute: '2-digit'}));

  const old = new Date(2019, 2, 5);
  assert.equal(timestamp(old).textContent,
    old.toLocaleDateString(undefined, {month: 'short', day: 'numeric', year: 'numeric'}));

  assert.equal(timestamp({toDate: () => thisYear}).textContent, el.textContent,
    'anything with toDate() — dayjs included — is accepted');
  assert.equal(timestamp('not a date').textContent, '');
  assert.equal(timestamp('not a date').title, '');
});

/* WO-V5 — the designer verbs that live outside the view: `ADD_FILTER` (VP-11) over a rendered
   FilterGroup double — the dialog OK'd, cancelled, and refused for want of a table — and the Open
   menu's sample list (VP-15), which carries the viewer sample only on a registry with the tags. */

const ORDERS = () => new DataFrame([{name: 'city', type: 'string'}, {name: 'total', type: 'double'}],
  [{city: 'Kyiv', total: 1240}], 'orders');

const dialogButton = (text) => [...document.querySelectorAll('.u2-dialog .u2-dialog-buttons button')]
  .find((b) => b.textContent === text);

function designer(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warning = shell.warning;
    const warnings = [];
    shell.warning = (message) => warnings.push(message);
    WidgetDescriptor.registry = [new WidgetDescriptor('Grid', []),
      new WidgetDescriptor('Filters', [new Property('filters', 'list', {subType: 'map'})])];
    const reg = new Registry();
    registerPlatformComponents(reg);
    const scope = new Scope();
    const filters = (node, df) => renderSpec({$schema: 'dg-ui/1', root: node},
      new SpecContext({data: df ? {orders: dfBindings(signal(df), scope)} : {}}), reg);
    try {
      await body({reg, scope, warnings, filters});
    } finally {
      scope.dispose();
      shell.warning = warning;
      WidgetDescriptor.registry = [];
      platform.reset();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

designer('ADD_FILTER: the picked column joins the filters as a set-prop, a histogram for a numeric one',
  async ({filters, warnings}) => {
    const df = ORDERS();
    const instance = filters({tag: 'u2-viewer-filters', name: 'filters', bind: {table: '$.orders'},
      props: {filters: [{type: 'categorical', column: 'city'}]}}, df);
    const node = instance.spec.root;
    const fg = instance.node('filters');
    assert.ok(fg instanceof FilterGroup);
    const subscriptions = df.liveSubscriptions();

    const pending = ADD_FILTER.produce(node, fg);
    const dialog = document.querySelector('.u2-dialog');
    assert.ok(dialog, 'the action asks first');
    const select = dialog.querySelector('select');
    assert.deepEqual([...select.querySelectorAll('option')].map((o) => o.value).filter((v) => v !== ''),
      ['total'], 'the frame\'s columns, minus the ones already filtered');
    assert.equal(dialogButton('OK').disabled, true, 'nothing picked yet');
    fire(select, 'keydown', {key: 'Enter'});
    assert.notEqual(document.querySelector('.u2-dialog'), null, 'Enter with nothing picked finishes nothing');
    select.value = 'total';
    fire(select, 'change');
    assert.equal(dialogButton('OK').disabled, false);
    fire(dialogButton('OK'), 'click');

    assert.deepEqual(await pending, {op: 'set-prop', name: 'filters',
      value: [{type: 'categorical', column: 'city'}, {type: 'histogram', column: 'total'}]});
    assert.deepEqual(node.props.filters, [{type: 'categorical', column: 'city'}],
      'pure: the node is the view\'s to patch');
    assert.equal(document.querySelector('.u2-dialog'), null);
    assert.equal(df.liveSubscriptions(), subscriptions, 'the picker let go of the frame');
    assert.deepEqual(warnings, []);
    instance.dispose();
  });

designer('ADD_FILTER: seeds from the live group — the platform builds panes the document does not name — and from ' +
  'the document where the live read is empty', async ({filters}) => {
  const df = ORDERS();
  const instance = filters({tag: 'u2-viewer-filters', name: 'filters', bind: {table: '$.orders'}}, df);
  const node = instance.spec.root;
  const fg = instance.node('filters');
  fg.props.filters = [{type: 'categorical', column: 'city'}];
  const add = async (column) => {
    const pending = ADD_FILTER.produce(node, fg);
    const select = document.querySelector('.u2-dialog select');
    select.value = column;
    fire(select, 'change');
    fire(dialogButton('OK'), 'click');
    return (await pending).value;
  };
  assert.deepEqual(await add('total'),
    [{type: 'categorical', column: 'city'}, {type: 'histogram', column: 'total'}], 'appended after the live panes');
  assert.equal(node.props?.filters, undefined, 'pure: the document is untouched');

  fg.props.filters = [];
  node.props = {filters: [{type: 'histogram', column: 'total'}]};
  assert.deepEqual(await add('city'),
    [{type: 'histogram', column: 'total'}, {type: 'categorical', column: 'city'}], 'the document\'s, as a fallback');
  instance.dispose();
});

designer('ADD_FILTER: Cancel resolves null and leaves nothing behind', async ({filters, warnings}) => {
  const df = ORDERS();
  const instance = filters({tag: 'u2-viewer-filters', name: 'filters', bind: {table: '$.orders'}}, df);
  const pending = ADD_FILTER.produce(instance.spec.root, instance.node('filters'));
  fire(dialogButton('CANCEL'), 'click');
  assert.equal(await pending, null);
  assert.equal(document.querySelector('.u2-dialog'), null);
  assert.deepEqual(warnings, []);
  instance.dispose();
});

designer('ADD_FILTER: Enter with nothing picked is not an OK — Escape afterwards resolves null', async ({filters}) => {
  const df = ORDERS();
  const instance = filters({tag: 'u2-viewer-filters', name: 'filters', bind: {table: '$.orders'}}, df);
  const pending = ADD_FILTER.produce(instance.spec.root, instance.node('filters'));
  const select = document.querySelector('.u2-dialog select');
  fire(select, 'keydown', {key: 'Enter'});
  assert.notEqual(document.querySelector('.u2-dialog'), null, 'still asking');
  fire(select, 'keydown', {key: 'Escape'});
  assert.equal(await pending, null);
  assert.equal(document.querySelector('.u2-dialog'), null);
  instance.dispose();
});

designer('ADD_FILTER: with no table bound it warns once and resolves null, asking nothing',
  async ({filters, warnings}) => {
    const instance = filters({tag: 'u2-viewer-filters', name: 'filters'});
    const node = instance.spec.root;
    assert.equal(await ADD_FILTER.produce(node, instance.nodes().get(node)), null);
    assert.equal(document.querySelector('.u2-dialog'), null);
    assert.deepEqual(warnings,
      ['Add filter for column…: no table to pick a column from — bind `table` to a source first']);
    instance.dispose();
  });

designer('the Open menu lists the viewer sample only where the registry has the viewer tags',
  ({reg, scope}) => {
    const samples = (registry) => {
      const host = {mode: document.createElement('div'), outlines: signal(true), editor: signal(undefined),
        effect: (fn) => scope.effect(fn), instance: () => registry && {registry}, revision: () => 0,
        actions: () => [], run() {}, open() {}, dump: () => null, report() {}, refocus() {}};
      const menu = Menu.popup();
      new Ribbon(host)._openMenu(menu);
      return menu._nodes.find((n) => n.kind === 'group' && n.label === 'Samples').children.map((n) => n.label);
    };
    const core = new Registry();
    registerAll(core);
    assert.deepEqual(samples(core), SAMPLES.map((s) => s.name));
    assert.equal(samples(core).length, 5);
    assert.deepEqual(samples(reg), [...SAMPLES.map((s) => s.name), 'Master–detail (grid + filters)']);
    assert.equal(samples(reg).length, 6);
    assert.equal(samples(undefined).length, 5, 'no instance yet: the core samples alone');
  });
