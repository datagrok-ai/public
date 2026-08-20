/* Headless smoke tests: construct → interact → dispose for every exported component, with the
   live-scope count asserted back to its baseline after each one (the leak detector as a test).
   Runs on the node DOM shim, so anything that needs real layout is only constructed, not measured. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {
  signal, computed, batch, untracked, Scope, Control, bindText, bindValue, AsyncSource,
  Overlay, OVERLAY_CLOSE_EVENT, Tooltip, VirtualList, VirtualTree, Combobox, TabStrip, Menu,
  Splitter, TextInput, TextArea, BoolInput, ChoiceInput, MultiChoiceInput, NumberInput, Dialog,
  Accordion, Breadcrumbs, Toolbar, Form, panel, div, span, h1, button, link,
} from '../src/index.js';
import {Perf} from '../src/core/perf.js';

/** Every test runs against a clean document and must leave the live-scope count where it was. */
function smoke(name, body) {
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

function mount(component) {
  document.body.append(component.root);
  return component;
}

function texts(root, selector) {
  return root.querySelectorAll(selector).map((el) => el.textContent);
}

smoke('signals: value, computed, batch and untracked', () => {
  const count = signal(1);
  const doubled = computed(() => count.value * 2);
  assert.equal(doubled.value, 2);

  const seen = [];
  const scope = new Scope();
  scope.effect(() => seen.push(doubled.value));
  count.value = 3;
  assert.deepEqual(seen, [2, 6]);

  batch(() => {
    count.value = 4;
    count.value = 5;
  });
  assert.deepEqual(seen, [2, 6, 10]);
  assert.equal(untracked(() => count.value), 5);

  scope.dispose();
  count.value = 9;
  assert.equal(seen.length, 3, 'effects die with their scope');
  assert.equal(scope.isDisposed, true);
  scope.dispose();
});

smoke('Scope: cleanups run once, and immediately when already disposed', () => {
  const scope = new Scope();
  let released = 0;
  scope.own(() => released++);
  scope.dispose();
  scope.dispose();
  assert.equal(released, 1);
  scope.own(() => released++);
  assert.equal(released, 2);
});

smoke('bindText/bindValue: live text and echo-suppressed two-way binding', () => {
  const scope = new Scope();
  const name = signal('world');
  const label = document.createElement('span');
  bindText(scope, label, computed(() => `Hello, ${name.value}!`));
  assert.equal(label.textContent, 'Hello, world!');
  name.value = 'u2';
  assert.equal(label.textContent, 'Hello, u2!');

  const input = document.createElement('input');
  let raw = '';
  let writes = 0;
  Object.defineProperty(input, 'value', {
    get: () => raw,
    set: (v) => {
      raw = v;
      writes++;
    },
  });
  bindValue(scope, input, name);
  assert.equal(input.value, 'u2');
  name.value = 'signal';
  assert.equal(input.value, 'signal');

  const before = writes;
  raw = 'typed';
  fire(input, 'input');
  assert.equal(name.value, 'typed');
  assert.equal(writes, before, 'the effect does not write the value back mid-keystroke');

  scope.dispose();
  name.value = 'dead';
  assert.equal(input.value, 'typed');
  assert.equal(label.textContent, 'Hello, typed!');
});

smoke('Control.build: everything built inside is disposed with the result', () => {
  let inner;
  let released = 0;
  const value = signal('x');
  const outer = Control.build(() => {
    inner = new TextInput({label: 'Name'});
    Scope.ambient.own(() => released++);
    return [inner, panel([span(value)])];
  });
  assert.equal(outer.root.contains(inner.root), true);
  assert.equal(inner.scope.isDisposed, false);
  assert.equal(outer.root.textContent.includes('x'), true);

  outer.dispose();
  assert.equal(inner.scope.isDisposed, true, 'ambient adoption cascades');
  assert.equal(released, 1);
});

smoke('elements: factories bind signals and release their listeners on dispose', () => {
  const scope = new Scope();
  const name = signal('world');
  let clicks = 0;
  const root = Scope.runWith(scope, () => panel([
    h1(computed(() => `Hello, ${name.value}!`)),
    div([span('static'), link('docs', 'https://datagrok.ai')]),
    button('Go', () => clicks++),
  ]));
  const heading = root.querySelector('h1');
  assert.equal(heading.textContent, 'Hello, world!');
  name.value = 'u2';
  assert.equal(heading.textContent, 'Hello, u2!');

  const go = root.querySelector('button');
  fire(go, 'click');
  assert.equal(clicks, 1);
  scope.dispose();
  fire(go, 'click');
  assert.equal(clicks, 1, 'listeners owned by the scope are gone');
  assert.throws(() => span(name), /needs an owner/);
});

smoke('VirtualList: virtualizes a million rows, selects by click and by key', () => {
  const list = mount(new VirtualList({
    render: (item, index, row) => {
      row.setAttribute('data-value', String(item));
      return span(`Row ${item}`);
    },
  }));
  list.root.clientHeight = 220;
  list.setItems(Array.from({length: 1000000}, (_, i) => i));
  assert.ok(list.renderedCount > 0 && list.renderedCount < 30, `rendered ${list.renderedCount} rows`);

  const row = list.root.querySelector('.u2-list-row[data-index="5"]');
  fire(row, 'click');
  assert.equal(list.selectedIndex.value, 5);
  assert.equal(row.getAttribute('aria-selected'), 'true');
  assert.equal(list.root.getAttribute('aria-activedescendant'), row.id);

  fire(list.root, 'keydown', {key: 'End'});
  assert.equal(list.selectedIndex.value, 999999);
  assert.ok(list.root.scrollTop > 0);
  assert.ok(list.renderedCount < 30, 'the tail is virtualized too');

  // the shim never clamps scrollTop when the content shrinks, as a real scroller would
  list.root.scrollTop = 0;
  list.setItems(['a', 'b']);
  assert.equal(list.selectedIndex.value, -1, 'a shorter list drops a stale selection');
  assert.equal(list.renderedCount, 2);
  list.dispose();
});

smoke('Combobox: opens on ArrowDown, filters, commits with Enter', async () => {
  const combobox = mount(new Combobox({items: ['alpha', 'beta', 'gamma'], placeholder: 'Pick'}));
  const input = combobox.root.querySelector('input');
  input.focus();
  assert.equal(combobox.machineState.value, 'focused');

  fire(input, 'keydown', {key: 'ArrowDown'});
  await flush();
  assert.equal(combobox.isOpen.value, true);
  assert.equal(combobox.machineState.value, 'open');
  assert.equal(input.getAttribute('aria-expanded'), 'true');
  const popup = document.body.querySelector('.u2-combobox-popup');
  assert.deepEqual(texts(popup, '.u2-combobox-option'), ['alpha', 'beta', 'gamma']);

  input.value = 'ta';
  fire(input, 'input');
  await flush();
  assert.deepEqual(texts(popup, '.u2-combobox-option'), ['beta']);

  fire(input, 'keydown', {key: 'ArrowDown'});
  fire(input, 'keydown', {key: 'Enter'});
  assert.equal(combobox.value.value, 'beta');
  assert.equal(combobox.isOpen.value, false);
  assert.equal(popup.isConnected, false);
  combobox.dispose();
});

smoke('TabStrip: activation, lazy content built once, close', () => {
  let builds = 0;
  const tabs = mount(new TabStrip({tabs: [{id: 'a', label: 'Alpha', content: span('A body')}]}));
  tabs.addTab({
    id: 'b',
    label: 'Beta',
    closable: true,
    content: () => {
      builds++;
      return span('B body');
    },
  });
  assert.equal(tabs.activeTab.value, 'a');
  assert.equal(builds, 0, 'inactive content stays unbuilt');

  tabs.activeTab.value = 'b';
  assert.equal(builds, 1);
  tabs.activeTab.value = 'a';
  tabs.activeTab.value = 'b';
  assert.equal(builds, 1, 'content is built once');

  const alpha = tabs.root.querySelector('.u2-tabs-tab[data-id="a"]');
  fire(alpha, 'pointerdown', {button: 0});
  assert.equal(tabs.activeTab.value, 'a');
  assert.equal(alpha.getAttribute('aria-selected'), 'true');

  tabs.closeTab('b');
  assert.equal(tabs.onTabClosed.value, 'b');
  assert.equal(tabs.root.querySelectorAll('.u2-tabs-tab').length, 1);
  tabs.dispose();
});

smoke('Menu: builds at coordinates, activates an item, closes', async () => {
  let clicked = '';
  const menu = new Menu()
    .item('Copy', () => clicked = 'copy')
    .separator()
    .group('More')
    .item('Deep', () => clicked = 'deep')
    .endGroup();

  menu.show({x: 40, y: 60});
  await flush();
  assert.equal(menu.isOpen.value, true);
  assert.equal(menu.machineState.value, 'open');
  const panelEl = document.body.querySelector('.u2-menu');
  assert.deepEqual(texts(panelEl, '.u2-menu-label'), ['Copy', 'More']);

  fire(panelEl, 'keydown', {key: 'ArrowDown'});
  assert.equal(panelEl.getAttribute('aria-activedescendant'), panelEl.querySelector('.u2-menu-item').id);

  fire(panelEl.querySelector('.u2-menu-label'), 'click');
  assert.equal(clicked, 'copy');
  assert.equal(menu.isOpen.value, false);
  assert.equal(menu.machineState.value, 'closed');
  assert.equal(panelEl.isConnected, false);

  menu.show({x: 10, y: 10});
  await flush();
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(menu.isOpen.value, false);
  menu.close();
});

smoke('Overlay: one layer host, idempotent close, close event', async () => {
  const scope = new Scope();
  const anchor = document.createElement('button');
  document.body.append(anchor);
  const content = document.createElement('div');
  let closed = 0;
  content.addEventListener(OVERLAY_CLOSE_EVENT, () => closed++);

  const close = Overlay.show(anchor, content, scope);
  await flush();
  assert.equal(content.isConnected, true);
  assert.equal(content.parentNode, Overlay.host);
  assert.ok(Number(content.style.zIndex) > 0);

  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(closed, 1);
  assert.equal(content.isConnected, false);
  close();
  assert.equal(closed, 1, 'close is idempotent');
  scope.dispose();
});

smoke('Splitter: keyboard resize clamps at minSize', () => {
  const splitter = mount(new Splitter([span('left'), span('right')],
    {direction: 'horizontal', minSize: 100}));
  splitter.root.clientWidth = 400;
  assert.deepEqual(splitter.sizes.value, [0.5, 0.5]);

  const sash = splitter.root.querySelector('.u2-sash');
  fire(sash, 'keydown', {key: 'ArrowRight'});
  assert.ok(Math.abs(splitter.sizes.value[0] - 0.52) < 1e-9);

  for (let i = 0; i < 40; i++)
    fire(sash, 'keydown', {key: 'ArrowRight', shiftKey: true});
  const [left, right] = splitter.sizes.value;
  assert.ok(left <= 0.75 + 1e-9, `left pane clamped at ${left}`);
  assert.ok(right >= 0.25 - 1e-9, 'the other pane never goes below minSize');
  assert.ok(Math.abs(left + right - 1) < 1e-9);

  fire(sash, 'dblclick');
  assert.deepEqual(splitter.sizes.value, [0.5, 0.5]);
  splitter.dispose();
});

smoke('VirtualTree: expands a branch and renames a node', () => {
  const renamed = [];
  const tree = mount(new VirtualTree({onRename: (node, label) => renamed.push([node.id, label])}));
  tree.root.querySelector('.u2-list').clientHeight = 220;
  tree.setRoots([
    {id: 'a', label: 'Alpha', children: [{id: 'a1', label: 'Alpha one'}]},
    {id: 'b', label: 'Beta'},
  ]);
  assert.deepEqual(texts(tree.root, '.u2-tree-label'), ['Alpha', 'Beta']);

  const twistie = tree.root.querySelector('.u2-list-row[data-index="0"] .u2-tree-twistie');
  fire(twistie, 'click');
  assert.equal(tree.expanded.value.has('a'), true);
  assert.deepEqual(texts(tree.root, '.u2-tree-label'), ['Alpha', 'Alpha one', 'Beta']);
  assert.equal(tree.root.querySelector('.u2-list-row[data-index="0"]').getAttribute('aria-expanded'), 'true');

  const label = tree.root.querySelector('.u2-list-row[data-index="0"] .u2-tree-label');
  fire(label, 'dblclick');
  const editor = tree.root.querySelector('.u2-tree-rename');
  assert.equal(editor.value, 'Alpha');
  editor.value = 'Renamed';
  fire(editor, 'keydown', {key: 'Enter'});
  assert.deepEqual(renamed, [['a', 'Renamed']]);
  assert.equal(tree.root.querySelector('.u2-list-row[data-index="0"] .u2-tree-label').textContent, 'Renamed');
  tree.dispose();
});

smoke('TextInput: validation flips validity, enabled toggles the editor', () => {
  const input = mount(new TextInput({label: 'Name'}));
  input.addValidator((v) => v.trim() ? null : 'Name is required');
  assert.equal(input.validity.value, 'Name is required');
  const editor = input.root.querySelector('input');
  assert.equal(editor.classList.contains('u2-invalid'), true);
  assert.equal(input.root.querySelector('.u2-input-error').textContent, 'Name is required');

  editor.value = 'Aspirin';
  fire(editor, 'input');
  assert.equal(input.value.value, 'Aspirin');
  assert.equal(input.validity.value, null);
  assert.equal(editor.classList.contains('u2-invalid'), false);

  input.enabled = false;
  assert.equal(editor.disabled, true);
  assert.equal(input.root.classList.contains('u2-input-disabled'), true);
  input.enabled = true;
  assert.equal(editor.disabled, false);
  input.dispose();
});

smoke('TextArea and BoolInput: bound editors and the switch variant', () => {
  const notes = mount(new TextArea({label: 'Notes', value: 'hello'}));
  const area = notes.root.querySelector('textarea');
  assert.equal(area.value, 'hello');
  area.value = 'edited';
  fire(area, 'input');
  assert.equal(notes.value.value, 'edited');

  const active = mount(new BoolInput({label: 'Active'}));
  const box = active.root.querySelector('input');
  box.checked = true;
  fire(box, 'change');
  assert.equal(active.value.value, true);
  active.value.value = false;
  assert.equal(box.checked, false);

  const toggle = mount(new BoolInput({label: 'Compact', switch: true}));
  const track = toggle.root.querySelector('.u2-input-switch');
  fire(track, 'keydown', {key: ' '});
  assert.equal(toggle.value.value, true);
  assert.equal(track.getAttribute('aria-checked'), 'true');

  notes.dispose();
  active.dispose();
  toggle.dispose();
});

smoke('NumberInput: parses, steps, and validates instead of clamping typed text', () => {
  const count = mount(new NumberInput({label: 'Count', mode: 'int', min: 0, max: 10, step: 2, spinner: true}));
  const editor = count.root.querySelector('input');

  editor.value = 'abc';
  fire(editor, 'input');
  assert.equal(count.validity.value, 'Not a number');
  assert.equal(count.value.value, null, 'the value keeps the last good number');

  editor.value = '4';
  fire(editor, 'input');
  assert.equal(count.value.value, 4);
  assert.equal(count.validity.value, null);

  fire(count.root.querySelector('.u2-number-spin'), 'click');
  assert.equal(count.value.value, 6);
  fire(editor, 'keydown', {key: 'ArrowUp', shiftKey: true});
  assert.equal(count.value.value, 10, 'stepping clamps at max');

  editor.value = '-5';
  fire(editor, 'input');
  fire(editor, 'blur');
  assert.equal(count.value.value, -5, 'out-of-range text is kept, as the platform input keeps it');
  assert.equal(editor.value, '-5');
  assert.equal(count.validity.value, 'Value must be at least 0');
  fire(count.root.querySelector('.u2-number-spin'), 'click');
  assert.equal(count.value.value, 0, 'stepping clamps back into range');
  assert.equal(count.validity.value, null);
  count.dispose();
});

smoke('ChoiceInput: setItems keeps the value only while it is still an item', () => {
  const series = mount(new ChoiceInput({label: 'Series', items: ['a', 'b', 'c'], value: 'b'}));
  const select = series.root.querySelector('select');
  assert.equal(select.value, 'b');
  assert.equal(select.children.length, 4, 'nullable adds the empty option');

  series.setItems(['b', 'd']);
  assert.equal(series.value.value, 'b');
  series.setItems(['x', 'y']);
  assert.equal(series.value.value, null);
  assert.equal(select.value, '');

  select.value = 'y';
  fire(select, 'change');
  assert.equal(series.value.value, 'y');
  series.dispose();
});

smoke('MultiChoiceInput: checked boxes rebuild the array, and follow it back', () => {
  const tags = mount(new MultiChoiceInput({label: 'Tags', items: ['a', 'b', 'c']}));
  const boxes = tags.root.querySelectorAll('.u2-input-checkbox');
  assert.deepEqual(tags.value.value, []);

  boxes[0].checked = true;
  boxes[2].checked = true;
  fire(boxes[0], 'change');
  assert.deepEqual(tags.value.value, ['a', 'c']);

  tags.value.value = ['b'];
  assert.deepEqual(boxes.map((b) => b.checked), [false, true, false]);
  tags.dispose();
});

smoke('Dialog: shows, runs OK and restores focus', () => {
  const trigger = document.createElement('button');
  document.body.append(trigger);
  trigger.focus();

  let confirmed = 0;
  const name = new TextInput({label: 'Name'});
  const dialog = Dialog.create('Rename').add(name).onOK(() => confirmed++);
  dialog.show({modal: true});
  assert.equal(dialog.isOpen.value, true);
  assert.equal(dialog.root.isConnected, true);
  assert.equal(dialog.root.getAttribute('aria-modal'), 'true');
  assert.equal(document.activeElement, name.root.querySelector('input'), 'focus moves into the dialog');

  const ok = dialog.root.querySelectorAll('button').find((b) => b.textContent === 'OK');
  fire(ok, 'click');
  assert.equal(confirmed, 1);
  assert.equal(dialog.isOpen.value, false);
  assert.equal(dialog.root.isConnected, false);
  assert.equal(document.activeElement, trigger, 'focus goes back to the opener');

  dialog.dispose();
  name.dispose();
});

smoke('Accordion: lazy pane content is built once and disposed with the pane', () => {
  let builds = 0;
  let released = 0;
  const accordion = mount(new Accordion());
  const pane = accordion.addPane('Details', () => {
    builds++;
    const content = new Control();
    content.own(() => released++);
    return content.root;
  });
  assert.equal(builds, 0);

  fire(accordion.root.querySelector('.u2-accordion-header'), 'click');
  assert.equal(pane.expanded.value, true);
  assert.equal(builds, 1);
  assert.equal(pane.root.querySelector('.u2-accordion-content').childNodes.length, 1);

  pane.expanded.value = false;
  pane.expanded.value = true;
  assert.equal(builds, 1, 'collapsed content is hidden, not rebuilt');

  accordion.removePane('Details');
  assert.equal(released, 1, 'removePane disposes what the builder made');
  assert.equal(accordion.getPane('Details'), undefined);
  accordion.dispose();
});

smoke('Breadcrumbs: renders a path, reports clicks, and setItems replaces it', () => {
  const clicks = [];
  const crumbs = mount(new Breadcrumbs({
    items: ['Home', 'Data', 'Table'],
    onClick: (index, item) => clicks.push([index, item]),
  }));
  assert.deepEqual(texts(crumbs.root, '.u2-breadcrumbs-item'), ['Home', '…', 'Data']);
  assert.equal(crumbs.root.querySelector('.u2-breadcrumbs-current').textContent, 'Table');

  fire(crumbs.root.querySelectorAll('.u2-breadcrumbs-item')[2], 'click');
  assert.deepEqual(clicks, [[1, 'Data']]);

  crumbs.setItems(['Projects', 'u2']);
  assert.deepEqual(crumbs.items.value, ['Projects', 'u2']);
  assert.deepEqual(texts(crumbs.root, '.u2-breadcrumbs-item'), ['Projects']);
  assert.equal(crumbs.root.querySelector('.u2-breadcrumbs-current').textContent, 'u2');
  crumbs.dispose();
});

smoke('Toolbar: a toggle mirrors its signal both ways', () => {
  let saves = 0;
  const bold = signal(false);
  const toolbar = mount(new Toolbar({ariaLabel: 'Format'}));
  toolbar.addButton('Save', () => saves++).addSeparator().addToggle('Bold', bold);

  const buttons = toolbar.root.querySelectorAll('.u2-toolbar-btn');
  const save = buttons.find((b) => b.textContent.includes('Save'));
  const toggle = buttons.find((b) => b.textContent.includes('Bold'));
  fire(save, 'click');
  assert.equal(saves, 1);

  fire(toggle, 'click');
  assert.equal(bold.value, true);
  assert.equal(toggle.getAttribute('aria-pressed'), 'true');
  assert.equal(toggle.classList.contains('u2-toolbar-btn-pressed'), true);

  bold.value = false;
  assert.equal(toggle.getAttribute('aria-pressed'), 'false');
  assert.equal(toggle.classList.contains('u2-toolbar-btn-pressed'), false);
  toolbar.dispose();
});

smoke('Form: aggregates validity, focuses the first invalid input, round-trips values', () => {
  const name = new TextInput({label: 'Name'});
  name.addValidator((v) => v.trim() ? null : 'Name is required');
  const age = new NumberInput({label: 'Age', mode: 'int'});
  const form = mount(new Form().addAll([name, age]));

  assert.equal(form.validity.value, 'Name is required');
  assert.equal(form.validate(), false);
  assert.equal(document.activeElement, name.root.querySelector('input'));

  name.value.value = 'Aspirin';
  assert.equal(form.validity.value, null);
  assert.equal(form.validate(), true);
  assert.deepEqual(form.getValues(), {Name: 'Aspirin', Age: null});

  form.setValues({Name: 'Ibuprofen', Age: 3});
  assert.equal(name.value.value, 'Ibuprofen');
  assert.equal(age.value.value, 3);
  assert.equal(form.inputs.length, 2);

  form.dispose();
  name.dispose();
  age.dispose();
});

smoke('Tooltip: bind shows on hover, hide clears it, dispose unbinds', async () => {
  const host = new Control();
  const target = document.createElement('button');
  document.body.append(target);
  Tooltip.bind(target, () => 'Compound name', host.scope);

  fire(target, 'pointerenter');
  assert.equal(Tooltip.isVisible, false, 'hover is delayed');
  Tooltip.hide();

  Tooltip.show('warm-up', target);
  await flush();
  assert.equal(Tooltip.isVisible, true);
  Tooltip.hide();
  assert.equal(Tooltip.isVisible, false);

  fire(target, 'pointerenter');
  assert.equal(Tooltip.isVisible, true, 'the warm window switches with no delay');
  assert.equal(Tooltip.root.textContent, 'Compound name');
  fire(target, 'pointerleave');
  assert.equal(Tooltip.isVisible, false);

  host.dispose();
  fire(target, 'pointerenter');
  assert.equal(Tooltip.isVisible, false, 'listeners died with the scope');
});

smoke('AsyncSource: debounces, aborts the superseded call, retries and reports empty', async () => {
  const calls = [];
  const source = new AsyncSource((query, abort) => {
    const call = {query, abort};
    const promise = new Promise((resolve, reject) => Object.assign(call, {resolve, reject}));
    calls.push(call);
    return promise;
  }, {debounceMs: 0});

  assert.equal(source.state.value.kind, 'idle');
  source.query('a');
  assert.equal(source.state.value.kind, 'loading');
  assert.equal(calls.length, 0, 'the fetch waits for the debounce');
  await flush();
  assert.equal(calls.length, 1);

  source.query('ab');
  await flush();
  assert.equal(calls.length, 2);
  assert.equal(calls[0].abort.aborted, true, 'the superseded call is aborted');
  calls[0].resolve(['stale']);
  await flush();
  assert.equal(source.state.value.kind, 'loading', 'an aborted result is ignored');

  calls[1].reject(new Error('boom'));
  await flush();
  assert.deepEqual(source.state.value, {kind: 'error', message: 'boom'});

  source.retry();
  await flush();
  assert.equal(calls.length, 3);
  assert.equal(calls[2].query, 'ab');
  calls[2].resolve(['beta']);
  await flush();
  assert.deepEqual(source.state.value, {kind: 'ready', items: ['beta']});

  source.query('zz');
  await flush();
  calls[3].resolve([]);
  await flush();
  assert.equal(source.state.value.kind, 'empty');
  source.dispose();
});

smoke('Perf: times construction, disposes the tree, accumulates the report', async () => {
  const before = Perf.report().length;
  let inner;
  const entry = await Perf.measure('two inputs', () => Control.build(() => {
    inner = new TextInput({label: 'Name'});
    return [inner, new BoolInput({label: 'Flag'})];
  }));
  assert.ok(entry.construct > 0);
  assert.equal(entry.paint, 0, 'paint needs a compositor');
  assert.equal(entry.total, entry.construct);
  assert.equal(inner.scope.isDisposed, true, 'measure disposes what it built');

  const report = Perf.report();
  assert.equal(report.length, before + 1);
  assert.equal(report[report.length - 1].label, 'two inputs');
});
