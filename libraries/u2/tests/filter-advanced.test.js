/* FilterBuilder, advanced mode: nested groups rendered recursively, `not`, the mode guard both
   ways, the drag gesture over the pure drop model, and what a template locks. Every builder is
   disposed in `finally` — an overlay or a drag left open would keep the shim's loop alive. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal, Filters, FilterBuilder} from '../src/index.js';
import {TYPE} from 'datagrok-api/u2core';

const SCHEMA = Filters.schema([
  {name: 'name', type: TYPE.STRING, friendlyName: 'Name'},
  {name: 'age', type: TYPE.INT, min: 0, max: 120},
  {name: 'mw', type: TYPE.FLOAT},
  {name: 'sex', type: TYPE.STRING, choices: ['F', 'M']},
]);

const mounted = [];

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      Filters.resetIds('');
      await body();
    } finally {
      for (const component of mounted.splice(0))
        component.dispose();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function mount(component) {
  document.body.append(component.root);
  mounted.push(component);
  return component;
}

/** root f5 (and) [ f1: age > 30, f4 (or) [ f2: sex = F, f3: name like "a" ] ] */
function nested() {
  return Filters.group('and', [
    Filters.cond('age', '>', 30),
    Filters.group('or', [Filters.cond('sex', '=', 'F'), Filters.cond('name', 'like', 'a')]),
  ]);
}

function nodes(fb) {
  return fb.root.querySelectorAll('[data-u2-node]');
}

function node(fb, id) {
  return fb.root.querySelector(`[data-u2-node="${id}"]`);
}

function part(el, name) {
  return el.querySelector(`[data-u2-part="${name}"]`);
}

/** The header part of the group `el` itself, not of a nested group. */
function own(el, name) {
  return el.children[0].querySelector(`[data-u2-part="${name}"]`);
}

function layout(fb, rects) {
  for (const [id, [x, y, w, h]] of Object.entries(rects))
    (id === 'rows' ? fb.root.querySelector('[data-u2-part="rows"]') : node(fb, id)).rect = new DOMRect(x, y, w, h);
}

function drag(handle, steps, commit = true) {
  fire(handle, 'pointerdown', {button: 0, clientX: 5, clientY: 5});
  for (const [x, y] of steps)
    fire(handle, 'pointermove', {clientX: x, clientY: y});
  if (commit)
    fire(handle, 'pointerup', {});
}

smoke('advanced mode renders the group hierarchy: data-u2-node on groups, a header per group, rows inside', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: nested()}));
  assert.deepEqual(nodes(fb).map((n) => n.dataset.u2Node), ['f1', 'f4', 'f2', 'f3'],
    'document order: the group before its children');
  const group = node(fb, 'f4');
  assert.equal(group.dataset.u2, 'filter-group');
  assert.equal(group.classList.contains('u2-fb-group'), true);
  for (const name of ['handle', 'connector', 'not', 'add', 'add-group', 'remove'])
    assert.notEqual(own(group, name), null, name);
  assert.deepEqual(part(group, 'rows').children.map((r) => r.dataset.u2Node), ['f2', 'f3']);
  assert.equal(node(fb, 'f2').parentElement, part(group, 'rows'), 'a child row lives in its group host');
  assert.equal(own(group, 'connector').querySelector('[data-id="or"]').classList.contains('u2-btn-group-selected'),
    true);
  const header = fb.root.querySelector('.u2-fb-header');
  assert.equal(part(header, 'not').hidden, false, 'the root has not in advanced mode');
  assert.equal(part(header, 'add-group').hidden, false);
  assert.equal(part(header, 'mode').textContent, 'Simple');
  const handle = part(node(fb, 'f1'), 'handle');
  assert.equal(handle.hidden, false, 'handles show in advanced mode');
  assert.deepEqual([handle.getAttribute('role'), handle.getAttribute('aria-label')], ['button', 'Drag to move']);
  assert.equal(fb.root.querySelector('[data-u2-part="hint"]').hidden, true);
  assert.equal(fb.root.querySelector('.u2-fb').classList.contains('u2-fb-advanced'), true);
  assert.equal(fb.getWidgetStatus().parts.mode, part(header, 'mode'));
});

smoke('not toggles group.not on a sub-group and on the root; the query shows it', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: nested()}));
  const not = own(node(fb, 'f4'), 'not');
  fire(not, 'click');
  assert.equal(Filters.find(fb.value.peek(), 'f4').not, true);
  assert.equal(not.classList.contains('u2-fb-not-on'), true);
  assert.equal(not.getAttribute('aria-pressed'), 'true');
  assert.equal(node(fb, 'f4').classList.contains('u2-fb-negated'), true);
  assert.equal(fb.query.peek(), 'age > 30 and not (sex = "F" or name like "a")');
  fire(not, 'click');
  assert.equal('not' in Filters.find(fb.value.peek(), 'f4'), false, 'off is the key gone, not false');
  assert.equal(fb.query.peek(), 'age > 30 and (sex = "F" or name like "a")');
  fire(part(fb.root.querySelector('.u2-fb-header'), 'not'), 'click');
  assert.equal(fb.value.peek().not, true);
  assert.equal(fb.query.peek(), Filters.format(fb.value.peek()));
  assert.equal(fb.query.peek().startsWith('not ('), true);
});

smoke('a group header edits its own group: and|or, + condition, + group, − removes the subtree', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: nested()}));
  const group = node(fb, 'f4');
  fire(own(group, 'connector').querySelector('[data-id="and"]'), 'click');
  assert.equal(Filters.find(fb.value.peek(), 'f4').op, 'and');
  assert.equal(fb.value.peek().op, 'and', 'the root op is untouched');
  fire(own(group, 'add'), 'click');
  let sub = Filters.find(fb.value.peek(), 'f4');
  assert.equal(sub.nodes.length, 3);
  assert.equal(sub.nodes[2].property, 'name', 'the first offered property');
  assert.equal(part(group, 'rows').children.length, 3, 'rendered inside the group');
  fire(own(group, 'add-group'), 'click');
  sub = Filters.find(fb.value.peek(), 'f4');
  assert.equal(Filters.isGroup(sub.nodes[3]), true);
  assert.equal(sub.nodes[3].op, 'or', 'the other connector than the parent');
  assert.equal(part(group, 'rows').children[3].dataset.u2, 'filter-group', 'nested two levels');
  assert.equal(part(group, 'rows').children[3].querySelector('.u2-fb-none'), null,
    'only the root shows the empty hint');
  fire(own(group, 'remove'), 'click');
  assert.deepEqual(nodes(fb).map((n) => n.dataset.u2Node), ['f1']);
  assert.equal(fb.root.querySelector('.u2-fb-hint').hidden, true);
});

smoke('rows keep their elements across edits inside a group, and a group keeps its element across child edits', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: nested()}));
  const group = node(fb, 'f4');
  const row = node(fb, 'f3');
  const editor = part(row, 'value');
  const input = editor.querySelector('input');
  input.value = 'ab';
  fire(input, 'input');
  assert.equal(Filters.find(fb.value.peek(), 'f3').value, 'ab');
  assert.equal(node(fb, 'f4'), group);
  assert.equal(node(fb, 'f3'), row);
  assert.equal(part(node(fb, 'f3'), 'value'), editor);
  fb.value.value = Filters.insert(fb.value.peek(), 'f4', 0, Filters.cond('mw', '<', 500));
  // f6 went to the builder's own default group; the inserted condition is f7
  assert.deepEqual(part(group, 'rows').children.map((r) => r.dataset.u2Node), ['f7', 'f2', 'f3']);
  assert.equal(node(fb, 'f3'), row, 'the reorder moved only the new element');
});

smoke('mode guard: the header link switches to advanced; back to simple is refused while nested, ' +
  'flatten from the hint', () => {
  const mode = signal('simple');
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode, value: nested()}));
  const hint = fb.root.querySelector('[data-u2-part="hint"]');
  const link = fb.root.querySelector('[data-u2-part="mode"]');
  assert.equal(hint.hidden, false, 'simple mode cannot show the tree');
  assert.equal(node(fb, 'f4').classList.contains('u2-fb-nested'), true, 'a summary row');
  assert.equal(link.textContent, 'Advanced');
  fire(link, 'click');
  assert.equal(mode.peek(), 'advanced');
  assert.equal(hint.hidden, true);
  assert.equal(node(fb, 'f4').classList.contains('u2-fb-group'), true, 'the summary became a real group');
  assert.equal(node(fb, 'f4').classList.contains('u2-fb-nested'), false);
  assert.equal(link.textContent, 'Simple');
  fire(link, 'click');
  assert.equal(mode.peek(), 'advanced', 'refused');
  assert.equal(hint.hidden, false);
  assert.equal(part(hint, 'flatten').hidden, true, 'or and and cannot flatten');
  fire(own(node(fb, 'f4'), 'connector').querySelector('[data-id="and"]'), 'click');
  fire(link, 'click');
  assert.equal(mode.peek(), 'advanced');
  assert.equal(part(hint, 'flatten').hidden, false, 'same connectors: flatten is offered');
  fire(part(hint, 'flatten'), 'click');
  assert.equal(mode.peek(), 'simple');
  assert.equal(hint.hidden, true);
  assert.deepEqual(nodes(fb).map((n) => n.dataset.u2Node), ['f1', 'f2', 'f3']);
  assert.equal(fb.root.querySelector('.u2-fb-header [data-u2-part="not"]').hidden, true);
  assert.equal(fb.root.querySelector('.u2-fb-header [data-u2-part="add-group"]').hidden, true);
  assert.equal(part(node(fb, 'f1'), 'handle').hidden, true, 'no handles in simple mode');
});

smoke('mode guard: a template with allowAdvanced false hides the switch and refuses advanced', () => {
  const template = {root: Filters.group('and', [Filters.cond('age', '>', 30)]), allowAdvanced: false};
  const fb = mount(new FilterBuilder({schema: SCHEMA, template, value: Filters.applyTemplate(template)}));
  assert.equal(fb.root.querySelector('[data-u2-part="mode"]').hidden, true);
  assert.equal(fb.setMode('advanced'), false);
  assert.equal(fb.mode.peek(), 'simple');
  const hint = fb.root.querySelector('[data-u2-part="hint"]');
  assert.equal(hint.hidden, false);
  assert.equal(part(hint, 'flatten').hidden, true);
  assert.equal(fb.setMode('simple'), true);
});

smoke('mode guard: allowAdvanced false renders simple even when the mode says advanced', () => {
  const template = {root: nested(), allowAdvanced: false};
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', template,
    value: Filters.applyTemplate(template)}));
  assert.equal(fb.root.querySelector('.u2-fb').classList.contains('u2-fb-advanced'), false);
  assert.equal(node(fb, 'f4').classList.contains('u2-fb-nested'), true, 'the sub-group is a summary row');
  assert.equal(part(node(fb, 'f1'), 'handle').hidden, true);
  assert.equal(fb.root.querySelector('.u2-fb-header [data-u2-part="not"]').hidden, true);
  assert.equal(fb.getWidgetStatus().mode, 'simple', 'the status says what is rendered');
  assert.equal(fb.root.querySelector('[data-u2-part="hint"]').hidden, false, 'a nested tree in simple mode');
});

/** The nested() tree laid out top to bottom, the group's rows inset. */
function layoutNested(fb) {
  layout(fb, {rows: [0, 0, 200, 140], f1: [0, 0, 200, 40], f4: [0, 40, 200, 100], f2: [10, 60, 190, 40],
    f3: [10, 100, 190, 40]});
}

smoke('drag: the handle moves a row into a group after its last child; the indicator follows the target', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: nested()}));
  layoutNested(fb);
  const indicator = fb.root.querySelector('.u2-fb-drop');
  assert.equal(indicator.hidden, true);
  const handle = part(node(fb, 'f1'), 'handle');
  fire(handle, 'pointerdown', {button: 0, clientX: 5, clientY: 20});
  fire(handle, 'pointermove', {clientX: 6, clientY: 22});
  assert.equal(node(fb, 'f1').classList.contains('u2-fb-dragging'), false, 'under the threshold nothing starts');
  fire(handle, 'pointermove', {clientX: 100, clientY: 45});
  assert.equal(node(fb, 'f1').classList.contains('u2-fb-dragging'), true);
  assert.equal(indicator.hidden, false);
  assert.equal(indicator.classList.contains('u2-fb-drop-into'), true, 'over the group header: into it');
  assert.equal(indicator.style.top, '40px');
  assert.equal(indicator.style.height, '100px');
  fire(handle, 'pointermove', {clientX: 100, clientY: 125});
  assert.equal(indicator.classList.contains('u2-fb-drop-line'), true, 'below f3\'s midpoint: a line after it');
  assert.equal(indicator.classList.contains('u2-fb-drop-into'), false);
  assert.equal(indicator.style.top, '140px');
  assert.equal(indicator.style.height, '2px');
  const before = fb.value.peek();
  fire(handle, 'pointerup', {});
  const after = fb.value.peek();
  assert.notEqual(after, before);
  assert.deepEqual(after.nodes.map((n) => n.id), ['f4']);
  assert.deepEqual(Filters.find(after, 'f4').nodes.map((n) => n.id), ['f2', 'f3', 'f1']);
  assert.equal(indicator.hidden, true);
  assert.equal(node(fb, 'f1').classList.contains('u2-fb-dragging'), false);
  assert.equal(node(fb, 'f1').parentElement, part(node(fb, 'f4'), 'rows'), 'the row element moved with the node');
  assert.equal(fb.query.peek(), Filters.format(after));
});

smoke('drag: Escape cancels, a refused target shows nothing, a locked group takes nothing, ' +
  'a press without motion is no drag', () => {
  const template = {root: nested()};
  Filters.find(template.root, 'f4').lock = 'value';
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', template,
    value: Filters.applyTemplate(template)}));
  layoutNested(fb);
  const indicator = fb.root.querySelector('.u2-fb-drop');
  assert.equal(part(node(fb, 'f2'), 'handle').hidden, true, 'a locked row has no handle');
  assert.equal(own(node(fb, 'f4'), 'handle').hidden, true, 'nor has the locked group');
  const handle = part(node(fb, 'f1'), 'handle');
  assert.equal(handle.hidden, false);
  const before = fb.value.peek();
  drag(handle, [[100, 45]], false);
  assert.equal(indicator.hidden, true, 'the locked group is refused as a target');
  fire(handle, 'pointermove', {clientX: 100, clientY: 10});
  assert.equal(indicator.hidden, true, 'before itself is nothing to do');
  fire(handle, 'pointermove', {clientX: 100, clientY: 125});
  assert.equal(indicator.hidden, true, 'into the locked group by a line is refused too');
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(node(fb, 'f1').classList.contains('u2-fb-dragging'), false);
  fire(handle, 'pointerup', {});
  assert.equal(fb.value.peek(), before, 'nothing moved');
  fire(handle, 'pointerdown', {button: 0, clientX: 5, clientY: 20});
  fire(handle, 'pointerup', {});
  assert.equal(fb.value.peek(), before);
  fire(handle, 'pointerdown', {button: 2, clientX: 5, clientY: 20});
  fire(handle, 'pointermove', {clientX: 100, clientY: 125});
  assert.equal(indicator.hidden, true, 'a secondary button never drags');
});

smoke('drag: a group moves with its subtree, and can leave a root it then empties into a sibling group', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: Filters.group('and', [
    Filters.group('or', [Filters.cond('age', '>', 30)]),
    Filters.group('or', [Filters.cond('sex', '=', 'F')]),
  ])}));
  // f1 in f2, f3 in f4, root f5
  layout(fb, {rows: [0, 0, 200, 160], f2: [0, 0, 200, 80], f1: [10, 20, 190, 40], f4: [0, 80, 200, 80],
    f3: [10, 100, 190, 40]});
  drag(own(node(fb, 'f2'), 'handle'), [[100, 30], [100, 85]]);
  const after = fb.value.peek();
  assert.deepEqual(after.nodes.map((n) => n.id), ['f4']);
  assert.deepEqual(Filters.find(after, 'f4').nodes.map((n) => n.id), ['f3', 'f2']);
  assert.deepEqual(Filters.find(after, 'f2').nodes.map((n) => n.id), ['f1'], 'the subtree came along');
  assert.equal(node(fb, 'f1').parentElement, part(node(fb, 'f2'), 'rows'));
});

smoke('templates: a group lock is inherited by its rows, allowAdd hides both adds everywhere, ' +
  'pickers are narrowed in nested rows', () => {
  const template = {
    root: Filters.group('and', [
      Filters.cond('age', '>', 30),
      Filters.group('or', [Filters.cond('sex', '=', 'F'), Filters.cond('name', 'like', 'a')], {lock: 'value'}),
      Filters.group('or', [Filters.cond('mw', '<', 500)]),
    ]),
    allowedProperties: ['age', 'sex', 'mw'],
    allowedOperators: {mw: ['<', '>']},
    allowAdd: false,
  };
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', template,
    value: Filters.applyTemplate(template)}));
  const header = fb.root.querySelector('.u2-fb-header');
  assert.equal(part(header, 'add').hidden, true);
  assert.equal(part(header, 'add-group').hidden, true);
  const locked = node(fb, 'f4');
  assert.equal(locked.classList.contains('u2-fb-locked'), true);
  assert.equal(own(locked, 'remove').hidden, true);
  assert.equal(own(locked, 'add').hidden, true);
  assert.equal(own(locked, 'add-group').hidden, true);
  assert.equal(own(locked, 'not').disabled, true);
  assert.equal(own(locked, 'connector').querySelector('[data-id="and"]').disabled, true);
  for (const id of ['f2', 'f3']) {
    const row = node(fb, id);
    assert.equal(row.classList.contains('u2-fb-locked'), true, `${id} inherits the group lock`);
    assert.equal(part(row, 'prop').hidden, true, 'a locked picker reads as text');
    assert.equal(part(row, 'op').hidden, true);
    assert.equal(part(row, 'prop-text').hidden, false);
    assert.equal(part(row, 'remove').hidden, true);
  }
  assert.equal(part(node(fb, 'f2'), 'value').querySelector('select').disabled, false,
    'value lock: the value stays editable');
  const free = node(fb, 'f6');
  assert.equal(free.classList.contains('u2-fb-locked'), false);
  assert.equal(own(free, 'add').hidden, true, 'allowAdd covers every group');
  assert.equal(own(free, 'remove').hidden, false);
  const mw = node(fb, 'f5');
  assert.deepEqual(part(mw, 'prop').querySelectorAll('option').map((o) => o.value), ['age', 'mw', 'sex'],
    'allowedProperties narrows the nested picker, schema order kept');
  assert.deepEqual(part(mw, 'op').querySelectorAll('option').map((o) => o.value), ['>', '<'], 'registry order');
  assert.deepEqual(Filters.validate(fb.value.peek(), SCHEMA, undefined, template), []);
  const diff = Filters.diff(template.root, fb.value.peek());
  assert.deepEqual(diff, {changed: [], added: [], removed: [], moved: []});
});

smoke('drag: allowAdd false keeps a row under its parent (a reorder stays), a locked row never leaves ' +
  'its group', () => {
  const template = {root: nested(), allowAdd: false};
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', template,
    value: Filters.applyTemplate(template)}));
  layoutNested(fb);
  const indicator = fb.root.querySelector('.u2-fb-drop');
  const before = fb.value.peek();
  drag(part(node(fb, 'f1'), 'handle'), [[100, 125]]);
  assert.equal(fb.value.peek(), before, 'into the sub-group: refused');
  assert.equal(indicator.hidden, true);
  drag(part(node(fb, 'f3'), 'handle'), [[100, 65]]);
  assert.deepEqual(Filters.find(fb.value.peek(), 'f4').nodes.map((n) => n.id), ['f3', 'f2'], 'a reorder is a drop');
  assert.deepEqual(Filters.validate(fb.value.peek(), SCHEMA, undefined, template), [], 'and the template agrees');

  Filters.resetIds('');
  const locked = {root: nested()};
  Filters.find(locked.root, 'f4').lock = 'value';
  const fb2 = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', template: locked,
    value: Filters.applyTemplate(locked)}));
  layoutNested(fb2);
  const frozen = fb2.value.peek();
  drag(part(node(fb2, 'f2'), 'handle'), [[100, 10]]);
  assert.equal(fb2.value.peek(), frozen, 'a row of the locked group stays, even pressed through its hidden handle');
});

function documentListeners(type) {
  return document._listeners.get(type)?.length ?? 0;
}

smoke('disposing mid-drag releases the gesture, every row scope, the indicator and the document listeners', () => {
  const base = Scope.liveCount;
  const listeners = ['pointermove', 'pointerup', 'pointercancel', 'keydown'].map(documentListeners);
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', value: nested()}));
  layoutNested(fb);
  const before = fb.value.peek();
  const handle = part(node(fb, 'f1'), 'handle');
  drag(handle, [[100, 125]], false);
  assert.equal(fb.root.querySelector('.u2-fb-drop').hidden, false);
  assert.equal(documentListeners('keydown'), listeners[3] + 1);
  fb.dispose();
  assert.equal(Scope.liveCount, base);
  fire(document, 'pointerup', {});
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(fb.value.peek(), before, 'nothing dropped');
  assert.equal(document.body.querySelector('.u2-fb-drop'), null, 'the indicator went with the builder');
  assert.deepEqual(['pointermove', 'pointerup', 'pointercancel', 'keydown'].map(documentListeners), listeners);
});

smoke('a row disposed mid-drag by the other builder on the signal ends the gesture on the document', () => {
  const tree = signal(nested());
  const a = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', bind: tree}));
  const b = mount(new FilterBuilder({schema: SCHEMA, mode: 'advanced', bind: tree}));
  layoutNested(a);
  const listeners = ['pointermove', 'pointerup', 'keydown'].map(documentListeners);
  drag(part(node(a, 'f1'), 'handle'), [[100, 125]], false);
  assert.equal(a.root.querySelector('.u2-fb-drop').hidden, false);
  fire(part(node(b, 'f1'), 'remove'), 'click');
  const removed = tree.peek();
  assert.equal(node(a, 'f1'), null, 'the dragged row is gone from both');
  // the browser lets the capture go with the element; what follows reaches the document only
  fire(document, 'pointermove', {clientX: 100, clientY: 45});
  fire(document, 'pointerup', {});
  assert.equal(tree.peek(), removed, 'no drop of a node that is gone');
  assert.equal(a.root.querySelector('.u2-fb-drop').hidden, true);
  assert.deepEqual(['pointermove', 'pointerup', 'keydown'].map(documentListeners), listeners);
});
