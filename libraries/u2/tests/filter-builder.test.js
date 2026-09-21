/* FilterBuilder, simple mode: rows keyed by node id, the immutable-edit contract, the and/or
   toggle, horizontal joins, locks, the per-kind value editors and the registration-free surface
   (query, problems, status). Popup-opening tests dispose the builder in `finally` — an overlay
   left open loops the animation-frame updater and OOMs node --test. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal, Filters, FilterBuilder, TextInput, Control} from '../src/index.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {TYPE} from 'datagrok-api/u2core';

const SCHEMA = Filters.schema([
  {name: 'name', type: TYPE.STRING, friendlyName: 'Name'},
  {name: 'age', type: TYPE.INT, min: 0, max: 120},
  {name: 'mw', type: TYPE.FLOAT},
  {name: 'big', type: TYPE.BIG_INT},
  {name: 'sex', type: TYPE.STRING, choices: ['F', 'M']},
  {name: 'active', type: TYPE.BOOL},
  {name: 'created', type: TYPE.DATE_TIME},
]);
const WITH_VALUES = Filters.schema([{name: 'status', type: TYPE.STRING}], {status: ['Open', 'Blocked', 'Closed']});

const mounted = [];

/** Every mounted control is disposed in `finally` — that closes any popup it opened. */
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

function rows(fb) {
  return fb.root.querySelectorAll('[data-u2-node]');
}

function part(row, name) {
  return row.querySelector(`[data-u2-part="${name}"]`);
}

function select(row, name) {
  return part(row, name).querySelector('select');
}

function pick(selectEl, value) {
  selectEl.value = value;
  fire(selectEl, 'change');
}

function type(input, text) {
  input.value = text;
  fire(input, 'input');
}

function tree(...conds) {
  return Filters.group('and', conds);
}

smoke('renders one row per condition, addressed by node id and part', () => {
  const fb = mount(new FilterBuilder({label: 'Criteria', schema: SCHEMA, name: 'fb',
    value: tree(Filters.cond('age', '>', 30), Filters.cond('name', 'like', 'an'))}));
  assert.equal(fb.root.dataset.u2, 'filter-builder');
  assert.equal(fb.root.dataset.u2Name, 'fb');
  const [a, b] = rows(fb);
  assert.deepEqual([a.dataset.u2Node, b.dataset.u2Node], ['f1', 'f2']);
  assert.equal(select(a, 'prop').value, 'age');
  assert.equal(select(a, 'op').value, '>');
  assert.equal(part(a, 'value').dataset.u2, 'number-input');
  assert.equal(part(a, 'value').querySelector('input').value, '30');
  assert.equal(part(a, 'remove').tagName, 'BUTTON');
  assert.equal(part(b, 'value').dataset.u2, 'text-input');
  assert.deepEqual(select(a, 'prop').querySelectorAll('option').map((o) => o.textContent),
    ['Name', 'age', 'mw', 'big', 'sex', 'active', 'created'], 'friendly names, no empty option');
  assert.deepEqual(select(a, 'op').querySelectorAll('option').map((o) => o.value),
    Filters.operators.for(SCHEMA.properties[1]).map((o) => o.id), 'operators of the property kind');
  assert.equal(fb.root.querySelector('.u2-fb-none'), null);
});

smoke('empty tree shows the empty hint (a link that adds); + adds a row on the first property, − removes it', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA}));
  assert.equal(rows(fb).length, 0);
  const none = fb.root.querySelector('.u2-fb-none');
  assert.equal(none.tagName, 'A');
  assert.equal(none.dataset.u2Part, 'none');
  fire(none, 'click');
  assert.equal(rows(fb).length, 1, 'the empty hint adds a condition');
  fb.removeNode(fb.value.peek().nodes[0].id);
  assert.equal(rows(fb).length, 0);
  const add = fb.root.querySelector('[data-u2-part="add"]');
  fire(add, 'click');
  assert.equal(rows(fb).length, 1);
  assert.equal(fb.root.querySelector('.u2-fb-none'), null);
  const cond = fb.value.peek().nodes[0];
  assert.equal(cond.property, 'name');
  assert.equal(cond.operator, Filters.operators.for(SCHEMA.properties[0])[0].id);
  assert.equal(cond.value, undefined);
  assert.deepEqual(fb.problems.peek().map((p) => p.code), ['missing-value']);
  assert.equal(fb.validity.peek(), fb.problems.peek()[0].message, 'validity is the first problem');
  assert.equal(rows(fb)[0].classList.contains('u2-fb-invalid'), true);
  assert.equal(rows(fb)[0].title, fb.problems.peek()[0].message);
  const problem = part(rows(fb)[0], 'problem');
  assert.equal(problem.hidden, false, 'the row tells its problem on a line of its own');
  assert.equal(problem.textContent, fb.problems.peek()[0].message);
  type(part(rows(fb)[0], 'value').querySelector('input'), 'x');
  assert.equal(problem.hidden, true, 'gone with the value');
  assert.equal(fb.validity.peek(), null);
  const id = fb.addCondition();
  assert.equal(rows(fb).length, 2);
  assert.equal(rows(fb)[1].dataset.u2Node, id);
  fire(part(rows(fb)[0], 'remove'), 'click');
  assert.deepEqual(rows(fb).map((r) => r.dataset.u2Node), [id]);
  fb.removeNode(id);
  assert.equal(rows(fb).length, 0);
  assert.equal(fb.validity.peek(), null);
});

smoke('allowAdd: false keeps the empty hint plain text', () => {
  const template = {root: Filters.group('and', []), allowAdd: false};
  const fb = mount(new FilterBuilder({schema: SCHEMA, template, value: Filters.applyTemplate(template)}));
  const none = fb.root.querySelector('.u2-fb-none');
  assert.equal(none.tagName, 'SPAN');
  fire(none, 'click');
  assert.equal(rows(fb).length, 0);
});

smoke('property pick resets the operator to the first applicable and clears the value; keeps it when ' +
  'applicable', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: tree(Filters.cond('age', '>', 30))}));
  const row = rows(fb)[0];
  pick(select(row, 'prop'), 'name');
  let cond = fb.value.peek().nodes[0];
  assert.equal(cond.property, 'name');
  assert.equal(cond.operator, '=', '> is not applicable to a string');
  assert.equal(cond.value, undefined);
  assert.equal(part(row, 'value').dataset.u2, 'text-input', 'the editor followed the kind');
  assert.equal(rows(fb)[0] === row, true, 'same row element');

  pick(select(row, 'prop'), 'mw');
  cond = fb.value.peek().nodes[0];
  assert.equal(cond.operator, '=', 'still applicable, kept');
  assert.equal(part(row, 'value').dataset.u2, 'number-input');
});

smoke('operator pick keeps the value across a same-shape swap and clears it across a shape change', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: tree(Filters.cond('age', '>', 30))}));
  const row = rows(fb)[0];
  const editor = part(row, 'value');
  pick(select(row, 'op'), '<=');
  assert.deepEqual(fb.value.peek().nodes[0], {id: 'f1', property: 'age', operator: '<=', value: 30});
  assert.equal(part(row, 'value') === editor, true, 'the number editor survived');
  pick(select(row, 'op'), 'between');
  assert.equal(fb.value.peek().nodes[0].value, undefined);
  assert.equal(part(row, 'value').dataset.u2, 'number-input');
  assert.equal(part(row, 'value2').dataset.u2, 'number-input');
  pick(select(row, 'op'), 'is null');
  assert.equal(part(row, 'value'), null, 'arity 0 has no editor');
  assert.equal(fb.query.peek(), 'age = null');
});

smoke('a value edit writes a new root, keeps the untouched sibling node and the row element', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA,
    value: tree(Filters.cond('age', '>', 30), Filters.cond('name', 'like', 'an'))}));
  const before = fb.value.peek();
  const [row] = rows(fb);
  const editor = part(row, 'value');
  type(editor.querySelector('input'), '40');
  const after = fb.value.peek();
  assert.notEqual(after, before, 'a new root');
  assert.equal(after.nodes[1], before.nodes[1], 'the sibling keeps its identity');
  assert.equal(after.nodes[0].value, 40);
  assert.equal(rows(fb)[0] === row, true, 'the row element is the same');
  assert.equal(part(rows(fb)[0], 'value') === editor, true, 'and so is the editor being typed into');
  assert.equal(fb.query.peek(), 'age > 40 and name like "an"');
  assert.equal(fb.query.peek(), Filters.format(after));
});

smoke('external tree writes reconcile the rows: same ids keep their elements, new ones appear, gone ones leave', () => {
  const bound = signal(tree(Filters.cond('age', '>', 30), Filters.cond('name', 'like', 'an')));
  const fb = mount(new FilterBuilder({schema: SCHEMA, bind: bound}));
  const [a, b] = rows(fb);
  const added = Filters.cond('sex', '=', 'F');
  bound.value = Filters.insert(Filters.remove(bound.peek(), 'f1'), bound.peek().id, 0, added);
  const now = rows(fb);
  assert.deepEqual(now.map((r) => r.dataset.u2Node), [added.id, 'f2']);
  assert.equal(now[1] === b, true, 'the surviving row kept its element');
  assert.equal(now[0] === a, false);
  assert.equal(select(now[0], 'value').value, 'F', 'choices render as a select');
  bound.value = Filters.update(bound.peek(), 'f2', {value: 'ol'});
  assert.equal(part(rows(fb)[1], 'value').querySelector('input').value, 'ol', 'the editor follows the model');
});

smoke('the and|or toggle writes the root op, and the root op drives the toggle', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: tree(Filters.cond('age', '>', 30))}));
  const connector = fb.root.querySelector('[data-u2-part="connector"]');
  assert.equal(connector.querySelector('[data-id="and"]').classList.contains('u2-btn-group-selected'), true);
  fire(connector.querySelector('[data-id="or"]'), 'click');
  assert.equal(fb.value.peek().op, 'or');
  fb.value.value = Filters.update(fb.value.peek(), fb.value.peek().id, {op: 'and'});
  assert.equal(connector.querySelector('[data-id="and"]').classList.contains('u2-btn-group-selected'), true);
});

smoke('horizontal: rows inline with a join chip between them that flips the op; + at the end', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, orientation: 'horizontal',
    value: tree(Filters.cond('age', '>', 30), Filters.cond('name', 'like', 'an'), Filters.cond('sex', '=', 'F'))}));
  const editor = fb.root.querySelector('.u2-fb');
  assert.equal(editor.classList.contains('u2-fb-horizontal'), true);
  assert.equal(fb.root.querySelector('.u2-fb-header').hidden, true);
  const host = fb.root.querySelector('[data-u2-part="rows"]');
  const kinds = host.children.map((el) => el.classList.contains('u2-fb-join') ? 'join' :
    el.dataset.u2Part === 'add' ? 'add' : 'row');
  assert.deepEqual(kinds, ['row', 'join', 'row', 'join', 'row', 'add']);
  const joins = host.querySelectorAll('.u2-fb-join');
  assert.deepEqual(joins.map((j) => j.textContent), ['and', 'and']);
  assert.equal(joins[0].tagName, 'BUTTON', 'a real button: Enter and Space click it natively');
  assert.equal(joins[0].getAttribute('aria-label'), 'Switch and/or');
  assert.equal(joins[0].title, '', 'the tooltip service, not title');
  fire(joins[0], 'click');
  assert.equal(fb.value.peek().op, 'or');
  assert.deepEqual(host.querySelectorAll('.u2-fb-join').map((j) => j.textContent), ['or', 'or']);
  fb.removeNode('f2');
  assert.deepEqual(host.children.map((el) => el.dataset.u2Node ?? el.className),
    ['f1', 'u2-btn u2-fb-join', 'f3', 'u2-btn u2-icon-btn u2-fb-add']);
});

smoke('locks: value keeps the editor and hides remove, all disables everything, allowAdd hides +', () => {
  const template = {
    root: Filters.group('and', [
      Filters.cond('sex', '=', 'F', {lock: 'value'}),
      Filters.cond('created', '>', {span: '-30d'}, {lock: 'all'}),
      Filters.cond('age', '>', 30),
    ]),
    allowedProperties: ['sex', 'created', 'age'],
    allowedOperators: {age: ['>', '<']},
    allowAdd: false,
  };
  const fb = mount(new FilterBuilder({schema: SCHEMA, template, value: Filters.applyTemplate(template)}));
  assert.equal(fb.root.querySelector('[data-u2-part="add"]').hidden, true);
  const [valueLocked, allLocked, free] = rows(fb);
  assert.equal(valueLocked.classList.contains('u2-fb-locked'), true);
  assert.equal(part(valueLocked, 'remove').hidden, true);
  assert.equal(part(valueLocked, 'prop').hidden, true, 'a locked picker gives way to text');
  assert.equal(part(valueLocked, 'op').hidden, true);
  assert.deepEqual([part(valueLocked, 'prop-text').textContent, part(valueLocked, 'op-text').textContent],
    ['sex', 'equals']);
  assert.equal(part(valueLocked, 'prop-text').hidden, false);
  assert.equal(part(free, 'prop-text').hidden, true);
  assert.equal(part(free, 'prop').hidden, false);
  assert.equal(select(valueLocked, 'value').disabled, false, 'the value stays editable');
  assert.equal(part(allLocked, 'remove').hidden, true);
  assert.equal(part(allLocked, 'value').querySelectorAll('input').every((i) => i.disabled), true);
  assert.equal(part(allLocked, 'value').dataset.u2, 'datetime-input');
  assert.equal(free.classList.contains('u2-fb-locked'), false);
  assert.equal(part(free, 'remove').hidden, false);
  assert.deepEqual(select(free, 'prop').querySelectorAll('option').map((o) => o.value), ['age', 'sex', 'created'],
    'allowedProperties narrows the picker, schema order kept');
  assert.deepEqual(select(free, 'op').querySelectorAll('option').map((o) => o.value), ['>', '<'],
    'allowedOperators narrows the operators');
  pick(select(valueLocked, 'value'), 'M');
  assert.deepEqual(Filters.validate(fb.value.peek(), SCHEMA, undefined, template), [], 'a value edit is allowed');
});

smoke('editors per kind: string, suggestions, int, float, bigint, choices, bool round-trip their values', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: tree(
    Filters.cond('name', '=', 'x'), Filters.cond('age', '=', 1), Filters.cond('mw', '=', 1.5),
    Filters.cond('big', '=', '12345678901234567890'), Filters.cond('sex', '=', 'F'),
    Filters.cond('active', '=', true))}));
  const [name, age, mw, big, sex, active] = rows(fb);
  assert.equal(part(name, 'value').dataset.u2, 'text-input');
  assert.equal(part(age, 'value').dataset.u2, 'number-input');
  assert.equal(part(mw, 'value').dataset.u2, 'number-input');
  assert.equal(part(big, 'value').dataset.u2, 'bigint-input');
  assert.equal(part(sex, 'value').dataset.u2, 'choice-input');
  assert.equal(part(active, 'value').dataset.u2, 'choice-input');
  type(part(name, 'value').querySelector('input'), 'y');
  type(part(age, 'value').querySelector('input'), '7');
  type(part(mw, 'value').querySelector('input'), '2.5');
  type(part(big, 'value').querySelector('input'), '99999999999999999999');
  pick(select(sex, 'value'), 'M');
  pick(select(active, 'value'), 'false');
  assert.deepEqual(fb.value.peek().nodes.map((n) => n.value), ['y', 7, 2.5, '99999999999999999999', 'M', false]);
  assert.deepEqual(fb.problems.peek(), []);
  type(part(name, 'value').querySelector('input'), '');
  assert.equal(fb.value.peek().nodes[0].value, undefined, 'a cleared editor clears the value');
  assert.deepEqual(fb.problems.peek().map((p) => p.code), ['missing-value']);
});

smoke('string over schema values gets a suggest input; the popup is closed with the builder', async () => {
  const fb = mount(new FilterBuilder({schema: WITH_VALUES, value: tree(Filters.cond('status', '=', ''))}));
  try {
    const editor = part(rows(fb)[0], 'value');
    assert.equal(editor.dataset.u2, 'suggest-input');
    type(editor.querySelector('input'), 'Op');
    await new Promise((resolve) => setTimeout(resolve, 200)); // the source's default debounce
    await flush();
    assert.equal(fb.value.peek().nodes[0].value, 'Op');
    assert.deepEqual(document.body.querySelectorAll('.u2-typeahead-option').map((r) => r.textContent), ['Open']);
  } finally {
    fb.dispose();
  }
  assert.equal(document.body.querySelector('.u2-typeahead-popup'), null);
});

smoke('the suggest editor opens the whole list on focus, on a click and on Space in an empty box', async () => {
  const fb = mount(new FilterBuilder({schema: WITH_VALUES, value: tree(Filters.cond('status', '=', 'Op'))}));
  const options = () => document.body.querySelectorAll('.u2-typeahead-option').map((r) => r.textContent);
  const settle = async () => {
    await new Promise((resolve) => setTimeout(resolve, 200));
    await flush();
  };
  try {
    const input = part(rows(fb)[0], 'value').querySelector('input');
    fire(input, 'focus');
    await settle();
    assert.deepEqual(options(), ['Open', 'Blocked', 'Closed'], 'focus over a held value: the unfiltered list');
    fire(input, 'keydown', {key: 'Escape'});
    assert.equal(options().length, 0);
    fire(input, 'click');
    await settle();
    assert.deepEqual(options(), ['Open', 'Blocked', 'Closed'], 'a click after Esc reopens');
    fire(input, 'keydown', {key: 'Escape'});
    type(input, '');
    await settle();
    fire(input, 'keydown', {key: 'Escape'});
    const typed = fire(input, 'keydown', {key: ' '});
    await settle();
    assert.equal(typed, false, 'the space is not typed');
    assert.deepEqual(options(), ['Open', 'Blocked', 'Closed'], 'Space in an empty box opens the list');
    assert.equal(fb.value.peek().nodes[0].value, undefined, 'the value stays empty');
    type(input, 'x');
    fire(input, 'keydown', {key: 'Escape'});
    assert.equal(fire(input, 'keydown', {key: ' '}), true, 'Space in a non-empty box types');
  } finally {
    fb.dispose();
  }
});

smoke('in → tags over the schema values; a typed value becomes a chip and a list value', async () => {
  const fb = mount(new FilterBuilder({schema: WITH_VALUES, value: tree(Filters.cond('status', 'in', ['Open']))}));
  try {
    const editor = part(rows(fb)[0], 'value');
    assert.equal(editor.dataset.u2, 'tags-input');
    assert.deepEqual(editor.querySelectorAll('.u2-tag-text').map((t) => t.textContent), ['Open']);
    const box = editor.querySelector('.u2-tags-input');
    type(box, 'Closed');
    fire(box, 'keydown', {key: 'Enter'});
    assert.deepEqual(fb.value.peek().nodes[0].value, ['Open', 'Closed']);
    assert.equal(fb.query.peek(), 'status in ("Open", "Closed")');
    fire(editor.querySelector('.u2-tag-remove'), 'click');
    assert.deepEqual(fb.value.peek().nodes[0].value, ['Closed']);
    fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: ['Blocked', 'Open']});
    assert.deepEqual(editor.querySelectorAll('.u2-tag-text').map((t) => t.textContent), ['Blocked', 'Open']);
  } finally {
    fb.dispose();
  }
});

smoke('between → two editors writing the bounds pair', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: tree(Filters.cond('age', 'between', [18, 65]))}));
  const row = rows(fb)[0];
  const lo = part(row, 'value').querySelector('input');
  const hi = part(row, 'value2').querySelector('input');
  assert.deepEqual([lo.value, hi.value], ['18', '65']);
  type(hi, '70');
  assert.deepEqual(fb.value.peek().nodes[0].value, [18, 70]);
  type(lo, '');
  assert.deepEqual(fb.value.peek().nodes[0].value, [null, 70]);
  assert.deepEqual(fb.problems.peek().map((p) => p.code), ['invalid-value']);
  fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: [20, 30]});
  assert.deepEqual([lo.value, hi.value], ['20', '30']);
});

smoke('datetime → one date-time box that takes an absolute date or a typed span', () => {
  const when = new Date(2026, 0, 15, 10, 30);
  const fb = mount(new FilterBuilder({schema: SCHEMA,
    value: tree(Filters.cond('created', '>', when), Filters.cond('created', '<', {span: '-1w'}))}));
  const [absolute, relative] = rows(fb);
  const box = part(absolute, 'value').querySelector('input');
  assert.equal(part(absolute, 'value').dataset.u2, 'datetime-input');
  assert.equal(box.value, '2026-01-15 10:30');
  assert.equal(box.placeholder, 'yyyy-MM-dd HH:mm or -1w');
  assert.equal(part(relative, 'value').querySelector('input').value, '-1w', 'a span shows as typed');

  type(box, '-2d');
  assert.deepEqual(fb.value.peek().nodes[0].value, {span: '-2d'});
  assert.equal(fb.query.peek(), 'created > -2d and created < -1w');
  type(box, 'soon');
  assert.deepEqual(fb.value.peek().nodes[0].value, {span: '-2d'}, 'unparseable text never reaches the model');
  assert.equal(part(absolute, 'value').querySelector('.u2-input-error').textContent, 'Not a date-time');
  type(box, '2026-02-01 08:00');
  assert.equal(fb.value.peek().nodes[0].value.getTime(), new Date(2026, 1, 1, 8).getTime());

  fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: {span: 'now'}});
  assert.equal(box.value, 'now', 'a span set from the model shows as its text');
  fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: when});
  assert.equal(box.value, '2026-01-15 10:30', 'a date set from the model shows as the date');
  fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: undefined});
  assert.equal(box.value, '');
});

smoke('a custom editor factory wins over the kind default and round-trips the value', () => {
  const calls = [];
  const editors = (prop, options) => {
    calls.push([prop.name, options.value]);
    if (prop.name !== 'name')
      return null;
    const input = new TextInput({...options, placeholder: 'custom'});
    input.root.dataset.custom = '1';
    return input;
  };
  const fb = mount(new FilterBuilder({schema: SCHEMA, editors,
    value: tree(Filters.cond('name', '=', 'a'), Filters.cond('age', '=', 1))}));
  const [name, age] = rows(fb);
  assert.deepEqual(calls, [['name', 'a'], ['age', 1]]);
  assert.equal(part(name, 'value').dataset.custom, '1');
  assert.equal(part(age, 'value').dataset.u2, 'number-input', 'null falls back to the kind default');
  type(part(name, 'value').querySelector('input'), 'b');
  assert.equal(fb.value.peek().nodes[0].value, 'b');
  fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: 'c'});
  assert.equal(part(name, 'value').querySelector('input').value, 'c');
});

smoke('problems follow validate; the target narrows what is expressible', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, target: 'dataframe',
    value: tree(Filters.cond('name', 'fuzzy', 'x'), Filters.cond('nope', '=', 1))}));
  assert.deepEqual(fb.problems.peek().map((p) => [p.nodeId, p.code]),
    [['f1', 'not-expressible'], ['f2', 'unknown-property']]);
  assert.deepEqual(rows(fb).map((r) => r.classList.contains('u2-fb-invalid')), [true, true]);
  assert.equal(select(rows(fb)[1], 'prop').value, 'nope', 'an unknown property is still shown');
});

smoke('a nested tree in simple mode: a summary row, the hint, flatten when the connectors allow', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: Filters.group('and', [
    Filters.cond('age', '>', 30),
    Filters.group('and', [Filters.cond('sex', '=', 'F')]),
  ])}));
  const hint = fb.root.querySelector('[data-u2-part="hint"]');
  assert.equal(hint.hidden, false);
  const nested = rows(fb)[1];
  assert.equal(nested.classList.contains('u2-fb-nested'), true);
  assert.equal(part(nested, 'value').textContent, '(sex = "F")');
  const flatten = hint.querySelector('[data-u2-part="flatten"]');
  assert.equal(flatten.hidden, false);
  const hintMode = hint.querySelector('[data-u2-part="hint-mode"]');
  assert.equal(hintMode.hidden, false, 'the hint offers advanced mode too');
  fire(hintMode, 'click');
  assert.equal(fb.mode.peek(), 'advanced');
  assert.equal(hint.hidden, true);
  assert.equal(fb.setMode('simple'), false);
  assert.equal(hintMode.hidden, true, 'already in advanced mode: only flatten is offered');
  assert.equal(flatten.hidden, false);
  fire(flatten, 'click');
  assert.equal(Filters.isFlat(fb.value.peek()), true);
  assert.equal(hint.hidden, true);
  assert.deepEqual(rows(fb).map((r) => r.dataset.u2Node), ['f1', 'f2'], 'the inlined condition keeps its id');
  assert.equal(rows(fb)[1].classList.contains('u2-fb-nested'), false, 'the condition got a real row');
});

smoke('setMode: advanced adopts, back to simple is refused with the hint while nested', () => {
  const mode = signal('advanced');
  const fb = mount(new FilterBuilder({schema: SCHEMA, mode, value: Filters.group('and', [
    Filters.cond('age', '>', 30),
    Filters.group('or', [Filters.cond('sex', '=', 'F'), Filters.cond('name', 'like', 'a')]),
  ])}));
  assert.equal(fb.mode === mode, true, 'a signal is adopted');
  const hint = fb.root.querySelector('[data-u2-part="hint"]');
  assert.equal(hint.hidden, true, 'no hint in advanced mode');
  assert.equal(fb.setMode('simple'), false);
  assert.equal(mode.peek(), 'advanced');
  assert.equal(hint.hidden, false);
  assert.equal(hint.querySelector('[data-u2-part="flatten"]').hidden, true, 'mixed connectors cannot flatten');
  fb.removeNode('f4');
  assert.equal(fb.setMode('simple'), true);
  assert.equal(mode.peek(), 'simple');
  assert.equal(hint.hidden, true);
  assert.equal(fb.addGroup().startsWith('f'), true);
  assert.equal(rows(fb)[1].classList.contains('u2-fb-nested'), true);
  assert.equal(fb.value.peek().nodes[1].op, 'or', 'the other connector');
});

smoke('orientation is ignored in advanced mode, with one warning', () => {
  const warnings = [];
  const original = console.warn;
  console.warn = (m) => warnings.push(m);
  try {
    const fb = mount(new FilterBuilder({schema: SCHEMA, orientation: 'horizontal', mode: 'advanced'}));
    fb.addCondition();
    fb.addCondition();
    assert.deepEqual(warnings, ['u2: FilterBuilder ignores orientation in advanced mode']);
  } finally {
    console.warn = original;
  }
});

smoke('a bound literal without ids is adopted with fresh ids (a system write)', () => {
  const bound = signal({op: 'and', nodes: [{property: 'age', operator: '>', value: 30}]});
  const fb = mount(new FilterBuilder({schema: SCHEMA, bind: bound}));
  assert.equal(typeof bound.peek().id, 'string');
  assert.equal(typeof bound.peek().nodes[0].id, 'string');
  assert.equal(rows(fb)[0].dataset.u2Node, bound.peek().nodes[0].id);
  fire(fb.root.querySelector('[data-u2-part="connector"] [data-id="or"]'), 'click');
  assert.equal(bound.peek().op, 'or', 'the root row addresses the adopted root, not the id it was built with');
  bound.value = undefined;
  assert.deepEqual(bound.peek().nodes, [], 'nothing is a fresh empty group');
});

smoke('rows come and go without leaving disposers behind on the builder', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, orientation: 'horizontal'}));
  const owned = fb.scope._disposers.length;
  for (let i = 0; i < 50; i++) {
    const a = fb.addCondition();
    const b = fb.addCondition();
    fb.removeNode(a);
    fb.removeNode(b);
  }
  assert.equal(rows(fb).length, 0);
  assert.equal(fb.scope._disposers.length, owned, 'a released row disowns its disposer');
});

smoke('a spec template without node ids seeds the value, locks its rows and never hangs the lock walk', () => {
  const reg = new Registry();
  registerAll(reg);
  const template = {root: {op: 'and', nodes: [
    {property: 'age', operator: '>', value: 18, lock: 'value'},
    {property: 'name', operator: 'like', value: 'a'},
  ]}, allowAdd: false};
  const schema = {properties: SCHEMA.properties};
  const instance = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-filter-builder', props: {schema, template}}},
    new SpecContext(), reg);
  document.body.append(instance.root);
  try {
    assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
    const fb = Control.forElement(instance.root.querySelector('[data-u2="filter-builder"]'));
    const [locked, free] = rows(fb);
    assert.equal(rows(fb).length, 2, 'the template is the starting value');
    assert.equal(typeof fb.value.peek().nodes[0].id, 'string');
    assert.equal(locked.classList.contains('u2-fb-locked'), true);
    assert.equal(part(locked, 'remove').hidden, true);
    assert.equal(free.classList.contains('u2-fb-locked'), false);
    assert.equal(fb.root.querySelector('[data-u2-part="add"]').hidden, true);
    assert.deepEqual(fb.problems.peek(), [], 'the seeded value is the template: nothing locked yet');
    type(part(locked, 'value').querySelector('input'), '21');
    assert.deepEqual(fb.problems.peek(), [], 'a value edit under a value lock is allowed');
  } finally {
    instance.dispose();
  }
  const both = renderSpec({$schema: 'dg-ui/1', root: {tag: 'u2-filter-builder', props: {schema, template,
    value: {op: 'and', nodes: [{property: 'age', operator: '>', value: 30}]}}}}, new SpecContext(), reg);
  document.body.append(both.root);
  try {
    assert.equal(both.root.querySelectorAll('[data-u2-node]').length, 1, 'an id-less value beside an id-less template');
  } finally {
    both.dispose();
  }
});

smoke('showQuery footer, getWidgetStatus and the onChanged hook', () => {
  const changes = [];
  const fb = mount(new FilterBuilder({schema: SCHEMA, showQuery: true,
    onChanged: (v) => changes.push(Filters.format(v)), value: tree(Filters.cond('age', '>', 30))}));
  const footer = fb.root.querySelector('[data-u2-part="query"]');
  assert.equal(footer.textContent, 'age > 30');
  fb.addCondition();
  assert.equal(footer.textContent, 'age > 30 and name =');
  const status = fb.getWidgetStatus();
  assert.equal(status.tree === fb.value.peek(), true);
  assert.equal(status.query, 'age > 30 and name =');
  assert.deepEqual(status.problems.map((p) => p.code), ['missing-value']);
  assert.equal(status.mode, 'simple');
  assert.deepEqual(Object.keys(status.parts).sort(), ['add', 'connector', 'editor', 'error', 'mode', 'query', 'rows']);
  assert.deepEqual(changes, ['age > 30 and name =']);
  fb.value.value = Filters.update(fb.value.peek(), 'f1', {value: 40});
  assert.deepEqual(changes, ['age > 30 and name =', 'age > 40 and name ='], 'the editor does not echo a model write');
  const plain = mount(new FilterBuilder({schema: SCHEMA, value: tree(Filters.cond('age', '>', 30))}));
  assert.equal(plain.root.querySelector('[data-u2-part="query"]'), null, 'the footer is off by default');
  assert.equal('query' in plain.getWidgetStatus().parts, false);
});

smoke('bindStep answers query (read-only) and mode (writable); bindProps advertises both', () => {
  const fb = mount(new FilterBuilder({schema: SCHEMA, value: tree(Filters.cond('age', '>', 30))}));
  fb.componentMeta = {tag: 'u2-filter-builder', props: [{name: 'value', type: TYPE.OBJECT},
    {name: 'query', type: TYPE.STRING}, {name: 'mode', type: TYPE.STRING}]};
  assert.equal(fb.bindStep('query').value, 'age > 30');
  assert.equal(fb.bindStep('mode') === fb.mode, true);
  assert.equal(fb.bindStep('') === fb.value, true);
  const props = Object.fromEntries(fb.bindProps().map((p) => [p.name, p.writable]));
  assert.deepEqual(props, {value: true, query: false, mode: true});
});

smoke('disposing the builder releases every row, editor and popup scope', () => {
  const base = Scope.liveCount;
  const fb = mount(new FilterBuilder({schema: WITH_VALUES, value: tree(
    Filters.cond('status', 'in', ['Open']), Filters.cond('status', 'like', 'x'))}));
  assert.equal(Scope.liveCount > base, true);
  fb.dispose();
  assert.equal(Scope.liveCount, base);
});
