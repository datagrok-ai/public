/* FilterQueryInput: the completion context at every expectation, picking with quoting and the
   re-open, commit into the tree, problems kept out of the tree, re-formatting from the bound
   signal, two controls on one signal, and the registration. Every mounted control is disposed
   in `finally` — an overlay left open loops the animation-frame updater and OOMs node --test. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope, signal, Filters, FilterBuilder, FilterQueryInput, Control} from '../src/index.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import {TYPE} from 'datagrok-api/u2core';

const PROPS = [
  {name: 'name', type: TYPE.STRING, friendlyName: 'Name'},
  {name: 'age', type: TYPE.INT, min: 0, max: 120},
  {name: 'sex', type: TYPE.STRING, choices: ['F', 'M']},
  {name: 'active', type: TYPE.BOOL},
  {name: 'created', type: TYPE.DATE_TIME},
  {name: 'mw', type: TYPE.FLOAT},
];
const VALUES = {name: ['Aspirin', 'Ibuprofen', 'Naproxen'], status: ['Open', 'Blocked']};
const SCHEMA = Filters.schema(PROPS, VALUES);

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
    assert.equal(document.body.querySelector('.u2-fq-popup'), null, 'popup gone with the control');
  });
}

function mount(component) {
  document.body.append(component.root);
  mounted.push(component);
  return component;
}

function box(q) {
  return q.root.querySelector('input');
}

/** Types the whole text — the shim has no caret, so completion runs at the end. */
async function type(q, text) {
  const input = box(q);
  input.focus();
  input.value = text;
  fire(input, 'input');
  await flush();
}

function rows() {
  return [...document.body.querySelectorAll('.u2-fq-option')].map((r) => r.textContent);
}

function key(q, name) {
  fire(box(q), 'keydown', {key: name});
}

smoke('renders as a combobox editor showing the formatted tree', () => {
  const q = mount(new FilterQueryInput({label: 'Query', schema: SCHEMA, name: 'q',
    value: Filters.group('and', [Filters.cond('age', '>', 30), Filters.cond('sex', '=', 'F')])}));
  assert.equal(q.root.dataset.u2, 'filter-query-input');
  assert.equal(q.root.dataset.u2Name, 'q');
  const input = box(q);
  assert.equal(input.dataset.u2Part, 'editor');
  assert.equal(input.getAttribute('role'), 'combobox');
  assert.equal(input.value, 'age > 30 and sex = "F"');
  assert.equal(q.text.value, 'age > 30 and sex = "F"');
  assert.equal(q.query.value, 'age > 30 and sex = "F"');
  assert.equal(q.isOpen.value, false);
  const status = q.getWidgetStatus();
  assert.deepEqual([status.mode, status.query, status.problems], ['text', 'age > 30 and sex = "F"', []]);
  assert.equal(status.tree, q.value.peek());
});

smoke('completion: properties with kind icons, then operators, then values, then connectors', async () => {
  const q = mount(new FilterQueryInput({schema: SCHEMA}));
  await type(q, '');
  assert.equal(q.isOpen.value, true);
  assert.deepEqual(rows(), ['Name' + 'name', 'age', 'sex', 'active', 'created', 'mw']);
  const glyph = document.body.querySelector('.u2-fq-option .u2-fq-icon');
  assert.equal(glyph.classList.contains('fa-font'), true, 'string kind → font glyph');
  assert.equal(box(q).getAttribute('aria-expanded'), 'true');

  await type(q, 'a');
  assert.deepEqual(rows(), ['age', 'active']);
  assert.equal(q.root.querySelector('input').getAttribute('aria-activedescendant') !== null, true,
    'a typed prefix auto-highlights the first row');

  await type(q, 'age ');
  const ops = rows();
  assert.equal(ops[0], '=' + 'equals');
  assert.equal(ops.includes('between' + 'between'), true);
  assert.equal(ops.includes('= null' + 'is empty'), true, 'the model-only operator is spelled the grammar way');
  assert.equal(ops.some((r) => r.startsWith('like')), false, 'int has no text operators');

  await type(q, 'age > ');
  assert.deepEqual(rows(), ['number 0–120', 'null'], 'a bounded number: its range as an info row, then null');
  assert.equal(q.isOpen.value, true);
  key(q, 'ArrowDown');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, 'age > ', 'the info row inserts nothing');
  await type(q, 'age > 3');
  assert.deepEqual(rows(), ['number 0–120'], 'the range stays while a number is typed');
  await type(q, 'mw > ');
  assert.equal(q.isOpen.value, false, 'a lone null is not worth a popup');
  key(q, 'ArrowDown');
  await flush();
  assert.equal(q.isOpen.value, false, 'not on request either');
  await type(q, 'active = ');
  assert.deepEqual(rows(), ['true', 'false', 'null']);
  await type(q, 'created > -');
  assert.deepEqual(rows(), ['-1d', '-1w', '-1m']);

  await type(q, 'name = ');
  assert.deepEqual(rows(), ['Aspirin', 'Ibuprofen', 'Naproxen', 'null'], 'schema values, then the hint');
  await type(q, 'name = "ib');
  assert.deepEqual(rows(), ['Ibuprofen']);

  await type(q, 'name = "Ibuprofen" ');
  assert.deepEqual(rows(), ['and', 'or']);
  await type(q, 'name = "Ibuprofen" o');
  assert.deepEqual(rows(), ['or']);

  await type(q, 'AG');
  assert.deepEqual(rows(), ['age'], 'a property prefix matches regardless of case');
  await type(q, 'AGE > ');
  assert.deepEqual(rows(), ['number 0–120', 'null'], 'a property typed in another case still drives the context');
  key(q, 'Enter');
  await type(q, 'AGE > 3');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, 'age > 3', 'the commit re-formats with the schema spelling');
  assert.equal(q.value.value.nodes[0].property, 'age');
});

smoke('a pick replaces the token, quotes strings, appends a space and re-opens for the next expectation', async () => {
  const q = mount(new FilterQueryInput({schema: SCHEMA}));
  await type(q, 'na');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, 'name ');
  assert.equal(q.isOpen.value, true, 're-opened');
  assert.equal(rows()[0], '=' + 'equals');
  // an empty prefix highlights nothing, so the first ArrowDown lands on row 0
  for (let i = rows().findIndex((r) => r.startsWith('like')); i >= 0; i--)
    key(q, 'ArrowDown');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, 'name like ');
  await type(q, 'name like nap');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, 'name like "Naproxen" ', 'a string value is quoted');
  assert.deepEqual(rows(), ['and', 'or']);
  key(q, 'Enter');
  assert.equal(q.isOpen.value, false, 'Enter with no active row commits');
  assert.equal(q.query.value, 'name like "Naproxen"', 'the committed text is re-formatted');
  assert.deepEqual(q.value.peek().nodes.map((n) => [n.property, n.operator, n.value]),
    [['name', 'like', 'Naproxen']]);

  await type(q, 'name like "Naproxen" and act');
  fire(document.body.querySelector('.u2-fq-option'), 'pointerdown');
  await flush();
  assert.equal(box(q).value, 'name like "Naproxen" and active ', 'a click picks too');
  assert.equal(q.isOpen.value, true);
  key(q, 'Escape');
  assert.equal(q.isOpen.value, false, 'Escape closes');
});

smoke('commit parses the text into the tree; a problem keeps the old tree and surfaces the message', async () => {
  const q = mount(new FilterQueryInput({schema: SCHEMA, value: Filters.group('and', [Filters.cond('age', '>', 30)])}));
  const before = q.value.peek();
  q.text.value = 'age > 5 and (sex = "F" or name like "an")';
  assert.equal(q.commit(), true);
  const root = q.value.peek();
  assert.notEqual(root, before);
  assert.equal(Filters.format(root), 'age > 5 and (sex = "F" or name like "an")');
  assert.equal(root.nodes.length, 2);
  assert.equal(Filters.isGroup(root.nodes[1]), true);
  assert.equal(q.validity.value, null);

  q.text.value = 'age >';
  assert.equal(q.commit(), false);
  assert.equal(q.value.peek(), root, 'the tree is untouched');
  assert.equal(q.problems.value[0].code, 'syntax');
  assert.match(q.validity.value, /^Expected /);
  assert.equal(box(q).classList.contains('u2-invalid'), true);

  q.text.value = 'height > 5';
  assert.equal(q.commit(), false);
  assert.equal(q.problems.value[0].code, 'unknown-property');
  assert.equal(q.value.peek(), root);

  q.text.value = '';
  assert.equal(q.commit(), true);
  assert.deepEqual(q.value.peek().nodes, [], 'blank text is the empty tree');
  assert.equal(q.validity.value, null);
});

smoke('blur commits a dirty text; a tree written elsewhere re-formats the text unless the box is focused ' +
  'and dirty', async () => {
  const tree = signal(Filters.group('and', [Filters.cond('age', '>', 30)]));
  const q = mount(new FilterQueryInput({schema: SCHEMA, bind: tree}));
  assert.equal(box(q).value, 'age > 30');

  tree.value = Filters.group('or', [Filters.cond('sex', '=', 'M'), Filters.cond('name', 'starts', 'A')]);
  assert.equal(box(q).value, 'sex = "M" or name starts "A"');

  await type(q, 'age >= 18');
  tree.value = Filters.group('and', [Filters.cond('active', '=', true)]);
  assert.equal(box(q).value, 'age >= 18', 'the draft under the user survives');
  box(q).blur();
  assert.equal(Filters.format(tree.value), 'age >= 18', 'blur commits the draft');
  assert.equal(q.isOpen.value, false);

  await type(q, 'age >=');
  box(q).blur();
  assert.equal(Filters.format(tree.value), 'age >= 18', 'a broken draft leaves the tree alone');
  assert.equal(q.problems.value[0].code, 'syntax');
  tree.value = Filters.group('and', [Filters.cond('active', '=', true)]);
  assert.equal(box(q).value, 'active = true', 'unfocused: the tree wins');
  assert.deepEqual(q.problems.value, [], 'and the problems of the last draft go with it');

  await type(q, 'active = true');
  tree.value = Filters.group('and', [Filters.cond('age', '<', 5)]);
  assert.equal(box(q).value, 'active = true', 'a typed text equal to the tree is still the user\'s draft');
  fire(box(q), 'focus');
  tree.value = Filters.group('and', [Filters.cond('age', '<', 9)]);
  assert.equal(box(q).value, 'active = true', 'focus with a text the tree does not spell keeps the draft');
  box(q).blur();
  q.text.value = Filters.format(tree.peek());
  fire(box(q), 'focus');
  tree.value = Filters.group('and', [Filters.cond('age', '<', 7)]);
  assert.equal(box(q).value, 'age < 7', 'focus over a text equal to the tree is not a draft: the tree wins');
});

smoke('a caret moved by the mouse refreshes the context: the pick lands under the caret, not at the end', async () => {
  const q = mount(new FilterQueryInput({schema: SCHEMA}));
  await type(q, 'age > 30 and na');
  assert.deepEqual(rows(), ['Name' + 'name'], 'the end of the text: a property prefix');
  const input = box(q);
  input.selectionStart = 1;
  fire(input, 'click');
  await flush();
  assert.deepEqual(rows(), ['age', 'active'], 'the click moved the context to the first token');
  key(q, 'ArrowDown');
  key(q, 'Enter');
  await flush();
  assert.equal(input.value, 'active > 30 and na', 'replaced the token under the caret, one space kept');
  input.selectionStart = 12;
  fire(input, 'select');
  await flush();
  assert.deepEqual(rows(), ['and', 'or'], 'select refreshes too');
  key(q, 'Escape');
  fire(input, 'click');
  assert.equal(q.isOpen.value, false, 'a click with the list closed opens nothing but at a value position');
});

smoke('a click where a value is expected opens the values; Space in an empty box opens the properties', async () => {
  const q = mount(new FilterQueryInput({schema: SCHEMA}));
  await type(q, 'name = ');
  key(q, 'Escape');
  assert.equal(q.isOpen.value, false);
  const input = box(q);
  fire(input, 'click');
  await flush();
  assert.deepEqual(rows(), ['Aspirin', 'Ibuprofen', 'Naproxen', 'null'], 'the click at the value position offers the values');
  key(q, 'Escape');
  input.selectionStart = 2;
  fire(input, 'click');
  assert.equal(q.isOpen.value, false, 'a click on the property token opens nothing');
  await type(q, '');
  key(q, 'Escape');
  const typed = fire(input, 'keydown', {key: ' '});
  await flush();
  assert.equal(typed, false, 'the space is not typed');
  assert.equal(rows().length, PROPS.length, 'Space in an empty box offers the properties');
  assert.equal(input.value, '');
  key(q, 'Escape');
  await type(q, 'a');
  key(q, 'Escape');
  assert.equal(fire(input, 'keydown', {key: ' '}), true, 'Space in a non-empty box types');
});

smoke('a commit of text equal to the tree re-formats without a new root: a bound builder keeps its rows', () => {
  const tree = signal(Filters.group('and', [Filters.cond('age', '>', 30)]));
  const fb = mount(new FilterBuilder({schema: SCHEMA, bind: tree}));
  const q = mount(new FilterQueryInput({schema: SCHEMA, bind: tree}));
  const root = tree.peek();
  const row = fb.root.querySelector('[data-u2-node]');
  q.text.value = '  age   >  30 ';
  assert.equal(q.commit(), true);
  assert.equal(tree.peek(), root, 'the same root');
  assert.equal(fb.root.querySelector('[data-u2-node]'), row, 'the same row element');
  assert.equal(box(q).value, 'age > 30', 're-formatted');
  assert.deepEqual(q.problems.value, []);
});

smoke('a value picked under a like-family operator is spelled as a literal pattern and round-trips', async () => {
  const schema = Filters.schema(PROPS, {name: ['50%', 'in_progress']});
  const q = mount(new FilterQueryInput({schema}));
  await type(q, 'name like ');
  assert.deepEqual(rows(), ['50%', 'in_progress', 'null']);
  key(q, 'ArrowDown');
  key(q, 'Enter');
  await flush();
  const spelled = (op, v) => Filters.format(Filters.group('and', [Filters.cond('name', op, v)]));
  assert.equal(box(q).value, `${spelled('like', '50%')} `, 'escaped the way the formatter spells it');
  key(q, 'Escape');
  key(q, 'Enter');
  const [cond] = q.value.peek().nodes;
  assert.deepEqual([cond.operator, cond.value, cond.options], ['like', '50%', undefined], 'a literal, not a wildcard');
  assert.equal(q.query.value, spelled('like', '50%'));
  await type(q, 'name starts in');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, `${spelled('starts', 'in_progress')} `);
  await type(q, 'name = 50');
  key(q, 'Enter');
  await flush();
  assert.equal(box(q).value, 'name = "50%" ', 'equality takes the value as it is');
});

smoke('a builder and a query input on one signal stay in sync both ways', async () => {
  const tree = signal(Filters.group('and', [Filters.cond('age', '>', 30)]));
  const fb = mount(new FilterBuilder({schema: SCHEMA, bind: tree}));
  const q = mount(new FilterQueryInput({schema: SCHEMA, bind: tree}));
  const op = fb.root.querySelector('[data-u2-part="op"] select');
  op.value = '<';
  fire(op, 'change');
  assert.equal(box(q).value, 'age < 30', 'builder edit → text');

  q.text.value = 'age < 30 and sex = "F"';
  assert.equal(q.commit(), true);
  assert.equal(fb.root.querySelectorAll('[data-u2-node]').length, 2, 'text commit → builder rows');
  assert.equal(fb.query.value, 'age < 30 and sex = "F"');
  assert.equal(fb.value, tree, 'both adopted the signal');
  assert.equal(q.value, tree);
});

smoke('a rejected value provider shows the error row; an aborted one is ignored', async () => {
  let reject;
  const schema = {properties: PROPS, values: () => new Promise((_, r) => reject = r)};
  const q = mount(new FilterQueryInput({schema}));
  await type(q, 'name = ');
  assert.notEqual(document.body.querySelector('.u2-fq-loading'), null);
  reject(new Error('offline'));
  await flush();
  assert.equal(document.body.querySelector('.u2-fq-error').textContent, 'offlineRetry');
  await type(q, 'age = ');
  assert.equal(document.body.querySelector('.u2-fq-error'), null, 'a new context replaces the failed one');
  assert.notEqual(document.body.querySelector('.u2-fq-loading'), null);
});

function spec(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      Filters.resetIds('');
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function registry() {
  const reg = new Registry();
  registerAll(reg);
  return reg;
}

spec('registration: the manifest entry states the grammar and the binding tiers; the example renders', () => {
  const reg = registry();
  const meta = reg.manifest().components.find((c) => c.tag === 'u2-filter-query-input');
  assert.equal(meta.category, 'Inputs');
  assert.match(meta.usage, /prop op value/);
  assert.match(meta.usage, /age > 30 and sex = "F"/);
  assert.match(meta.usage, /mw between 200 and 500/);
  assert.match(meta.usage, /-1w 2d now/);
  const props = Object.fromEntries(meta.props.map((p) => [p.name, p]));
  assert.deepEqual([props.value.type, props.value.bindable, props.value.twoWay], [TYPE.OBJECT, true, true]);
  assert.deepEqual([props.query.type, props.query.bindable, props.query.twoWay], [TYPE.STRING, true, undefined]);
  assert.equal(props.schema.bindable, undefined, 'schema is re-render tier');
  assert.deepEqual(meta.events, ['change']);

  const instance = renderSpec({$schema: 'dg-ui/1', root: reg.get('u2-filter-query-input').example},
    new SpecContext(), reg);
  document.body.append(instance.root);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  assert.equal(instance.root.querySelector('input').value, 'age > 30');
  assert.equal(instance.root.querySelector('input').placeholder, 'age > 30 and sex = "F"');
  instance.dispose();
});

spec('registration: value binds two-way as the context signal itself; query reads through a bind path', () => {
  const reg = registry();
  const ctx = new SpecContext({data: {criteria: Filters.group('and', [Filters.cond('age', '>', 30)])}});
  const instance = renderSpec({
    $schema: 'dg-ui/1',
    root: {tag: 'u2-div-v', children: [
      {tag: 'u2-filter-query-input', name: 'q', props: {schema: {properties: PROPS, values: VALUES}},
        bind: {value: '$.criteria'}},
      {tag: 'h3', name: 'echo', bind: {text: '$.q.query'}},
    ]},
  }, ctx, reg);
  document.body.append(instance.root);
  assert.equal(instance.root.querySelectorAll('.u2-spec-error').length, 0);
  const q = Control.forElement(instance.root.querySelector('[data-u2="filter-query-input"]'));
  assert.equal(q.value === ctx.data.criteria, true, 'the context signal itself');
  const echo = instance.root.querySelector('h3');
  assert.equal(echo.textContent, 'age > 30');
  q.text.value = 'age > 30 and sex = "M"';
  q.commit();
  assert.equal(ctx.data.criteria.peek().nodes.length, 2, 'commit reaches the context');
  assert.equal(echo.textContent, 'age > 30 and sex = "M"');
  ctx.data.criteria.value = Filters.group('and', [Filters.cond('name', 'like', 'a')]);
  assert.equal(echo.textContent, 'name like "a"', 'a context write re-formats the text');
  const {writable} = instance.resolveBinding('$.q.query');
  assert.equal(writable, false);
  instance.dispose();
});
