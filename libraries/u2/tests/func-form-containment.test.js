import {test, testOwner} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {DataFrame, Func, FuncParam, strayTables} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);

const grok = await import('datagrok-api/grok');
await import('datagrok-api/dg');
const {funcForm} = await import('../src/dg/funcs/func-form.js');
const {Editors} = await import('../src/dg/forms/editors.js');

const shell = grok.shell;

function contained(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      Func.registry = [];
      shell.dart.tables = strayTables();
      shell.dart.t = null;
      shell.dart.tableAdded.subs.clear();
      shell.dart.tableRemoved.subs.clear();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const fp = (name, type, data = {}) => new FuncParam(name, type, data);
const callOf = (...params) => new Func('containment', {inputs: params}).prepare({});

function columnP(name, data = {}, type = 'column') {
  const p = fp(name, type, data);
  p.dart.parentTableParamName = 'df';
  return p;
}

function mount(f) {
  document.body.append(f.root);
  return f;
}

/** A platform editor that cannot be built — for a param of semType `Broken`. */
function breaking(message) {
  return Editors.register({match: (p) => p.semType === 'Broken', create: () => {
    throw new Error(message);
  }});
}

/** The form's rows as `name` (an input) or `!name` (an error row), in document order. */
function rows(f) {
  return [...f.root.querySelectorAll('[data-u2-name]')]
    .filter((el) => el.dataset.u2 === 'readonly-field' || el.classList.contains('u2-input-root'))
    .map((el) => el.dataset.u2 === 'readonly-field' ? `!${el.dataset.u2Name}` : el.dataset.u2Name);
}

contained('an editor that fails to build is an error row in its place and a problem; the other fields work', () => {
  const off = breaking('Cannot read properties of null (reading \'label\')');
  try {
    const call = callOf(fp('name', 'string'), fp('structure', 'string', {semType: 'Broken', friendlyName: 'Core'}),
      fp('count', 'int'));
    const f = mount(funcForm(call));
    assert.deepEqual(rows(f), ['name', '!structure', 'count']);
    const row = f.root.querySelector('[data-u2="readonly-field"]');
    assert.equal(row.querySelector('.u2-input-label').textContent, 'Core');
    assert.equal(row.querySelector('.u2-param-source-error').textContent,
      'Could not show Core: Cannot read properties of null (reading \'label\')');
    assert.deepEqual(f.problems.value, [{name: 'structure',
      message: 'Cannot read properties of null (reading \'label\')'}]);
    assert.deepEqual(f.inputs.map((i) => i.name), ['name', 'count'], 'the failed field is no input of the form');
    assert.equal(f.getInput('structure'), undefined);
    const count = f.getInput('count').root.querySelector('.u2-number-box input');
    count.value = '3';
    fire(count, 'input');
    fire(count, 'blur');
    assert.equal(call.inputParams.values().find((p) => p.name === 'count').value, 3, 'the others still write');
    assert.deepEqual(f.missingRequired.value, []);
    f.dispose();
  }
  finally {
    off();
  }
});

contained('a param whose value cannot be read is an error row too, and the form still builds', () => {
  const call = callOf(fp('a', 'string'), fp('b', 'int'), fp('c', 'bool'));
  Object.defineProperty(call.inputParams.values()[1], 'value', {get() {
    throw new Error('the value is a disposed handle');
  }});
  const f = mount(funcForm(call));
  assert.deepEqual(rows(f), ['a', '!b', 'c']);
  assert.deepEqual(f.problems.value.map((p) => p.name), ['b']);
  f.dispose();
});

contained('a failed field under a category stays in its section; a failed table orphans its column fields', () => {
  const off = breaking('no table picker');
  try {
    const call = callOf(fp('df', 'dataframe', {semType: 'Broken', category: 'Data'}),
      columnP('col', {category: 'Data'}), fp('n', 'int', {category: 'Options'}));
    const f = mount(funcForm(call));
    const sections = f.root.querySelectorAll('.u2-section');
    assert.equal(sections.length, 2);
    assert.equal(sections[0].querySelector('[data-u2="readonly-field"]').dataset.u2Name, 'df');
    assert.deepEqual(f.problems.value.map((p) => p.name), ['df']);
    assert.equal(f.getInput('col').disabled === null, false, 'a column field with no table field is orphaned, not broken');
    f.dispose();
  }
  finally {
    off();
  }
});

/** One param of every kind funcForm routes: scalars, editor hints, the dynamic routes, the table routes, files and
 * the unsupported ones. */
function everyKind() {
  return [
    fp('s', 'string'), fp('i', 'int'), fp('d', 'double'), fp('fl', 'float'), fp('nu', 'num'), fp('big', 'bigint'),
    fp('q', 'qnum'), fp('b', 'bool'), fp('dt', 'datetime'),
    fp('area', 'string', {options: {editor: 'textarea'}}), fp('pw', 'string', {options: {editor: 'password'}}),
    fp('sw', 'bool', {options: {editor: 'switch'}}), fp('sl', 'double', {options: {editor: 'slider'}, min: 0, max: 1}),
    fp('pick', 'string', {choices: ['x', 'y']}), fp('dyn', 'string', {options: {choices: 'Nowhere:choices'}}),
    fp('sug', 'string', {options: {suggestions: 'Nowhere:suggest'}}),
    fp('many', 'list', {options: {choices: 'Nowhere:choices'}}), fp('tags', 'list', {choices: ['p', 'q']}),
    fp('file', 'file'), fp('df', 'dataframe'), columnP('col'), columnP('num', {options: {type: 'numerical'}}),
    columnP('cols', {}, 'column_list'), fp('m', 'map'), fp('o', 'object'), fp('g', 'graphics'),
    fp('none', 'string', {options: {editor: 'none'}}),
  ];
}

contained('over a param of every kind, with the hostile workspace open, nothing throws and nothing is a problem', async () => {
  const hidden = shell.addTable(new DataFrame([{name: 'Molecule', type: 'string', semType: 'Molecule'},
    {name: '~height.cache', type: 'int'}], [{Molecule: 'CCO', '~height.cache': 1}], 'hostile'));
  hidden.name = null;
  shell.addTable(new DataFrame([{name: 'age', type: 'int'}], [{age: 3}], 'demog'));
  const params = everyKind();
  for (const options of [{}, {editors: 'auto'}, {showTableSelectors: false, skipTableAutoFill: true},
    {skipDefaultInit: true, condensed: true}]) {
    const f = mount(testOwner().build(() => funcForm(callOf(...params), options)));
    await f.settled;
    assert.deepEqual(f.problems.value, [], JSON.stringify(options));
    assert.ok(f.inputs.length > 20, `${f.inputs.length} fields`);
    f.dispose();
  }
});

contained('the table field holds the call\'s own frame: of two open tables under one name, the one in the call', () => {
  const first = shell.addTable(new DataFrame([{name: 'age', type: 'int'}], [{age: 1}], 'demog'));
  const second = shell.addTable(new DataFrame([{name: 'age', type: 'int'}], [{age: 2}], 'other'));
  second.name = 'demog';
  const call = callOf(fp('df', 'dataframe'));
  call.setParamValue('df', second);
  const f = mount(funcForm(call));
  const input = f.getInput('df');
  const select = input.root.querySelector('select');
  assert.equal(select.value, 'demog');
  assert.equal(input.table, second, 'the field holds the frame the call holds, not the first of that name');
  select.value = '';
  fire(select, 'change');
  assert.equal(call.inputParams.values()[0].value, null);
  select.value = 'demog';
  fire(select, 'change');
  assert.equal(call.inputParams.values()[0].value, second, 'picked again, the call gets the frame back');
  assert.notEqual(first, second);
  f.dispose();
});

contained('a column field whose table lists no matching column says so, a hidden `~` column not counted', () => {
  const table = shell.addTable(new DataFrame([{name: 'sex', type: 'string'}, {name: '~weight.cache', type: 'int'}],
    [{sex: 'M', '~weight.cache': 1}], 'hidden'));
  shell.t = table;
  const f = mount(funcForm(callOf(fp('df', 'dataframe'), columnP('num', {nullable: false, options: {type: 'numerical'}}))));
  assert.equal(f.getInput('num').validity.peek(), 'No numerical columns in \'hidden\'');
  f.dispose();
});
