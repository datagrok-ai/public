/* The W3 dataframe/column/column_list contract of plan-w3.md (FP-W3-2…FP-W3-9) over the FuncCall
   and Shell doubles: routing with association demotion, table auto-fill, objects-into-the-call /
   names-in-the-input conversion, dependent retargeting on table change, column auto-pick,
   ColumnFilter parity, the closed-table prune (#11) and the W2 source interplay. Written to the
   plan's contracts, blind to the implementation (the zero-reconciliation method).
   Pickers-parity (plan-pickers.md FP-P-4/FP-P-5): the column fields are the grid combo, so
   display-side reads go through ColumnInput.names() and the native popup — the call-side
   assertions stand verbatim as W3's contract. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Column, DataFrame, Func, FuncParam, Shell} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);

const grok = await import('datagrok-api/grok');
// loading the DG stub installs the feature-detected globals (grok_Property_Get among them)
await import('datagrok-api/dg');
const {funcForm} = await import('../src/dg/funcs/func-form.js');
const {TextInput} = await import('../src/components/inputs/text-input.js');

const shell = grok.shell;

/** The w2 harness plus shell hygiene: every test leaves no open tables, no current table, no
 * follow subscriptions and the live-scope count where it was. */
function w3(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      Func.registry = [];
      shell.dart.tables = [];
      shell.dart.t = null;
      shell.dart.tableAdded.subs.clear();
      shell.dart.tableRemoved.subs.clear();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** Doubles-only pins. */
function dbl(name, body) {
  test(name, async () => {
    try {
      await body();
    } finally {
      Func.registry = [];
      resetDom();
      await flush();
    }
  });
}

function fp(name, propertyType, data = {}) {
  return new FuncParam(name, propertyType, data);
}

function callOf(...params) {
  return new Func('fceTables', {inputs: params}).prepare({});
}

function tableP(name = 'df', data = {}) {
  return fp(name, 'dataframe', data);
}

/** The annotation linker's field, set the way the platform sets it (P-W3-1): directly on the
 * property handle, where `grok_Property_Get` answers it. */
function columnP(name, data = {}, parent = 'df') {
  const p = fp(name, 'column', data);
  p.dart.parentTableParamName = parent;
  return p;
}

function columnsP(name, data = {}, parent = 'df') {
  const p = fp(name, 'column_list', data);
  p.dart.parentTableParamName = parent;
  return p;
}

function demog(name = 'demog') {
  return new DataFrame([
    {name: 'age', type: 'int'}, {name: 'height', type: 'int'},
    {name: 'weight', type: 'int'}, {name: 'sex', type: 'string'},
  ], [{age: 25, height: 175, weight: 68, sex: 'M'}, {age: 30, height: 162, weight: 55, sex: 'F'}],
  name);
}

function alt(name = 'alt') {
  return new DataFrame([{name: 'weight', type: 'int'}, {name: 'label', type: 'string'}],
    [{weight: 70, label: 'a'}, {weight: 80, label: 'b'}], name);
}

function open(frame) {
  return shell.addTable(frame);
}

/** Counts every write, so the one-write-per-column pins read a plain ledger. */
function wired(call, options = {}) {
  const sets = [];
  const setParamValue = call.setParamValue.bind(call);
  call.setParamValue = (name, value) => {
    sets.push([name, value]);
    setParamValue(name, value);
  };
  const f = mount(funcForm(call, options));
  return {f, sets};
}

function mount(f) {
  document.body.append(f.root);
  return f;
}

function editorOf(f, name, selector = 'input') {
  return f.getInput(name).root.querySelector(selector);
}

function selectOf(f, name) {
  const select = editorOf(f, name, 'select');
  assert.ok(select, `${name} renders a select`);
  return select;
}

/** The items on offer, the nullable empty option dropped. */
function items(f, name) {
  return selectOf(f, name).children.map((o) => o.value).filter((v) => v !== '');
}

function summaryOf(f, name) {
  return f.getInput(name).root.querySelector('.u2-columns-summary')?.textContent;
}

/** What a column-combo field offers — the ColumnInput answers its own filtered names (FP-P-4). */
function colNames(f, name) {
  const input = f.getInput(name);
  assert.equal(typeof input.names, 'function', `${name} renders the column combo`);
  return input.names();
}

/** Picks a column through the combo's own native popup — the user path (FP-P-4). */
async function pickColumnVia(f, name, colName) {
  fire(f.getInput(name).root.querySelector('.u2-columns'), 'click');
  const list = document.body.querySelector('.u2-overlay .u2-list');
  assert.ok(list, `${name} opens the column popup`);
  // the shim lays nothing out: give the list a viewport — and rewind the highlight scroll,
  // which the layoutless shim never clamps — so every row is rendered
  list.clientHeight = 400;
  list.scrollTop = 0;
  fire(list, 'scroll');
  const row = list.querySelectorAll('.u2-list-row')
    .find((el) => el.querySelector('.u2-column-name')?.textContent === colName);
  assert.ok(row, `${colName} is on offer`);
  fire(row, 'click');
  await flush();
}

function pickTable(f, name, value) {
  const select = selectOf(f, name);
  select.value = value ?? '';
  fire(select, 'change');
}

function provider(name, inputs, run) {
  const fn = new Func(name, {inputs, run});
  Func.registry.push(fn);
  return fn;
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function until(cond, what, timeoutMs = 1500) {
  const start = Date.now();
  while (!cond()) {
    if (Date.now() - start > timeoutMs)
      assert.fail(`timed out waiting for ${what}`);
    await flush();
  }
}

// --- the doubles themselves (FP-W3-9; green with or without the implementation) ---

dbl('doubles: the shell stores frames — tables, t, derived tableNames, dedup, close by frame or name', () => {
  const sh = new Shell();
  assert.deepEqual(Object.keys({...sh}), ['dart'], 'spread yields dart only');
  const events = [];
  sh.dart.tableAdded.subscribe(() => events.push('+'));
  sh.dart.tableRemoved.subscribe(() => events.push('-'));

  const a = sh.addTable(demog());
  assert.equal(sh.tables[0], a, 'the frame itself is stored');
  assert.deepEqual(sh.tableNames, ['demog'], 'names derive from the frames');
  assert.equal(sh.tableByName('demog'), a, 'tableByName answers the frame');
  assert.equal(sh.tableByName('nope'), null);
  const a2 = sh.addTable(demog());
  assert.equal(a2.name, 'demog (2)', 'name dedup kept');
  assert.equal(sh.t, null);
  sh.t = a2;
  assert.equal(sh.t, a2);

  sh.closeTable(a);
  assert.deepEqual(sh.tableNames, ['demog (2)'], 'close by frame');
  sh.closeTable('demog (2)');
  assert.deepEqual(sh.tables, [], 'close by name');
  assert.equal(sh.t, null, 'the current table never answers a closed frame');
  assert.deepEqual(events, ['+', '+', '-', '-']);

  sh.addTable(demog());
  sh.dart.tableNames = [];
  assert.deepEqual(sh.tables, [], 'the resetShell name-list write clears the same backing store');
});

dbl('doubles: Column isCategorical/categories and FuncParam columnTypeFilter, all getter-backed', () => {
  const frame = demog();
  const sex = frame.columns.byName('sex');
  assert.equal(sex.isCategorical, true);
  assert.equal(frame.columns.byName('age').isCategorical, false, 'int is not categorical');
  assert.deepEqual(sex.categories, ['M', 'F'], 'distinct values off the frame rows');
  assert.deepEqual(new Column('loose', 'string').categories, [], 'no frame, no categories');
  assert.deepEqual(Object.keys({...sex}), ['dart']);

  const p = fp('c', 'column', {options: {type: 'numerical'}});
  assert.equal(p.columnTypeFilter, 'numerical');
  assert.equal(fp('c2', 'column', {options: {columns: 'categorical'}}).columnTypeFilter,
    'categorical', 'the info ?? field union reads the columns key too');
  assert.equal(fp('c3', 'column').columnTypeFilter, null);
  assert.deepEqual(Object.keys({...p}), ['dart']);

  p.dart.parentTableParamName = 'df';
  assert.equal(globalThis.grok_Property_Get(p.dart, 'parentTableParamName'), 'df',
    'the generic metadata door the dg-stub installs');
  assert.equal(globalThis.grok_Property_Get(p.dart, 'maxCategories'), null, 'absent keys answer null');
});

// --- routing (FP-W3-2) ---

w3('routing: dataframe, column and column_list route to fields, in declaration order (#10)', () => {
  open(demog());
  const call = callOf(tableP(), columnP('age'), fp('note', 'string'), columnsP('cols'));
  const f = mount(funcForm(call));
  assert.deepEqual(f.inputs.map((i) => i.name), ['df', 'age', 'note', 'cols'],
    'declared order — never table-clustered as the Dart form renders (P-W3-3)');
  assert.deepEqual([...f.unsupported], []);
  assert.equal(f.getInput('df').root.dataset.u2, 'choice-input', 'the table field is tableInput');
  assert.equal(f.getInput('age').root.dataset.u2, 'column-combo',
    'the column field is the grid combo (FP-P-4/FP-P-5)');
  assert.equal(f.getInput('cols').root.dataset.u2, 'columns-input');
  f.dispose();
});

w3('routing: an editor-carrying dataframe/column/column_list param stays unsupported (W5 territory)', () => {
  open(demog());
  const call = callOf(
    tableP(),
    tableP('t', {options: {editor: 'grid'}}),
    columnP('c', {options: {editor: 'Chem:col'}}),
    columnsP('cs', {options: {editor: 'custom'}}));
  const f = mount(funcForm(call));
  assert.deepEqual([...f.unsupported], ['t', 'c', 'cs'],
    'in Dart a non-null editor wins over the type branches (fpe:497)');
  assert.ok(f.getInput('df'), 'the editor-less table param still routes');
  f.dispose();
});

w3('routing: editor: columnsMap (type map) stays unsupported — platform defect #9', () => {
  open(demog());
  const call = callOf(tableP(), fp('m', 'map', {options: {editor: 'columnsMap'}}));
  const f = mount(funcForm(call));
  assert.deepEqual([...f.unsupported], ['m'], 'the shipped Dart form crashes here (P-W3-7)');
  assert.equal(f.getInput('m'), undefined);
  f.dispose();
});

w3('routing: an unassociated column is LISTED unsupported, never silently dropped (#13)', () => {
  open(demog());
  const bare = fp('bare', 'column');
  const call = callOf(
    tableP(),
    columnP('good'),
    bare,
    columnP('ghost', {}, 'nope'),
    columnP('wrong', {}, 'label'),
    fp('label', 'string'));
  const f = mount(funcForm(call));
  assert.ok(f.getInput('good'), 'the linked column renders');
  assert.deepEqual([...f.unsupported], ['bare', 'ghost', 'wrong'],
    'no link, a missing param, and a non-dataframe param all demote — visibly');
  assert.equal(f.getInput('bare'), undefined);
  assert.equal(f.getInput('wrong'), undefined);
  f.dispose();
});

w3('routing: the explicit options.table fallback resolves without the metadata key', () => {
  const frame = open(demog());
  shell.t = frame;
  const call = callOf(tableP(), fp('cat', 'column', {options: {table: 'df', type: 'categorical'}}));
  const f = mount(funcForm(call));
  assert.deepEqual([...f.unsupported], [], 'the signature-registered decoration path (P-W3-1)');
  assert.deepEqual(colNames(f, 'cat'), ['sex']);
  assert.equal(call.inputs.cat, frame.columns.byName('sex'), 'auto-picked through the fallback link');
  f.dispose();
});

w3('routing: column_filter, row_filter, dataframe_list and list<column> stay unsupported', () => {
  open(demog());
  const call = callOf(tableP(),
    fp('cf', 'column_filter'), fp('rf', 'row_filter'), fp('dfl', 'dataframe_list'),
    fp('lc', 'list', {subType: 'column'}));
  const f = mount(funcForm(call));
  assert.deepEqual([...f.unsupported], ['cf', 'rf', 'dfl', 'lc']);
  f.dispose();
});

// --- table auto-fill (FP-W3-3) ---

w3('auto-fill: the current table wins, one write into the call, and the field opens showing it', () => {
  open(demog());
  const b = open(alt());
  shell.t = b;
  const call = callOf(tableP());
  const {f, sets} = wired(call);
  assert.deepEqual(sets, [['df', b]], 'a real write, exactly one (fpe:62-64, P-W3-6)');
  assert.equal(call.inputs.df, b, 'the OBJECT lands in the call');
  assert.equal(f.getInput('df').value.value, 'alt', 'the input holds the NAME');
  assert.equal(selectOf(f, 'df').value, 'alt', 'auto-fill ran before the field seeded');
  f.dispose();
});

w3('auto-fill: no current table falls back to the first open one; zero tables stay null, empty select', () => {
  const a = open(demog());
  open(alt());
  const first = callOf(tableP());
  const f1 = mount(funcForm(first));
  assert.equal(first.inputs.df, a);
  f1.dispose();

  shell.dart.tables = [];
  const none = callOf(tableP());
  const {f, sets} = wired(none);
  assert.deepEqual(sets, [], 'nothing to fill is no write');
  assert.equal(none.inputs.df ?? null, null);
  const placeholder = selectOf(f, 'df').children;
  assert.equal(placeholder.length, 1, 'no crash, no disable (P-W3-5)');
  assert.equal(placeholder[0].value, '');
  assert.equal(placeholder[0].disabled, true, 'a disabled placeholder, never a pickable blank');
  assert.equal(placeholder[0].textContent, 'No open tables — open or import one');
  f.dispose();
});

w3('skipTableAutoFill suppresses the write; a pre-valued param is never overwritten', () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = b;
  const call = callOf(tableP());
  const {f, sets} = wired(call, {skipTableAutoFill: true});
  assert.deepEqual(sets, [], 'the forms.ts:226 option, honored');
  assert.equal(call.inputs.df ?? null, null);
  f.dispose();

  const valued = callOf(tableP());
  valued.setParamValue('df', a);
  const f2 = mount(funcForm(valued));
  assert.equal(valued.inputs.df, a, 'a held value beats the current-table preference');
  assert.equal(f2.getInput('df').value.value, 'demog');
  f2.dispose();
});

w3('rebind re-runs auto-fill for the new call\'s null table params only', () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const f = mount(funcForm(callOf(tableP())));
  const callB = callOf(tableP());
  f.source = callB;
  assert.equal(callB.inputs.df, a, 'a null param on the new call is filled');

  const callC = callOf(tableP());
  callC.setParamValue('df', b);
  f.source = callC;
  assert.equal(callC.inputs.df, b, 'a valued one is untouched');
  assert.equal(f.getInput('df').value.value, 'alt', 'and the field follows it');
  f.dispose();
});

w3('rebind: a dependent orphans with its missing table param, and revives when it returns', async () => {
  const a = open(demog());
  shell.t = a;
  const withTable = () => callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}));
  const f = mount(funcForm(withTable()));
  assert.equal(f.getInput('weight').enabled, true);

  const noTable = new Func('fceNoTable',
    {inputs: [columnP('weight', {options: {type: 'numerical'}})]}).prepare({});
  f.source = noTable;
  assert.equal(f.getInput('df').enabled, false, 'the table field orphans');
  assert.equal(f.getInput('weight').enabled, false, 'and the dependent orphans WITH it');
  f.getInput('weight').value.value = 'age';
  await flush();
  assert.equal(noTable.inputs.weight ?? null, null,
    'no write resolved against the OLD call\'s table');

  const back = withTable();
  f.source = back;
  assert.equal(f.getInput('weight').enabled, true, 'revives');
  assert.deepEqual(colNames(f, 'weight'), ['age', 'height', 'weight']);
  assert.equal(back.inputs.weight, a.columns.byName('weight'), 'auto-picked against the new call');
  f.dispose();
});

// --- conversion (FP-W3-7: objects into the call, names in the input) ---

w3('a table pick writes the frame OBJECT into the call, never a string', async () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const call = callOf(tableP());
  const {f, sets} = wired(call);
  sets.length = 0;

  pickTable(f, 'df', 'alt');
  await flush();
  assert.equal(sets.length, 1);
  assert.equal(sets[0][1], b, 'dart identity — resolved via the open-tables list');
  assert.notEqual(typeof call.inputs.df, 'string',
    'a string would come back as a Resolve* FuncCall (P-W3-2)');
  f.dispose();
});

w3('readback: a frame shows its name; a Resolve-shaped value ({name, func}) reads as null, no echo', async () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const call = callOf(tableP());
  const {f, sets} = wired(call);
  call.setParamValue('df', b);
  await flush();
  assert.equal(f.getInput('df').value.value, 'alt', 'an external frame write shows its name');
  sets.length = 0;

  call.setParamValue('df', {name: 'demog', func: {}});
  sets.length = 0;
  await flush();
  assert.equal(f.getInput('df').value.value, null,
    'a pending resolver never leaks a bogus name (the \'func\' in v guard)');
  assert.deepEqual(sets.filter(([n]) => n === 'df'), [], 'and the refresh writes nothing back');
  f.dispose();
});

w3('readback: column object → name; array and ColumnList-shaped column_list values render their names', async () => {
  const a = open(demog());
  shell.t = a;
  const call = callOf(tableP(), columnP('num', {options: {type: 'numerical'}}), columnsP('metrics'));
  const f = mount(funcForm(call));
  call.setParamValue('num', a.columns.byName('height'));
  await flush();
  assert.equal(f.getInput('num').value.value, 'height');

  call.setParamValue('metrics', [a.columns.byName('age'), a.columns.byName('height')]);
  await flush();
  assert.deepEqual(f.getInput('metrics').value.value, ['age', 'height']);
  assert.equal(summaryOf(f, 'metrics'), '(2) age, height', 'the Dart summary shape (P-W3-9)');

  call.setParamValue('metrics', [a.columns.byName('age'), {name: 'pending', func: {}}]);
  await flush();
  assert.deepEqual(f.getInput('metrics').value.value, ['age'],
    'a resolver-FuncCall element never leaks its func name (the asNamed guard)');

  const listShaped = new DataFrame([{name: 'age', type: 'int'}], [], 'listShape').columns;
  call.setParamValue('metrics', listShaped);
  await flush();
  assert.deepEqual(f.getInput('metrics').value.value, ['age'], 'a ColumnList reads through names()');
  f.dispose();
});

// --- dependent rewiring (FP-W3-4) ---

w3('a table pick retargets dependents and re-defaults them — ONE write per column, the re-default wins', async () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const call = callOf(tableP(),
    columnP('weight', {options: {type: 'numerical'}}),
    columnP('sex', {options: {type: 'categorical'}}));
  const {f, sets} = wired(call);
  assert.equal(call.inputs.weight, a.columns.byName('weight'), 'auto-picked into the call (P-W3-3)');
  assert.equal(call.inputs.sex, a.columns.byName('sex'));
  assert.deepEqual(colNames(f, 'weight'), ['age', 'height', 'weight']);
  assert.deepEqual(colNames(f, 'sex'), ['sex']);
  sets.length = 0;

  pickTable(f, 'df', 'alt');
  await flush();
  assert.deepEqual(colNames(f, 'weight'), ['weight'], 'the offer follows the new frame through the filter');
  assert.deepEqual(colNames(f, 'sex'), ['label']);
  assert.equal(call.inputs.weight, b.columns.byName('weight'), 're-defaulted into the call (P-W3-4)');
  assert.equal(call.inputs.sex, b.columns.byName('label'));
  const sexWrites = sets.filter(([n]) => n === 'sex');
  assert.equal(sets.filter(([n]) => n === 'weight').length, 1,
    'exactly one write per dependent column per table change — no prune/re-default double write');
  assert.equal(sexWrites.length, 1);
  assert.equal(sexWrites[0][1], b.columns.byName('label'),
    'the re-default, never the prune\'s null, wins when a match exists');
  f.dispose();
});

w3('no-match on the new table re-defaults to NULL in the call (P-W3-10 — Dart parity)', async () => {
  const a = open(demog());
  open(new DataFrame([{name: 'x', type: 'int'}], [], 'nums'));
  shell.t = a;
  const call = callOf(tableP(), columnP('cat', {options: {type: 'categorical'}}));
  const {f, sets} = wired(call);
  assert.equal(call.inputs.cat, a.columns.byName('sex'));
  sets.length = 0;

  pickTable(f, 'df', 'nums');
  await flush();
  assert.deepEqual(colNames(f, 'cat'), [], 'the empty combo — enabled, no special disable');
  assert.equal(call.inputs.cat, null);
  assert.deepEqual(sets.filter(([n]) => n === 'cat'), [['cat', null]], 'the null IS written');
  f.dispose();
});

w3('a table with no passing columns names the cause instead of the generic required error', async () => {
  const a = open(demog());
  open(new DataFrame([{name: 'x', type: 'int'}], [], 'nums'));
  shell.t = a;
  const call = callOf(tableP(), columnP('cat', {nullable: false, options: {type: 'categorical'}}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('cat').validity.peek(), null, 'a passing auto-pick keeps the field valid');

  pickTable(f, 'df', 'nums');
  await flush();
  assert.equal(f.getInput('cat').validity.peek(), 'No categorical columns in \'nums\'');
  assert.equal(f.getInput('cat').root.querySelector('.u2-input-error').textContent,
    'No categorical columns in \'nums\'', 'the cause is on the field, not the generic message');

  pickTable(f, 'df', 'demog');
  await flush();
  assert.equal(f.getInput('cat').validity.peek(), null, 'a table with matches re-defaults and clears it');
  f.dispose();
});

w3('the no-match message derives its kind from the filter: plain when unfiltered, semType when set', () => {
  const plain = open(new DataFrame([{name: 's', type: 'string'}], [], 'plain'));
  shell.t = plain;
  const empty = open(new DataFrame([], [], 'bare'));
  const call = callOf(tableP(),
    columnP('num', {nullable: false, options: {type: 'numerical'}}),
    columnP('mol', {nullable: false, semType: 'Molecule'}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('num').validity.peek(), 'No numerical columns in \'plain\'');
  assert.equal(f.getInput('mol').validity.peek(), 'No Molecule columns in \'plain\'');

  call.setParamValue('df', empty);
  assert.equal(f.getInput('num').validity.peek(), 'No numerical columns in \'bare\'');
  assert.equal(f.getInput('mol').validity.peek(), 'No Molecule columns in \'bare\'');
  f.dispose();
});

w3('dependent required errors are suppressed while the table param is null, and resume with a table', async () => {
  const call = callOf(tableP('df', {nullable: false}),
    columnP('cat', {nullable: false, options: {type: 'categorical'}}),
    columnsP('metrics', {nullable: false}));
  const f = mount(funcForm(call));
  assert.equal(f.getInput('df').validity.peek(), 'Value can\'t be empty',
    'the one root cause stays on the table field');
  assert.equal(f.getInput('cat').validity.peek(), null, 'no dependent echo of it');
  assert.equal(f.getInput('metrics').validity.peek(), null);

  const nums = open(new DataFrame([{name: 'x', type: 'int'}], [], 'nums'));
  call.setParamValue('df', nums);
  await flush();
  assert.equal(f.getInput('cat').validity.peek(), 'No categorical columns in \'nums\'',
    'suppression lifts with the table — the null-table state never masks a real cause');
  assert.equal(f.getInput('metrics').validity.peek(), 'Value can\'t be empty',
    'a satisfiable requirement resumes as the generic message');

  const b = open(demog());
  call.setParamValue('df', b);
  await flush();
  assert.equal(f.getInput('cat').validity.peek(), null, 'auto-picked and valid again');
  f.dispose();
});

w3('an external table write retargets dependents too (#12); a same-frame write causes zero retargets', async () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}));
  const {f, sets} = wired(call);
  sets.length = 0;

  call.setParamValue('df', b);
  await flush();
  assert.deepEqual(colNames(f, 'weight'), ['weight'],
    'the binding rides the param stream, not the input (FE-7 #2/#3)');
  assert.equal(call.inputs.weight, b.columns.byName('weight'));
  assert.equal(f.getInput('df').value.value, 'alt', 'the field follows the external write too');
  sets.length = 0;

  const frameSubs = b.liveSubscriptions();
  call.setParamValue('df', b);
  call.setParamValue('df', b);
  await flush();
  assert.deepEqual(sets.filter(([n]) => n === 'weight'), [],
    'same-frame writes are suppressed at the source (P-W3-2) — zero re-defaults');
  assert.equal(b.liveSubscriptions(), frameSubs, 'and zero frame re-subscriptions');
  f.dispose();
});

w3('an override on the TABLE param keeps its dependents driven through the param stream', async () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const custom = new TextInput({name: 'df'});
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}),
    columnsP('metrics'));
  const f = mount(funcForm(call, {overrides: {df: {input: custom}}}));
  assert.equal(f.getInput('df'), custom, 'the custom editor owns the table field');
  assert.deepEqual([...f.unsupported], [], 'the dependents stay supported');
  assert.equal(call.inputs.df, a, 'auto-fill still writes the frame');
  assert.deepEqual(colNames(f, 'weight'), ['age', 'height', 'weight'], 'the dependent is populated');
  assert.equal(call.inputs.weight, a.columns.byName('weight'), 'and auto-picked into the call');
  assert.equal(f.getInput('weight').enabled, true);

  call.setParamValue('df', b);
  await flush();
  assert.deepEqual(colNames(f, 'weight'), ['weight'],
    'an external write retargets — the binding rides the param, not the input');
  assert.equal(call.inputs.weight, b.columns.byName('weight'));
  assert.equal(summaryOf(f, 'metrics'), '(0)');

  custom.value.value = 'demog';
  await flush();
  assert.deepEqual(colNames(f, 'weight'), ['age', 'height', 'weight'],
    'the custom editor\'s own write drives them too (a name resolves via tableByName)');
  assert.equal(call.inputs.weight, a.columns.byName('weight'));
  f.dispose();
  custom.dispose();
});

// --- the columns field ---

w3('columns: a table change clears the selection to [] written into the call; a pre-set value survives start()', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnsP('metrics', {options: {type: 'numerical'}}));
  call.setParamValue('df', a);
  call.setParamValue('metrics', [a.columns.byName('age'), a.columns.byName('height')]);
  const {f, sets} = wired(call);
  assert.equal(summaryOf(f, 'metrics'), '(2) age, height', 'the pre-set value survives start()');
  assert.deepEqual(sets.filter(([n]) => n === 'metrics'), [], 'no spurious clear at open');

  pickTable(f, 'df', 'alt');
  await flush();
  assert.equal(summaryOf(f, 'metrics'), '(0)', 'changeTable clears the selection (columns.ts:97)');
  const writes = sets.filter(([n]) => n === 'metrics');
  assert.equal(writes.length, 1);
  assert.deepEqual(writes[0][1], [], 'the clear flows into the call');
  assert.deepEqual(call.inputs.metrics, []);
  f.dispose();
});

w3('nullable: false on a column_list flags an empty selection invalid', async () => {
  const a = open(demog());
  shell.t = a;
  const call = callOf(tableP(), columnsP('metrics', {nullable: false}));
  const f = mount(funcForm(call));
  assert.equal(f.isValid, false, 'ObjectForm.isEmpty([]) is false — the empty-array check must trip');
  call.setParamValue('metrics', [a.columns.byName('age')]);
  await flush();
  assert.equal(f.isValid, true);
  f.dispose();
});

w3('a rename remaps the held name in column and columns fields without a stale write', async () => {
  const a = open(demog());
  shell.t = a;
  const call = callOf(tableP(),
    columnP('weight', {options: {type: 'numerical'}}), columnsP('metrics'));
  call.setParamValue('metrics', [a.columns.byName('age')]);
  const {f, sets} = wired(call);
  const picked = call.inputs.weight;
  assert.equal(picked, a.columns.byName('weight'));
  sets.length = 0;

  a.columns.byName('weight').name = 'wt';
  await flush();
  assert.equal(f.getInput('weight').value.value, 'wt', 'the value moves first (pickers.ts:22-30)');
  assert.equal(call.inputs.weight, picked, 'the call still holds the same column');
  assert.ok(!sets.some(([, v]) => typeof v === 'string'), 'no name string ever reaches the call');
  assert.ok(!sets.some(([, v]) => v === null), 'and no stale-prune null either');

  a.columns.byName('age').name = 'years';
  await flush();
  assert.deepEqual(f.getInput('metrics').value.value, ['years'], 'the columns field remaps too');
  f.dispose();
});

w3('dropping the picked column prunes to null IN THE CALL — one write, through the combo\'s own follow', async () => {
  const a = open(demog());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}));
  const {f, sets} = wired(call);
  assert.equal(call.inputs.weight, a.columns.byName('weight'));
  sets.length = 0;

  a.columns.remove('weight');
  await flush();
  assert.equal(f.getInput('weight').value.value, null);
  assert.deepEqual(sets.filter(([n]) => n === 'weight'), [['weight', null]],
    'the prune the frame leaves behind is written, exactly once (FP-P-4 — the follow moved into the input)');
  assert.deepEqual(colNames(f, 'weight'), ['age', 'height'], 'the offer follows the frame');
  f.dispose();
});

// --- column-type filtering (FP-W3-5) ---

w3('columnTypeFilter: the magic strings, a concrete type and no filter all shape the offer', () => {
  const frame = open(new DataFrame([
    {name: 'i', type: 'int'}, {name: 'd', type: 'double'}, {name: 's', type: 'string'},
    {name: 'b', type: 'bool'}, {name: 't', type: 'datetime'},
  ], [], 'mixed'));
  shell.t = frame;
  const call = callOf(tableP(),
    columnP('num', {options: {type: 'numerical'}}),
    columnP('nodt', {options: {type: 'numerical_no_datetime'}}),
    columnP('cat', {options: {type: 'categorical'}}),
    columnP('dt', {options: {type: 'datetime'}}),
    columnP('catdt', {options: {type: 'categorical_or_datetime'}}),
    columnP('exact', {options: {type: 'bool'}}),
    columnP('any'));
  const f = mount(funcForm(call));
  assert.deepEqual(colNames(f, 'num'), ['i', 'd', 't'], 'datetime is numerical, as in ddt');
  assert.deepEqual(colNames(f, 'nodt'), ['i', 'd']);
  assert.deepEqual(colNames(f, 'cat'), ['s', 'b']);
  assert.deepEqual(colNames(f, 'dt'), ['t']);
  assert.deepEqual(colNames(f, 'catdt'), ['s', 'b', 't']);
  assert.deepEqual(colNames(f, 'exact'), ['b'], 'anything else non-null is a concrete type match');
  assert.deepEqual(colNames(f, 'any'), ['i', 'd', 's', 'b', 't'], 'null passes everything');
  f.dispose();
});

w3('semType and maxCategories constrain the offer', () => {
  const frame = open(new DataFrame([
    {name: 'a', type: 'string', semType: 'Sex'}, {name: 'v', type: 'string'},
    {name: 'g', type: 'string'}, {name: 'n', type: 'int'},
  ], [{a: 'M', v: 'x', g: 'p', n: 1}, {a: 'F', v: 'y', g: 'p', n: 2}], 'typed'));
  shell.t = frame;
  const capped = columnP('cap');
  capped.dart.maxCategories = 1;
  const call = callOf(tableP(), columnP('bySem', {semType: 'Sex'}), capped);
  const f = mount(funcForm(call));
  assert.deepEqual(colNames(f, 'bySem'), ['a'], 'semType must match when the prop carries one');
  assert.deepEqual(colNames(f, 'cap'), ['g', 'n'],
    'a categorical over the cap is out; a non-categorical is never capped');
  f.dispose();
});

// --- column auto-pick (fpe:162-179) ---

w3('auto-pick: a semType pick prefers the molecule-ish names, else the first match', () => {
  const frame = open(new DataFrame([
    {name: 'foo', type: 'string', semType: 'Molecule'},
    {name: 'structure', type: 'string', semType: 'Molecule'},
    {name: 'sex', type: 'string', semType: 'Sex'},
  ], [], 'chem'));
  shell.t = frame;
  const call = callOf(tableP(), columnP('mol', {semType: 'Molecule'}), columnP('who', {semType: 'Sex'}));
  const f = mount(funcForm(call));
  assert.equal(call.inputs.mol, frame.columns.byName('structure'),
    'the first named structure|smiles|canonical_smiles|molecule');
  assert.equal(call.inputs.who, frame.columns.byName('sex'), 'else the first match');
  f.dispose();
});

w3('auto-pick: the nearest name wins (1-char and longer alike); no candidates → null', () => {
  const frame = open(new DataFrame([
    {name: 'ax', type: 'int'}, {name: 'xyz', type: 'int'},
    {name: 'agex', type: 'string'}, {name: 'height', type: 'string'},
  ], [], 'names'));
  shell.t = frame;
  const call = callOf(tableP(),
    columnP('x', {options: {type: 'numerical'}}),
    columnP('age', {options: {type: 'categorical'}}),
    columnP('none', {options: {type: 'datetime'}}));
  const f = mount(funcForm(call));
  assert.equal(call.inputs.x, frame.columns.byName('ax'),
    'Levenshtein-nearest for the 1-char param name');
  assert.equal(call.inputs.age, frame.columns.byName('agex'), 'Jaro-Winkler-nearest otherwise');
  assert.equal(call.inputs.none ?? null, null, 'nothing passes the filter — no pick, no write');
  assert.deepEqual(colNames(f, 'none'), []);
  f.dispose();
});

w3('auto-fill and auto-pick ride over skipDefaultInit; column lists are never auto-picked', () => {
  const a = open(demog());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}), columnsP('metrics'));
  const f = mount(funcForm(call, {skipDefaultInit: true}));
  assert.equal(call.inputs.df, a, 'the Dart guard covers only options.default (fpe:558)');
  assert.equal(call.inputs.weight, a.columns.byName('weight'), 'auto-pick is unguarded');
  assert.equal(call.inputs.metrics ?? null, null, 'a column list is never auto-picked (P-W3-3)');
  f.dispose();
});

// --- the closed value table (#11) ---

w3('closing the value table prunes to null IN THE CALL and empties the dependents', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}));
  const {f, sets} = wired(call);
  sets.length = 0;

  shell.closeTable(a);
  await flush();
  assert.equal(call.inputs.df, null, 'never the Dart stale-frame-plus-lying-select (P-W3-12)');
  assert.equal(f.getInput('df').value.value, null);
  assert.deepEqual(items(f, 'df'), ['alt'], 'the remaining table stays on offer, unchosen');
  assert.deepEqual(colNames(f, 'weight'), [], 'dependents retarget to the null table');
  assert.equal(call.inputs.weight, null);
  f.dispose();
});

// --- options (FP-W3-8) ---

w3('showTableSelectors: false hides the table field root; auto-fill and rewiring still run', async () => {
  const a = open(demog());
  const b = open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}));
  const f = mount(funcForm(call, {showTableSelectors: false}));
  assert.equal(f.getInput('df').root.hidden, true, 'the Dart htmlSetDisplay parity (fpe:73)');
  assert.notEqual(f.getInput('weight').root.hidden, true, 'only the table selectors hide');
  assert.equal(call.inputs.df, a, 'auto-fill still ran');
  assert.equal(call.inputs.weight, a.columns.byName('weight'));

  call.setParamValue('df', b);
  await flush();
  assert.deepEqual(colNames(f, 'weight'), ['weight'], 'wiring stays live behind the hidden selector');
  assert.equal(call.inputs.weight, b.columns.byName('weight'));
  f.dispose();
});

// --- the import action ---

w3('the import action writes the imported frame into the call', async () => {
  const call = callOf(tableP());
  const {f, sets} = wired(call);
  assert.deepEqual(sets, [], 'nothing open, nothing filled');
  const df = f.getInput('df');
  assert.ok(df.root.querySelector('.u2-input-options [aria-label="Open file"]'),
    'the folder-open import rides the table field (table_input.dart:37-58)');
  let picker = df.root.querySelector('input[type=file]');
  if (picker === null) {
    // the action API makes the picker transient (FP-P-2): the icon mints it
    fire(df.root.querySelector('.u2-input-options [aria-label="Open file"]'), 'click');
    await flush();
    picker = df.root.querySelector('input[type=file]') ?? document.querySelector('input[type=file]');
  }

  picker.files = [{name: 'plates.csv', text: async () => 'a,b\n1,2',
    arrayBuffer: async () => new Uint8Array([1]).buffer}];
  fire(picker, 'change');
  await until(() => (call.inputs.df ?? null) !== null, 'the imported table to reach the call');
  const imported = shell.tables.find((t) => t.name === 'plates');
  assert.ok(imported, 'the file joined the workspace');
  assert.equal(call.inputs.df, imported, 'and flowed into the call as the OBJECT');
  assert.equal(f.getInput('df').value.value, 'plates');
  f.dispose();
});

w3('the form\'s table field carries the platform icons when the doors exist (FP-P-5)', () => {
  globalThis.grok_UI_PickTableFromFiles = async () => null;
  globalThis.grok_UI_PickTableFromQuery = async () => null;
  try {
    open(demog());
    const f = mount(funcForm(callOf(tableP())));
    assert.deepEqual(
      f.getInput('df').root.querySelectorAll('.u2-input-options [data-u2="icon-button"]')
        .map((el) => el.getAttribute('aria-label')),
      ['Open file', 'Add file from Files', 'Query database'],
      'feature detection reaches the form field for free');
    f.dispose();
  } finally {
    delete globalThis.grok_UI_PickTableFromFiles;
    delete globalThis.grok_UI_PickTableFromQuery;
  }
});

// --- W2 interplay (ParamSource over fielded table/column deps) ---

w3('a choices source with a dataframe dep re-evals exactly once per table pick', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  let runs = 0;
  provider('colsOf', [{name: 'df', propertyType: 'dataframe'}], async (i) => {
    runs++;
    return i.df?.columns.names() ?? [];
  });
  const call = callOf(tableP(), fp('city', 'string', {options: {choices: 'colsOf(df)'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(runs, 1);
  assert.deepEqual(items(f, 'city'), ['age', 'height', 'weight', 'sex'],
    'the provider read the auto-filled frame (P-W3-8)');

  pickTable(f, 'df', 'alt');
  await until(() => items(f, 'city').includes('label'), 'the dependent source to re-read the new frame');
  await sleep(250);
  assert.equal(runs, 2, 'one table pick is one re-eval through the debounce');
  f.dispose();
});

w3('a table pick cascades through the column re-default to its dependent source exactly once', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  let runs = 0;
  provider('byCol', [{name: 'weight', propertyType: 'column'}], async (i) => {
    runs++;
    return [i.weight?.name ?? 'none'];
  });
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}),
    fp('tag', 'string', {options: {choices: 'byCol(weight)'}}));
  const f = mount(funcForm(call));
  await f.settled;
  assert.equal(runs, 1);

  pickTable(f, 'df', 'alt');
  await until(() => runs === 2, 'the re-default write to re-trigger the column-dep source');
  await sleep(250);
  assert.equal(runs, 2, 'the cascade terminates — table pick → one column write → one re-eval');
  f.dispose();
});

w3('auto-fill lands before settled\'s computed defaults, and the table form settles (R5)', async () => {
  const a = open(demog());
  shell.t = a;
  provider('sizeOf', [{name: 'df', propertyType: 'dataframe'}], async (i) => i.df?.rowCount ?? -1);
  const call = callOf(tableP(), fp('count', 'int', {options: {default: 'sizeOf(df)'}}));
  const f = mount(funcForm(call));
  assert.equal(call.inputs.df, a, 'the fill is sync — in the call before any default eval reads it');
  assert.ok(await Promise.race([f.settled.then(() => true), sleep(1500).then(() => false)]),
    'settled resolves — W3 adds nothing async to it');
  assert.equal(call.inputs.count, 2, 'the default eval saw the auto-filled frame');
  f.dispose();
});

// --- the cleared-selection notice and the auto affordance (W3 final round) ---

w3('a table switch that clears a non-empty column selection says so on the field, until the next interaction', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnsP('metrics'));
  call.setParamValue('metrics', [a.columns.byName('age'), a.columns.byName('height')]);
  const f = mount(funcForm(call));
  const root = f.getInput('metrics').root;
  assert.equal(root.querySelector('.u2-param-source-notice'), null, 'no notice at open');

  pickTable(f, 'df', 'alt');
  await flush();
  const notice = root.querySelector('.u2-param-source-notice');
  assert.ok(notice, 'the cleared picks are announced on the field');
  assert.equal(notice.textContent, 'Selection cleared — columns belonged to \'demog\'');
  assert.deepEqual(call.inputs.metrics, [], 'the [] write itself stands (ruled parity)');

  fire(root, 'pointerdown');
  await flush();
  assert.equal(root.querySelector('.u2-param-source-notice'), null,
    'the next interaction with the field dismisses it');
  f.dispose();
});

w3('an already-empty column selection switches tables silently', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnsP('metrics'));
  const f = mount(funcForm(call));
  pickTable(f, 'df', 'alt');
  await flush();
  assert.equal(f.getInput('metrics').root.querySelector('.u2-param-source-notice'), null,
    'nothing was lost, nothing is said');
  f.dispose();
});

w3('auto-picked columns and the auto-filled table wear the auto badge; column lists never do', () => {
  const a = open(demog());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}), columnsP('metrics'));
  const f = mount(funcForm(call));
  assert.ok(f.getInput('df').root.querySelector('.u2-param-auto'), 'the auto-filled table is marked');
  assert.ok(f.getInput('weight').root.querySelector('.u2-param-auto'), 'the auto-pick is marked');
  assert.equal(f.getInput('metrics').root.querySelector('.u2-param-auto'), null,
    'a column list is never auto-picked, so never marked');
  f.dispose();
});

w3('the auto badge survives a retarget re-default and is gone for good after a user pick', async () => {
  const a = open(demog());
  open(alt());
  shell.t = a;
  const call = callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}));
  const f = mount(funcForm(call));
  const weightRoot = f.getInput('weight').root;

  pickTable(f, 'df', 'alt');
  await flush();
  assert.ok(weightRoot.querySelector('.u2-param-auto'), 'a binding re-default keeps the guess marked');
  assert.equal(f.getInput('df').root.querySelector('.u2-param-auto'), null,
    'the table pick was the user\'s own — its badge goes');

  await pickColumnVia(f, 'weight', 'weight');
  assert.equal(weightRoot.querySelector('.u2-param-auto'), null, 'a user pick clears it');

  pickTable(f, 'df', 'demog');
  await flush();
  assert.equal(weightRoot.querySelector('.u2-param-auto'), null,
    'and no later re-default brings it back');
  f.dispose();
});

w3('the auto badge never outlives its value: a no-match re-default and a closed table both clear it', async () => {
  const a = open(demog());
  const nums = open(new DataFrame([{name: 'x', type: 'int'}], [], 'nums'));
  shell.t = a;
  const call = callOf(tableP(), columnP('cat', {options: {type: 'categorical'}}));
  const f = mount(funcForm(call));
  assert.ok(f.getInput('cat').root.querySelector('.u2-param-auto'));

  call.setParamValue('df', nums);
  await flush();
  assert.equal(f.getInput('cat').root.querySelector('.u2-param-auto'), null,
    'nothing picked, nothing marked');
  assert.ok(f.getInput('df').root.querySelector('.u2-param-auto'),
    'an external table write is not a user touch');

  shell.closeTable(nums);
  await flush();
  assert.equal(f.getInput('df').root.querySelector('.u2-param-auto'), null,
    'the pruned table field no longer claims a value');
  f.dispose();
});

// --- scope hygiene ---

w3('dispose and rebind release every frame and param subscription', () => {
  const a = open(demog());
  shell.t = a;
  const addedBase = shell.dart.tableAdded.count;
  const removedBase = shell.dart.tableRemoved.count;
  const make = () => callOf(tableP(), columnP('weight', {options: {type: 'numerical'}}),
    columnsP('metrics'));
  const callA = make();
  const f = mount(funcForm(callA));
  const dfParam = callA.inputParams['df'];
  assert.equal(dfParam.onChanged.count, 2,
    'the field refresh AND the binding — one subscription each, never accumulated');
  const frameSubs = a.liveSubscriptions();
  assert.ok(frameSubs > 0, 'dependents follow the frame');

  const callB = make();
  f.source = callB;
  assert.equal(dfParam.onChanged.count, 0, 'the old call is fully released');
  assert.equal(callB.inputParams['df'].onChanged.count, 2, 'the new one is live');
  assert.equal(a.liveSubscriptions(), frameSubs, 'rebinding never accumulates frame subscriptions');

  f.dispose();
  assert.equal(callB.inputParams['df'].onChanged.count, 0);
  assert.equal(a.liveSubscriptions(), 0, 'frame subscriptions back to baseline');
  assert.equal(shell.dart.tableAdded.count, addedBase, 'the open-tables follow is gone');
  assert.equal(shell.dart.tableRemoved.count, removedBase);
});
