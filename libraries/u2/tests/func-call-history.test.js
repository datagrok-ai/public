/* FuncCallHistoryBrowser + FuncCallInput + the FuncForm history icon: the dapi paging chain,
   selection coherence across page appends and name changes, popup commits, and applyHistory's
   value copy with dataframe materialization. `DG`/`grok` come from tests/dg-stub.mjs; the dapi
   `functions`/`tables` sources are per-test fakes assigned onto the stub's `dapi`. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {signal} from '../src/core/signals.js';
import {Registry} from '../src/spec/registry.js';
import {Func, FuncCall, Property} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const {FuncCallHistoryBrowser, funcCallHistoryBrowser} =
  await import('../src/dg/entities/func-call-history-browser.js');
const {FuncCallInput} = await import('../src/dg/inputs/func-call-input.js');
const {FuncCallForm} = await import('../src/dg/funcs/func-form.js');
const {applyHistory} = await import('../src/dg/funcs/func-history.js');
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');

const sin = () => new Func('Sin', {inputs: [new Property('x', 'double')]});

/** A stored run double: id + scalar input values. */
function run(id, x) {
  const call = new FuncCall(sin(), {x});
  call.id = id;
  return call;
}

/** A recording `grok.dapi.functions.calls` chain; `pages` answers list() by pageNumber,
 * `requests.saved` collects what save() was handed. */
function fakeCalls(pages, byId = {}) {
  const requests = [];
  requests.saved = [];
  grok.dapi.functions = {
    get calls() {
      const rec = {};
      const src = {
        allPackageVersions: () => src,
        filter(q) { rec.filter = q; return src; },
        include(s) { rec.include = s; return src; },
        order(f, desc) { rec.order = f; rec.desc = desc; return src; },
        list(options) {
          requests.push({...rec, ...options});
          return Promise.resolve(pages[options.pageNumber] ?? []);
        },
        find: (id) => Promise.resolve(byId[id]),
        save(call) {
          requests.saved.push(call);
          return Promise.resolve(call);
        },
      };
      return src;
    },
  };
  return requests;
}

function ui(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      delete grok.dapi.functions;
      delete grok.dapi.tables;
      grok.shell.dart.writes.length = 0;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function rows(root = document.body) {
  return root.querySelectorAll('.u2-list-row[data-u2-funccall]');
}

ui('an empty functionName asks for a function and requests nothing; a name pages the dapi chain', async () => {
  const requests = fakeCalls({0: [run('a', 1), run('b', 2)]});
  const browser = funcCallHistoryBrowser({});
  document.body.append(browser.root);
  await flush();
  assert.equal(requests.length, 0);
  assert.equal(browser.root.querySelector('.u2-fch-empty-message').textContent,
    'Pick a function to see its runs.');

  browser.functionName.value = 'Chem:Sin';
  await flush();
  assert.equal(requests.length, 1);
  assert.equal(requests[0].filter, 'func.name="Sin"', 'short name, not the nqName');
  assert.equal(requests[0].include, 'session.user,func.package,func.params,inputs,options');
  assert.deepEqual([requests[0].order, requests[0].desc], ['started', true]);
  assert.equal(requests[0].pageNumber, 0);
  const list = browser.root.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  assert.deepEqual([...rows()].map((r) => r.getAttribute('data-u2-funccall')), ['a', 'b']);
  assert.equal(browser.call('b').id, 'b');
  browser.dispose();
});

ui('quotes and backslashes are sanitized out of the filter value', async () => {
  const requests = fakeCalls({});
  const browser = funcCallHistoryBrowser({functionName: 'Bad"\\Name'});
  document.body.append(browser.root);
  await flush();
  assert.equal(requests[0].filter, 'func.name="BadName"');
  browser.dispose();
});

ui('selection lands in the signal and the shell; a name change resets and clears it once', async () => {
  const changes = [];
  fakeCalls({0: [run('a', 1), run('b', 2)]});
  const browser = funcCallHistoryBrowser({functionName: 'Sin', onChanged: (c) => changes.push(c?.id ?? null)});
  document.body.append(browser.root);
  await flush();
  const list = browser.root.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  fire(rows()[1], 'click');
  assert.equal(browser.selected.value.id, 'b');
  assert.deepEqual(changes, ['b']);
  assert.equal(grok.shell.o.id, 'b', 'selection becomes the current object');

  fakeCalls({0: [run('c', 3)]});
  browser.functionName.value = 'Cos';
  await flush();
  assert.equal(browser.selected.value, null);
  assert.deepEqual(changes, ['b', null], 'cleared exactly once, no phantom emission');
  browser.dispose();
});

ui('setCurrentObject: false keeps the shell out of it; activate fires on Enter and dblclick', async () => {
  const activated = [];
  fakeCalls({0: [run('a', 1)]});
  const browser = funcCallHistoryBrowser({functionName: 'Sin', setCurrentObject: false,
    onActivate: (c) => activated.push(c.id)});
  document.body.append(browser.root);
  await flush();
  const list = browser.root.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  fire(rows()[0], 'click');
  assert.equal(grok.shell.dart.writes.length, 0);
  fire(rows()[0], 'dblclick');
  fire(list, 'keydown', {key: 'Enter'});
  assert.deepEqual(activated, ['a', 'a']);
  browser.dispose();
});

ui('bindProps carry a writable functionName and a read-only selected id step', async () => {
  fakeCalls({0: [run('a', 1)]});
  const browser = new FuncCallHistoryBrowser({functionName: 'Sin'});
  document.body.append(browser.root);
  await flush();
  const props = Object.fromEntries(browser.bindProps().map((p) => [p.name, p]));
  assert.equal(props.functionName.writable, true);
  assert.equal(props.selected.writable, false);
  const list = browser.root.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  fire(rows()[0], 'click');
  assert.equal(browser.bindStep('selected').value, 'a');
  fakeCalls({});
  browser.bindStep('functionName').value = 'Tan';
  assert.equal(browser.functionName.value, 'Tan');
  browser.dispose();
});

ui('FuncCallInput: a row click commits the id and closes; dismissals never write', async () => {
  fakeCalls({0: [run('a', 1), run('b', 2)]});
  const changes = [];
  const input = new FuncCallInput({label: 'Run', functionName: 'Sin', onChanged: (v) => changes.push(v)});
  document.body.append(input.root);
  const control = input.root.querySelector('.u2-func-call-input');
  assert.equal(input.root.dataset.u2, 'func-call-input');
  assert.equal(input.root.querySelector('.u2-func-call-input-name').textContent, 'Pick a run…');
  fire(control, 'click');
  await flush();
  const popup = document.body.querySelector('.u2-func-call-input-popup');
  assert.notEqual(popup.querySelector('[data-u2="func-call-history-browser"]'), null);
  const list = popup.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  fire(rows(popup)[1], 'click');
  assert.equal(input.value.value, 'b');
  assert.deepEqual(changes, ['b']);
  assert.equal(document.body.querySelector('.u2-func-call-input-popup'), null);
  assert.equal(input.root.querySelector('.u2-func-call-input-name').textContent,
    'Sin', 'the picked run labels the editor');

  fire(control, 'click');
  await flush();
  fire(document, 'keydown', {key: 'Escape'});
  assert.equal(document.body.querySelector('.u2-func-call-input-popup'), null);
  assert.equal(input.value.value, 'b', 'Escape closed without a write');
  input.enabled = false;
  fire(control, 'click');
  assert.equal(document.body.querySelector('.u2-func-call-input-popup'), null, 'disabled never opens');
  input.dispose();
});

ui('FuncCallInput resolves a preset id in the background and clears on Backspace', async () => {
  const stored = run('a', 1);
  fakeCalls({}, {a: stored});
  const input = new FuncCallInput({label: 'Run', value: 'a'});
  document.body.append(input.root);
  assert.equal(input.root.querySelector('.u2-func-call-input-name').textContent, 'a', 'raw id until resolved');
  await flush();
  await flush();
  assert.equal(input.root.querySelector('.u2-func-call-input-name').textContent, 'Sin');
  fire(input.root.querySelector('.u2-func-call-input'), 'keydown', {key: 'Backspace'});
  assert.equal(input.value.value, '');
  input.dispose();
});

ui('applyHistory copies materialized input values into the CURRENT call only where names match', async () => {
  const table = {name: 'materialized'};
  grok.dapi.tables = {getTable: async (id) => {
    if (id === 'gone')
      throw new Error('missing');
    return table;
  }};
  const func = new Func('F', {inputs: [
    new Property('x', 'double'), new Property('df', 'dataframe'), new Property('bad', 'dataframe')]});
  const form = new FuncCallForm(new FuncCall(func), {showTableSelectors: true});
  document.body.append(form.root);
  const picked = new FuncCall(func);
  picked.setParamValue('x', 42);
  picked.dart.params.find((p) => p.name === 'df').dart.value = 'table-id';
  picked.dart.params.find((p) => p.name === 'bad').dart.value = 'gone';
  const extra = new FuncCall(new Func('G', {inputs: [new Property('y', 'double')]}), {y: 7});
  await applyHistory(form, picked);
  assert.equal(form.source.inputs['x'], 42);
  assert.equal(form.source.inputs['df'], table, 'the TableInfo id string was materialized');
  assert.equal(form.source.inputs['bad'], undefined, 'a vanished table skips the one param');
  await applyHistory(form, extra);
  assert.equal(form.source.inputs['x'], 42, 'another function\'s params are a no-op');
  form.dispose();
});

ui('the FuncForm history icon lives in the buttons row, toggles with showHistory, and applies a pick', async () => {
  fakeCalls({0: [run('a', 5)]});
  const show = signal(true);
  const form = new FuncCallForm(new FuncCall(sin()), {showHistory: show});
  document.body.append(form.root);
  const icon = form.root.querySelector('[data-u2="ff-history-icon"]');
  assert.notEqual(icon, null);
  assert.equal(icon.closest('.u2-form-buttons') !== null, true, 'the icon rides the buttons row');
  show.value = false;
  assert.equal(icon.hidden, true);
  show.value = true;

  fire(icon, 'click');
  await flush();
  const popup = document.body.querySelector('.u2-func-form-history-popup');
  assert.notEqual(popup, null);
  const list = popup.querySelector('.u2-list');
  list.clientHeight = 220;
  fire(list, 'scroll');
  fire(rows(popup)[0], 'click');
  await flush();
  assert.equal(document.body.querySelector('.u2-func-form-history-popup'), null, 'the pick closed the popup');
  assert.equal(form.source.inputs['x'], 5, 'the picked run\'s value landed in the call');
  assert.equal(form.getInput('x').value.value, 5, 'and the input followed');

  const plain = new FuncCallForm(new FuncCall(sin()), {});
  assert.equal(plain.root.querySelector('[data-u2="ff-history-icon"]'), null, 'no option, no icon');
  plain.dispose();
  form.dispose();
});

ui('the Run button is gated on required fields, runs the call and saves it through dapi', async () => {
  const requests = fakeCalls({});
  const ran = [];
  const show = signal(true);
  const form = new FuncCallForm(new FuncCall(sin()), {showRun: show, onRun: (c) => ran.push(c)});
  document.body.append(form.root);
  const run = form.root.querySelector('[data-u2="ff-run"]');
  assert.notEqual(run, null);
  assert.equal(run.closest('.u2-form-buttons') !== null, true, 'Run rides the buttons row');
  assert.equal(run.disabled, true);
  assert.equal(run.title.startsWith('Missing: '), true, `title="${run.title}"`);

  form.getInput('x').value.value = 0.5;
  await flush();
  assert.equal(run.disabled, false);
  assert.equal(run.title, '');
  fire(run, 'click');
  await flush();
  await flush();
  assert.equal(requests.saved.length, 1, 'the run was saved through the dapi chain');
  assert.equal(requests.saved[0], form.source);
  assert.equal(ran.length, 1, 'onRun fired after run + save');

  show.value = false;
  assert.equal(run.hidden, true);
  const plain = new FuncCallForm(new FuncCall(sin()), {});
  assert.equal(plain.root.querySelector('[data-u2="ff-run"]'), null, 'no option, no button');
  plain.dispose();
  form.dispose();
});

ui('registerPlatformComponents carries the three new metas; a bound callId populates u2-func-form', async () => {
  fakeCalls({}, {a: run('a', 9)});
  Func.registry = [sin()];
  try {
    const reg = new Registry();
    registerPlatformComponents(reg);
    const empty = reg.get('u2-func-form').create({});
    assert.equal(empty.inputs.length, 0, 'no function — the empty call renders nothing');
    empty.dispose();

    const cid = signal(null);
    const form = reg.get('u2-func-form').create({functionName: 'Sin', callId: cid});
    document.body.append(form.root);
    assert.notEqual(form.getInput('x'), undefined, 'the name prepared a call');
    cid.value = 'a';
    await flush();
    await flush();
    assert.equal(form.getInput('x').value.value, 9, 'the bound callId copied the run\'s values in');

    const input = reg.get('u2-func-call-input').create({label: 'Run', functionName: 'Sin'});
    document.body.append(input.root);
    assert.equal(input.root.dataset.u2, 'func-call-input');
    const browser = reg.get('u2-func-call-history-browser').create({functionName: 'Sin'});
    document.body.append(browser.root);
    assert.equal(browser.root.dataset.u2, 'func-call-history-browser');
    browser.dispose();
    input.dispose();
    form.dispose();
  } finally {
    Func.registry = [];
  }
});
