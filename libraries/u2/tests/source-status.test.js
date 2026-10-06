/* The P4 acceptance gap: a source keeps its own clock, and nothing on screen said whether it had
   run. What is pinned here — the status a source reports through the bind protocol, the wording it
   becomes, the balloon an explicit Refresh answers with, and the panel's Status section, which is
   the one LIVE part of a panel whose form fields must keep refreshing in place. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {SpecEditor} from '../src/spec/editor.js';
import {registerAll} from '../src/spec/registrations.js';
import {backends} from '../src/sources/backends.js';
import {DataFrame, Property, Stream} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {SpecNodeRef} = await import('../src/dg/designer/node-ref.js');
const {SpecNodeHandler, disposePanel} = await import('../src/dg/designer/handler.js');
const {sourceStatus, statusText, refreshSource} = await import('../src/dg/designer/source-status.js');
const {shell} = await import('datagrok-api/grok');

const FRAMES = {name: 'Frames', kind: 'function', inputs: [],
  outputs: [new Property('orders', 'dataframe')]};

const frame = (rows) => new DataFrame([{name: 'customer', type: 'string'}], rows);

const SPEC = {
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-func-source', name: 'orders', props: {func: 'Frames', debounce: 0}}],
  root: {tag: 'u2-form', name: 'form', children: [
    {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'}},
  ]},
};

/** One designed document per test, with the runner the source calls faked. */
function status(name, run, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const warn = console.warn;
    const info = shell.info;
    const error = shell.error;
    const said = [];
    console.warn = () => {};
    shell.info = (message) => said.push(['info', message]);
    shell.error = (message) => said.push(['error', message]);
    backends.funcRunner = {find: (n) => n === 'Frames' ? FRAMES : null, run};
    const reg = new Registry();
    registerAll(reg);
    const spec = JSON.parse(JSON.stringify(SPEC));
    const instance = renderSpec(spec, new SpecContext(), reg, {designTime: true});
    const editor = new SpecEditor(instance);
    const source = instance.nodes().get(spec.components[0]);
    try {
      await body({instance, source, said, handler: new SpecNodeHandler(),
        panel: () => new SpecNodeHandler()
          .renderProperties(new SpecNodeRef(instance, spec.components[0], editor))});
    } finally {
      disposePanel();
      instance.dispose();
      delete backends.funcRunner;
      shell.info = info;
      shell.error = error;
      console.warn = warn;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

/** The Status section as `{row: value}` — an h3 followed by what the handler put under it. */
function section(panel, title) {
  const kids = [...panel.children];
  const i = kids.findIndex((el) => el.tagName === 'H3' && el.textContent === title);
  if (i < 0)
    return null;
  const shown = {};
  for (const tr of kids[i + 1].querySelectorAll('tr'))
    shown[tr.children[0]?.textContent.trim()] = tr.children[1]?.textContent.trim();
  return shown;
}

const ok = () => Promise.resolve({orders: frame([{customer: 'Bayer'}, {customer: 'Roche'}])});

status('a source that has not run says so, and reports its rows once it has', ok,
  async ({source}) => {
    assert.deepEqual(sourceStatus(source), {state: 'idle', count: 0, counts: 'rows', error: null},
      'a plain function never auto-runs at design time — that is what has to be legible');
    assert.equal(statusText(sourceStatus(source)), 'not run yet — use Refresh');

    await source.getFunctions().find((f) => f.name === 'refresh').apply();
    await flush();
    assert.deepEqual(sourceStatus(source), {state: 'ready', count: 2, counts: 'rows', error: null});
    assert.equal(statusText(sourceStatus(source)), 'ready · 2 rows');
  });

status('a failure is the status, not silence', () => Promise.reject(new Error('connection refused')),
  async ({source, said}) => {
    await refreshSource(source, 'orders');
    assert.equal(sourceStatus(source).error, 'connection refused');
    assert.equal(statusText(sourceStatus(source)), 'connection refused');
    assert.deepEqual(said, [['error', 'orders: connection refused']]);
  });

status('an explicit refresh runs the source and names what came back', ok, async ({source, said}) => {
  await refreshSource(source, 'orders');
  assert.deepEqual(said, [['info', 'orders: ready · 2 rows']],
    'no balloon, no spinner, no change is how a refresh that did nothing reads');
});

status('the panel section reports state, rows and the failure — and follows the source', ok,
  async ({source, panel}) => {
    const shown = panel();
    assert.deepEqual(section(shown, 'Status'),
      {State: 'idle — not run yet; use Refresh', Rows: '0'});

    await source.getFunctions().find((f) => f.name === 'refresh').apply();
    await flush();
    assert.deepEqual(section(shown, 'Status'), {State: 'ready', Rows: '2'},
      'the same panel, redrawn in place — the fields around it were never rebuilt');
  });

/* A table source and an entity source have no `designData` at all — by design, they are always
   live. Its ABSENCE from the panel is what reads as a missing feature, so the section says it. */
test('a source with no design-time policy says so instead of showing nothing', async () => {
  const live = Scope.liveCount;
  const df = frame([{customer: 'Bayer'}]);
  backends.workspace = {table: () => df, tableNames: () => ['orders'],
    onTablesChanged: new Stream()};
  const reg = new Registry();
  registerAll(reg);
  const spec = {$schema: 'dg-ui/1',
    components: [{tag: 'u2-table-source', name: 'orders', props: {table: 'orders'}}],
    root: {tag: 'u2-form', name: 'form'}};
  const instance = renderSpec(spec, new SpecContext(), reg, {designTime: true});
  const editor = new SpecEditor(instance);
  try {
    const shown = new SpecNodeHandler()
      .renderProperties(new SpecNodeRef(instance, spec.components[0], editor));
    assert.deepEqual(section(shown, 'Status'),
      {Rows: '1', 'Design data': 'always live — this source has no design-time policy'});
  } finally {
    disposePanel();
    instance.dispose();
    delete backends.workspace;
    resetDom();
    await flush();
  }
  assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
});
