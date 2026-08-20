/* P4.5 WO-1 — the platform drag channel, read without a drag. The payloads are `DG.FileInfo` /
   `DG.Func` doubles of tests/platform-doubles.mjs, whose fields are prototype getters exactly as the
   real entities' are: an implementation that spread one of them would read `undefined` here, which
   is the defect the P4 binding picker shipped with. The wiring half is driven through the stub's
   `ui.drops` registry — payload → patch, with no browser. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import './dom-shim.js';
import {registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';

register('./dg-stub.mjs', import.meta.url);
const DG = await import('datagrok-api/dg');
const ui = await import('datagrok-api/ui');
const {OPEN_FILE, dropNode, funcRef, makeDesignerDroppable, readDrop, tabularExtensions} =
  await import('../src/dg/designer/drop.js');

registerAll(registry);

const file = (path) => new DG.FileInfo(path);
const folder = (path) => new DG.FileInfo(path, {directory: true});
/** Names in the order they are handed out — what `SpecEditor.uniqueNames()` does for real. */
const counter = () => {
  const taken = new Map();
  return (base) => {
    const n = (taken.get(base) ?? 0) + 1;
    taken.set(base, n);
    return `${base}${n}`;
  };
};

function withFuncs(funcs, body) {
  const before = DG.Func.registry;
  DG.Func.registry = funcs;
  try {
    return body();
  } finally {
    DG.Func.registry = before;
  }
}

test('a dropped file becomes an OpenFile source over its own path', () => {
  const {items, refused} = readDrop(file('System:Demo/demog.csv'));
  assert.deepEqual(refused, []);
  assert.deepEqual(items, [{kind: 'file', ref: OPEN_FILE,
    params: {fullPath: 'System:Demo/demog.csv'}, label: 'demog', runs: true}]);
});

test('a file nothing opens as a table, and a folder, are refused BY NAME', () => {
  const reading = readDrop([file('System:Demo/report.pdf'), folder('System:Demo/tables')]);
  assert.deepEqual(reading.items, []);
  assert.equal(reading.refused.length, 2);
  assert.match(reading.refused[0], /report\.pdf.*\.pdf/);
  assert.match(reading.refused[1], /tables.*folder/);
});

test('a multi-selection drags an array — every member is read, in order', () => {
  const items = readDrop([file('a/one.csv'), file('a/two.tsv')]).items;
  assert.deepEqual(items.map((i) => i.label), ['one', 'two']);
  assert.deepEqual(items.map((i) => i.params.fullPath), ['a/one.csv', 'a/two.tsv']);
});

test('a mixed array reads both halves, and anything else is ignored', () => withFuncs(
  [new DG.DataQuery('orders')], () => {
    const items = readDrop([file('a/one.csv'), DG.Func.registry[0], {some: 'object'}, null]).items;
    assert.deepEqual(items.map((i) => [i.kind, i.ref]), [['file', OPEN_FILE], ['func', 'orders']]);
  }));

test('a dropped query is a func item, and the drop never runs it', () => {
  const table = new DG.TableQuery('visual', {friendlyName: 'NW orders'});
  const sql = new DG.DataQuery('handWritten');
  withFuncs([table, sql], () => {
    const items = readDrop([table, sql]).items;
    assert.deepEqual(items.map((i) => [i.label, i.ref, i.runs, Object.keys(i.params).length]),
      [['NW orders', 'visual', false, 0], ['handWritten', 'handWritten', false, 0]]);
  });
});

test('funcRef writes the bare name, and qualifies only what the bare name would not name', () => {
  const query = new DG.DataQuery('orders', {namespace: 'admin'});
  withFuncs([query], () => assert.equal(funcRef(query), 'orders'));
  withFuncs([query, new DG.Func('orders')], () =>
    assert.equal(funcRef(query), 'admin:orders', 'a query has a namespace, not a package'));
});

test('tabularExtensions counts the text formats and every file handler\'s own list', () => {
  const handler = new DG.Func('importSdf', {options: {role: 'file-handler', ext: 'sdf, MOL'}});
  const plain = tabularExtensions();
  assert.equal(plain.has('csv') && plain.has('tsv') && plain.has('txt') && plain.has('d42'), true);
  assert.equal(plain.has('sdf'), false);
  withFuncs([handler], () => {
    const exts = tabularExtensions();
    assert.equal(exts.has('sdf') && exts.has('mol'), true, 'trimmed and lower-cased');
    assert.deepEqual(readDrop(file('a/lib.sdf')).refused, [], 'so an sdf is accepted');
  });
});

test('dropNode: a named func-source, no params key where there are none, no designData ever', () => {
  const unique = counter();
  const [csv, query] = [
    {kind: 'file', ref: OPEN_FILE, params: {fullPath: 'System:Demo/demog.csv'},
      label: 'demog', runs: true},
    {kind: 'func', ref: 'admin:orders', params: {}, label: 'NW orders', runs: false},
  ].map((item) => dropNode(item, registry, unique));

  assert.deepEqual(csv, {tag: 'u2-func-source', name: 'demog1',
    props: {func: OPEN_FILE, params: {fullPath: 'System:Demo/demog.csv'}}});
  // `specName` is the tray's own table-name sanitizer: it camel-joins and leaves the first word
  // as written, so a picked table and a dropped query are named the same way
  assert.deepEqual(query, {tag: 'u2-func-source', name: 'NWOrders1',
    props: {func: 'admin:orders'}});
  assert.equal('params' in query.props, false);
  assert.equal(dropNode({kind: 'func', ref: 'x', params: {}, label: '42', runs: false},
    registry, unique).name, 'funcSource1', 'a label that cannot start an identifier falls back');
});

test('the droppable is wired without the body overlay, and every callback asks active() first', () => {
  let active = true;
  const dropped = [];
  makeDesignerDroppable({element: document.createElement('div'), active: () => active,
    onDragActive: () => {}, onDrop: (reading) => dropped.push(reading)});
  const wired = ui.drops[ui.drops.length - 1];
  assert.equal(wired.dropIndication, false,
    'the default overlay is body-level and owns the only mouseup listener');

  const args = {dragObject: file('a/one.csv'), handled: false};
  assert.equal(wired.acceptDrag(args), true);
  wired.doDrop(args);
  assert.equal(args.handled, true, 'the platform veto flag — but never stopPropagation');
  assert.equal(dropped.length, 1);

  active = false;
  const missed = {dragObject: file('a/two.csv'), handled: false};
  assert.equal(wired.acceptDrag(missed), false);
  wired.doDrop(missed);
  assert.deepEqual([dropped.length, missed.handled], [1, false],
    'a dart registration outlives the designer: a dead one takes nothing');
});

test('a drag carrying nothing we know is not ours, but one we must refuse is', () => {
  makeDesignerDroppable({element: document.createElement('div'), active: () => true,
    onDragActive: () => {}, onDrop: () => {}});
  const wired = ui.drops[ui.drops.length - 1];
  assert.equal(wired.acceptDrag({dragObject: {something: 'else'}}), false);
  assert.equal(wired.acceptDrag({dragObject: file('a/report.pdf')}), true,
    'refusing at acceptDrag would make the drop fall through in silence');
});
