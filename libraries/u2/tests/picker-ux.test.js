/* P4 acceptance — what the two pickers do to a query. The binding tree used to filter matches into
   collapsed branches and reset every expansion the moment the box was cleared, and the function list
   matched on text no row drew. Both are model questions, so both are answered here. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import './dom-shim.js';
import {bindGroups, bindRows} from '../src/dg/designer/bind-picker.js';
import {Func} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {funcEntries, filterFuncs} = await import('../src/dg/designer/func-picker.js');

/** A source with a walkable step under it — the three levels the acceptance pass had to open by
 * hand after every keystroke. */
const model = () => [
  {label: 'reagent : string ⇄', path: '$.reagent'},
  {label: 'demog (u2-table-source)', path: '$.demog', children: () => [
    {label: 'df : dataframe', path: '$.demog.df'},
    {label: 'currentRow : object', path: null, children: () => [
      {label: 'subj : int ⇄', path: '$.demog.currentRow.subj'},
      {label: 'site : string ⇄', path: '$.demog.currentRow.site'},
    ]},
  ]},
];

const labels = (rows) => rows.map((row) => row.label);

test('bindRows: a match is revealed, not merely kept', () => {
  const {rows, expand} = bindRows(model(), 'subj');
  assert.deepEqual(labels(rows), ['demog (u2-table-source)']);
  const source = rows[0];
  assert.ok(Array.isArray(source.children), 'a branch kept for a descendant carries its own rows');
  assert.deepEqual(labels(source.children), ['currentRow : object']);
  assert.deepEqual(labels(source.children[0].children), ['subj : int ⇄']);
  assert.deepEqual(expand, [source.id, source.children[0].id],
    'and every branch between the roots and the match is opened');
});

test('bindRows: a match shows its subtree, not an empty expander', async () => {
  const {rows} = bindRows(model(), 'demog');
  assert.deepEqual(labels(rows), ['demog (u2-table-source)']);
  assert.equal(typeof rows[0].children, 'function', 'the match keeps its lazy branch');
  const children = await rows[0].children();
  assert.deepEqual(labels(children), ['df : dataframe', 'currentRow : object'],
    'the whole subtree, filter or no filter — the match is what was asked for');
});

test('bindRows: an id is a place in the tree, so an expansion survives the query', async () => {
  const all = bindRows(model(), '');
  const filtered = bindRows(model(), 'subj');
  assert.equal(all.rows[1].id, filtered.rows[0].id, 'the same node, the same id');
  const open = await all.rows[1].children();
  assert.equal(open[1].id, filtered.rows[0].children[0].id, 'and so is every branch under it');
});

test('bindRows: nothing matches, nothing is invented', () => {
  assert.deepEqual(bindRows(model(), 'u2Record'), {rows: [], expand: []});
});

test('bindGroups: every root says where it came from, and a component wins its name', () => {
  const instance = {spec: {components: [{tag: 'u2-state', name: 'draft'}]},
    ctx: {data: {reagent: 1, draft: 2}}};
  const roots = [{label: 'reagent : string ⇄', path: '$.reagent'},
    {label: 'draft (u2-state)', path: '$.draft'},
    {label: 'nameInput (u2-text-input)', path: '$.nameInput'}];
  assert.deepEqual(bindGroups(instance, roots).map((g) => [g.title, labels(g.nodes)]), [
    ['App data', ['reagent : string ⇄']],
    ['Data sources', ['draft (u2-state)']],
    ['Form controls', ['nameInput (u2-text-input)']],
  ]);
});

test('bindGroups: a heading with nothing under it is not shown', () => {
  const instance = {spec: {}, ctx: {data: {reagent: 1}}};
  assert.deepEqual(bindGroups(instance, [{label: 'reagent : string ⇄', path: '$.reagent'}])
    .map((g) => g.title), ['App data']);
});

test('funcEntries: what no spec can call is not offered', () => {
  const entries = funcEntries([new Func('/molecule.json'), new Func('Chem/detect.js'),
    new Func('demoOrders', {description: 'Orders placed within N days'})]);
  assert.deepEqual(entries.map((entry) => entry.name), ['demoOrders']);
});

test('funcEntries: what the search matches is on the entry, so the row can show it', () => {
  const entries = funcEntries([
    new Func('findAssignable', {friendlyName: 'Find users assignable to issues',
      description: 'Records the assignment'}),
    new Func('demoOrders', {description: 'Orders placed within N days'})]);
  const found = filterFuncs(entries, 'record');
  assert.deepEqual(found.map((entry) => entry.label), ['Find users assignable to issues']);
  assert.equal(found[0].description, 'Records the assignment',
    'the description the match came from is what the row draws');
  assert.deepEqual(filterFuncs(entries, 'u2Record'), [], 'and a miss is a miss');
});
