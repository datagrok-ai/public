/* WO-16 — the function picker's model and what a pick becomes: the search over what the client
   registry knows, the name a spec writes for a function, the arguments a fired event carries, and
   the node the tray's `+` wizard inserts. All pure — the dialog itself is e2e's to drive; `DG` and
   `grok` come from tests/dg-stub.mjs, which the tray and the picker reach for at import. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {Registry} from '../src/spec/registry.js';
import {registerAll} from '../src/spec/registrations.js';
import {Func, Package, Property} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {funcEntries, filterFuncs, paramProps, paramValues, eventEntry, eventPick} =
  await import('../src/dg/designer/func-picker.js');
const {funcSourceNode, sourceNode, specName} = await import('../src/dg/designer/tray.js');

const FUNCS = [
  new Func('demoOrders', {friendlyName: 'Demo orders', description: 'Orders placed within N days'}),
  new Func('Sin', {description: 'Sine of a number'}),
  new Func('save', {package: new Package('Chem')}),
  new Func('save', {package: new Package('Bio')}),
];

const reg = new Registry();
registerAll(reg);

const named = (entries) => entries.map((entry) => entry.name);

test('funcEntries: the bare name, qualified only where it would name more than one function', () => {
  const entries = funcEntries(FUNCS);
  assert.deepEqual(named(entries).sort(), ['Bio:save', 'Chem:save', 'Sin', 'demoOrders']);
  const demo = entries.find((entry) => entry.name === 'demoOrders');
  assert.equal(demo.label, 'Demo orders', 'the caption is what the row shows');
  assert.equal(demo.func, FUNCS[0], 'and the platform object is what the chip renders');
  assert.deepEqual(entries.map((entry) => entry.label), ['Demo orders', 'save', 'save', 'Sin'],
    'rows come sorted by what they show');
});

test('filterFuncs: one enumeration, a substring over name, caption and description', () => {
  const entries = funcEntries(FUNCS);
  assert.equal(filterFuncs(entries, '').length, 4, 'an empty query filters nothing');
  assert.deepEqual(named(filterFuncs(entries, 'demo')), ['demoOrders'], 'by name');
  assert.deepEqual(named(filterFuncs(entries, 'demo orders')), ['demoOrders'], 'by caption');
  assert.deepEqual(named(filterFuncs(entries, 'sine')), ['Sin'], 'by description');
  assert.deepEqual(named(filterFuncs(entries, 'chem')), ['Chem:save'], 'and by the package qualifier');
  assert.deepEqual(filterFuncs(entries, '  NOTHING '), []);
});

test('paramProps: a typed editor per input, and the bound ones showing the path they follow', () => {
  const inputs = [
    new Property('days', 'int', {description: 'How many days'}),
    new Property('city', 'string', {choices: ['Kyiv', 'Lviv']}),
    new Property('tags', 'string_list'),
    new Property('table', 'dataframe'),
  ];
  const values = {days: 30};
  const props = paramProps(inputs, values, {city: '$.pick.value'});
  assert.deepEqual(props.map((p) => [p.name, p.type]),
    [['days', 'int'], ['city', 'string'], ['tags', 'list'], ['table', 'dataframe']],
    'string_list edits through the list editor; a type with no editor keeps its own name');

  const [days, city] = props;
  assert.equal(days.get(values), 30);
  days.set(values, 45);
  assert.deepEqual(values, {days: 45}, 'the editor writes the value record, nothing else');
  assert.equal(city.set, undefined, 'a bound param is the binding\'s to change');
  assert.equal(city.get(values), '$.pick.value');
  assert.equal(props[1].description, 'Bound to $.pick.value');
  assert.equal(props[3].choices, null, 'unset choices pass through as the platform answers them');
});

test('paramValues: an empty editor is not a value, but zero and false are', () => {
  assert.deepEqual(paramValues({days: 0, on: false, text: '', missing: undefined, none: null}),
    {days: 0, on: false});
});

test('eventEntry: the string form while nothing is passed, the structured one otherwise', () => {
  assert.equal(eventEntry({name: 'save', args: {}, binds: {}}), 'cmd:save');
  assert.deepEqual(eventEntry({name: 'U2demo:log', args: {text: 'hi', skipped: ''},
    binds: {row: '$.orders.currentRowIdx'}}),
  {cmd: 'cmd:U2demo:log', args: {text: 'hi', row: '$.orders.currentRowIdx'}});
  assert.deepEqual(eventEntry({name: 'echo', args: {text: '$.notAPath'}, binds: {}}),
    {cmd: 'cmd:echo', args: {text: '$$.notAPath'}}, 'a literal that looks like a path escapes once');
});

test('eventPick: what a wired event tells the picker, and the tiers it cannot represent', () => {
  assert.deepEqual(eventPick('cmd:save'), {name: 'save', args: {}, binds: {}});
  assert.deepEqual(eventPick({cmd: 'cmd:echo', args: {text: '$$.notAPath', row: '$.orders.rowCount'}}),
    {name: 'echo', args: {text: '$.notAPath'}, binds: {row: '$.orders.rowCount'}},
    'an escaped literal is unescaped once, a path is a bound argument');
  assert.equal(eventPick(undefined), undefined);
  assert.equal(eventPick('cmd:#save'), undefined, 'code-behind');
  assert.equal(eventPick('cmd:orders.refresh'), undefined, 'a component\'s own function');
  assert.equal(eventPick('save'), undefined, 'and anything that is not a command at all');

  const pick = {name: 'U2demo:log', args: {text: 'hi'}, binds: {row: '$.orders.currentRowIdx'}};
  assert.deepEqual(eventPick(eventEntry(pick)), pick, 'the two are each other\'s inverse');
});

test('the wizard seeds a func-source: params as literals, param binds as dotted keys, no designData',
  () => {
    const node = funcSourceNode(reg.get('u2-func-source'), 'funcSource1',
      {name: 'demoOrders', args: {days: 7, city: ''}, binds: {customer: '$.pickInput.value'}});
    assert.deepEqual(node, {tag: 'u2-func-source', name: 'funcSource1',
      props: {func: 'demoOrders', params: {days: 7}},
      bind: {'params.customer': '$.pickInput.value'}});
    assert.equal('designData' in node.props, false,
      'the source computes the policy from the function\'s kind — seeding it would freeze it');

    const bare = funcSourceNode(reg.get('u2-func-source'), 'funcSource2',
      {name: 'demoOrders', args: {}, binds: {}});
    assert.deepEqual(bare, {tag: 'u2-func-source', name: 'funcSource2', props: {func: 'demoOrders'}},
      'nothing to pass carries no empty containers');
  });

test('the wizard seeds the other sources off the tag\'s defaults, its own choices on top', () => {
  assert.deepEqual(sourceNode(reg.get('u2-entity-source'), 'u2-entity-source', 'entitySource1',
    {entity: 'groups', pageSize: 50, filter: 'admin'}),
  {tag: 'u2-entity-source', name: 'entitySource1',
    props: {entity: 'groups', pageSize: 50, filter: 'admin'}});
  assert.deepEqual(sourceNode(reg.get('u2-entity-source'), 'u2-entity-source', 'entitySource2', {}),
    {tag: 'u2-entity-source', name: 'entitySource2', props: {entity: 'users', pageSize: 20}},
    'the defaults the registry declares');
  assert.deepEqual(sourceNode(reg.get('u2-state'), 'u2-state', 'state1', {}),
    {tag: 'u2-state', name: 'state1', props: {type: 'string'}});
  assert.deepEqual(sourceNode(reg.get('u2-table-source'), 'u2-table-source', 'demog1', {table: 'demog'}),
    {tag: 'u2-table-source', name: 'demog1', props: {table: 'demog'}});
});

test('a table name becomes a name a bind path can address', () => {
  assert.equal(specName('demog', 'tableSource'), 'demog');
  assert.equal(specName('Cars 2024.csv', 'tableSource'), 'Cars2024Csv');
  assert.equal(specName('2024', 'tableSource'), 'tableSource', 'nothing a path could address');
  assert.equal(specName('%%%', 'tableSource'), 'tableSource');
});
