/* The contract of tests/platform-doubles.mjs: every double keeps its fields on prototype getters
   over the dart handle, as the real entity does — so a spread copies nothing but `dart` — and the
   loader hooks serve these very classes, never a second copy. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {BitSet, Column, DataFrame, DataQuery, Entity, FileInfo, Func, Package, Property, Script,
  Shell, TableQuery, UnreadableFileInfo, User} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const DG = await import('datagrok-api/dg');
const grok = await import('datagrok-api/grok');

const demog = () => new DataFrame([{name: 'age', type: 'int', semType: 'Age'}, {name: 'sex'}],
  [{age: 40, sex: 'F'}, {age: 55, sex: 'M'}], 'demog');

const doubles = () => {
  const frame = demog();
  return [
    new Property('days', 'int', {description: 'How many'}),
    new Func('orders'), new DataQuery('sql'), new Script('py'), new TableQuery('visual'),
    new FileInfo('Demo:Files/geo/cities.csv'), new UnreadableFileInfo('keys/server.pem'),
    frame, frame.columns, frame.columns.byName('age'), frame.selection, new BitSet(3),
    new Entity('alice'), new User('alice'), new Package('Chem'), new Shell(),
  ];
};

// stricter than the real DataFrame, which also keeps columns/rows/filter/temp/tags as own fields — none read by u2
test('spreading any double copies the dart handle and nothing else', () => {
  for (const double of doubles()) {
    assert.deepEqual(Object.keys({...double}), ['dart'], double.constructor.name);
    assert.deepEqual(Object.keys(double), ['dart'], double.constructor.name);
  }
});

test('a property answers what it was created with, and what fromOptions refined it with', () => {
  const get = () => 1;
  const p = Property.create('days', 'int', get, null, 7).fromOptions({description: 'How many',
    choices: ['1', '7'], caption: 'Days', nullable: false, category: 'Misc', semType: 'Span'});
  assert.deepEqual([p.name, p.propertyType, p.type, p.get, p.set, p.defaultValue],
    ['days', 'int', 'int', get, null, 7]);
  assert.deepEqual([p.description, p.choices, p.caption, p.nullable, p.category, p.semType],
    ['How many', ['1', '7'], 'Days', false, 'Misc', 'Span']);
  assert.equal(new Property('x', 'string').choices, null, 'unset choices are null, as the Dart field is');
  assert.equal(new Property('x', 'string').caption, undefined, 'what was not given is not there');
});

test('a function answers its declaration, its kind and the call it prepares', async () => {
  const chem = new Package('Chem');
  const inputs = [new Property('smiles', 'string')];
  const f = new Func('save', {package: chem, friendlyName: 'Save it', description: 'Saves',
    inputs, options: {role: 'file-handler'}, run: async (params) => ({ok: params.smiles})});
  assert.deepEqual([f.name, f.friendlyName, f.description, f.package, f.inputs, f.outputs, f.options],
    ['save', 'Save it', 'Saves', chem, inputs, [], {role: 'file-handler'}]);
  assert.deepEqual(await f.prepare({smiles: 'CCO'}).call(), {outputs: {ok: 'CCO'}});
  assert.equal(new Func('orders', {namespace: 'admin'}).nqName, 'admin:orders');
  assert.equal(new Func('orders').nqName, 'orders');

  const query = new TableQuery('visual');
  assert.ok(query instanceof DataQuery && query instanceof Func && query instanceof Entity);
  query.limit = 100;
  assert.equal(query.limit, 100);
  assert.deepEqual(Object.keys(query), ['dart'], 'the limit went behind the handle');

  Func.registry = [f, new Func('save'), new Func('other', {options: {role: 'file-handler'}})];
  assert.equal(Func.find({name: 'save'}).length, 2);
  assert.deepEqual(Func.find({name: 'save', package: 'Chem'}), [f]);
  assert.equal(Func.find({meta: {role: 'file-handler'}}).length, 2);
  Func.registry = [];
});

test('a file reads its names off the path, and its share off the connection', () => {
  const home = new Entity('Files', {nqName: 'Demo:Files'});
  const file = new FileInfo('Demo:Files/geo/cities.csv', {connection: home});
  assert.deepEqual([file.fullPath, file.path, file.fileName, file.extension, file.name],
    ['Demo:Files/geo/cities.csv', 'geo/cities.csv', 'cities.csv', 'csv', 'cities.csv']);
  assert.deepEqual([file.connection, file.isFile, file.isDirectory, file.data], [home, true, false, null]);

  const folder = new FileInfo('a/tables', {directory: true});
  assert.deepEqual([folder.path, folder.extension, folder.isFile, folder.isDirectory], ['a/tables', '', false, true]);
  assert.equal(new FileInfo('Demo:Files').path, '', 'the share itself has no path under it');
  const local = FileInfo.fromBytes('server.pem', new Uint8Array([1]));
  assert.deepEqual([local.connection, local.data.length], [null, 1]);
  assert.throws(() => new UnreadableFileInfo('keys/server.pem').connection, RangeError);
});

test('a frame answers its cells, counts its reads and writes, and fires the platform events', () => {
  const df = demog();
  assert.deepEqual([df.name, df.rowCount, df.currentRowIdx, df.columns.names()], ['demog', 2, 0, ['age', 'sex']]);
  assert.ok(df.selection instanceof BitSet && df.filter instanceof BitSet);
  assert.equal(df.selection.length, 2);
  const age = df.columns.byName('age');
  assert.ok(age instanceof Column);
  assert.deepEqual([age.name, age.type, age.semType, age.isNumerical], ['age', 'int', 'Age', true]);
  assert.deepEqual([df.columns.byName('sex').type, df.columns.byName('sex').semType], ['string', null]);
  assert.equal(df.columns.byName('gone'), null);
  assert.equal(df.columns.contains('sex'), true);
  assert.equal(new Column('when', 'datetime').isNumerical, true, 'as in ddt');

  const seen = [];
  const watched = ['onCurrentRowChanged', 'onValuesChanged', 'onColumnsChanged', 'onColumnNameChanged'];
  const subs = watched.map((event) => df[event].subscribe((x) => seen.push([event, x])));
  assert.equal(df.liveSubscriptions(), 4);

  assert.equal(df.get('age', 1), 55);
  assert.equal(df.dart.reads, 1);
  df.set('age', 1, 56);
  assert.deepEqual(df.dart.writes, [['age', 1, 56]]);
  assert.equal(df.dart.rows[1].age, 56);
  df.currentRowIdx = 0;
  df.currentRowIdx = 1;
  assert.equal(df.dart.idxWrites, 2, 'every write counts');
  age.name = 'years';
  df.columns.remove('sex');
  df.columns.add(new Column('flag', 'bool'));
  assert.deepEqual(seen, [['onValuesChanged', {column: 'age', row: 1}], ['onCurrentRowChanged', 1],
    ['onColumnNameChanged', {args: {oldName: 'age', newName: 'years'}}],
    ['onColumnsChanged', undefined], ['onColumnsChanged', undefined]],
  'a rename fires onColumnNameChanged alone; an unchanged row index fires nothing');
  assert.deepEqual(df.columns.names(), ['years', 'flag']);
  df.columns.remove('gone');
  assert.deepEqual(df.columns.names(), ['years', 'flag'], 'an unknown name removes nothing');
  subs[0].unsubscribe();
  assert.equal(df.liveSubscriptions(), 3, 'an unsubscribed stream is let go');
  const bytes = new Uint8Array([1, 2]);
  const read = DataFrame.fromByteArray(bytes);
  assert.deepEqual([read.rowCount, read.dart.bytes], [0, bytes], 'the bytes stay behind the handle');
});

test('an entity row answers its names; a user adds an email', () => {
  const alice = new User('alice', {id: 'u1', friendlyName: 'Alice M', email: 'alice@datagrok.ai'});
  assert.deepEqual([alice.id, alice.name, alice.friendlyName, alice.nqName, alice.email],
    ['u1', 'alice', 'Alice M', 'alice', 'alice@datagrok.ai']);
  assert.equal(new Entity('bob').friendlyName, undefined,
    'no fallback here: the platform derives camelCaseToWords(name), and every u2 read goes friendlyName ?? name');
});

test('the shell records every assignment of the current object, and names tables apart', () => {
  const shell = new Shell();
  const ref = {};
  shell.o = ref;
  shell.o = null;
  assert.deepEqual(shell.dart.writes, [ref, null]);
  assert.equal(shell.o, null);
  assert.equal(shell.windows.showContextPanel, false);

  const added = [];
  shell.dart.tableAdded.subscribe(() => added.push(shell.tableNames));
  assert.equal(shell.addTable(demog()).name, 'demog');
  assert.equal(shell.addTable(demog()).name, 'demog (2)');
  assert.deepEqual(added, [['demog'], ['demog', 'demog (2)']]);
  shell.closeTable('demog');
  assert.deepEqual(shell.tableNames, ['demog (2)']);
  assert.doesNotThrow(() => shell.info('x') ?? shell.warning('x') ?? shell.error('x'));
});

test('the loader hook serves the same classes, not a second copy', () => {
  assert.equal(DG.Property, Property);
  assert.equal(DG.Func, Func);
  assert.equal(DG.DataQuery, DataQuery);
  assert.equal(DG.Script, Script);
  assert.equal(DG.TableQuery, TableQuery);
  assert.equal(DG.FileInfo, FileInfo);
  assert.equal(DG.DataFrame, DataFrame);
  assert.ok(grok.shell instanceof Shell);
  assert.ok(new DG.FileInfo('a.csv') instanceof FileInfo, 'so an instanceof in the module under test holds');
});
