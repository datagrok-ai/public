/* The dg filter adapters against the platform doubles: the three schema constructors (columns,
   domain table, entity type), the platform value-editor factory, and the BitSet delivery of a
   mask. `DG`/`grok` come from tests/dg-stub.mjs; the `grok.meta` reflection the stub does not
   carry is layered over its `grok` by a local hook, so `forEntityType` runs headless too. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {TextInput, TextArea} from '../src/components/inputs/text-input.js';
import {Filters, FilterError} from '../src/core/filter/index.js';
import {BitSet, DataFrame, Entity, Property} from './platform-doubles.mjs';

const META_HOOK = `
const GROK = 'u2-test:datagrok-api-grok-meta';
export async function resolve(specifier, context, next) {
  if (specifier === 'datagrok-api/grok' && context.parentURL !== GROK)
    return {url: GROK, format: 'module', shortCircuit: true};
  return next(specifier, context);
}
export async function load(url, context, next) {
  if (url === GROK) {
    return {format: 'module', shortCircuit: true,
      source: "export * from 'datagrok-api/grok'; export const meta = {propertiesOf: async () => null};"};
  }
  return next(url, context);
}`;

register('./dg-stub.mjs', import.meta.url);
register(`data:text/javascript,${encodeURIComponent(META_HOOK)}`, import.meta.url);
const grok = await import('datagrok-api/grok');
const {FilterSchemas, frameLike, toBitSet} = await import('../src/dg/filter/index.js');
const {RefInput} = await import('../src/index.js');
const {Editors} = await import('../src/dg/forms/editors.js');
const {inputForProperty} = await import('../src/dg/forms/object-form.js');
const {FilterBuilder} = await import('../src/components/filter/filter-builder.js');

const abort = new AbortController().signal;

function smoke(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      delete grok.dapi.domains;
      delete grok.dapi.groups;
      grok.meta.propertiesOf = async () => null;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function propOf(schema, name) {
  return schema.properties.find((p) => p.name === name);
}

function people() {
  return new DataFrame(
    [{name: 'name', type: 'string'}, {name: 'age', type: 'int'}, {name: 'weight', type: 'double'},
      {name: 'smiles', type: 'string', semType: 'Molecule'}],
    [{name: 'Ann', age: 34, weight: 1.5, smiles: 'CCO'}, {name: 'Bob', age: 28, weight: 2.5, smiles: 'CC'},
      {name: 'Cid', age: 45, weight: 0.5, smiles: 'C'}]);
}

// --- forDataFrame ---

smoke('forDataFrame: one property per column with kind and a molecule renderer; the data sets no bounds', async () => {
  const df = people();
  const schema = FilterSchemas.forDataFrame(df);
  assert.deepEqual(schema.properties.map((p) => p.name), ['name', 'age', 'weight', 'smiles']);
  assert.deepEqual(schema.properties.map((p) => Filters.kindOf(p)), ['string', 'int', 'float', 'string']);
  assert.deepEqual(propOf(schema, 'name'), {name: 'name', type: 'string'}, 'categories are values, not choices');
  assert.deepEqual(propOf(schema, 'age'), {name: 'age', type: 'int'}, 'the frame min/max are not bounds');
  assert.equal(propOf(schema, 'smiles').semType, 'Molecule');
  assert.ok(schema.renderer(propOf(schema, 'smiles')), 'molecule columns render through the sketcher depiction');
  assert.equal(schema.renderer(propOf(schema, 'name')), undefined);
  assert.deepEqual(Filters.validate(Filters.group('and', [Filters.cond('age', '<', 20), Filters.cond('name', '=', 'Dan')]),
    schema, 'dataframe'), [], 'a filter may name what the frame does not hold');

  assert.deepEqual(await schema.values(propOf(schema, 'name'), 'b', abort), [{value: 'Bob', label: 'Bob'}]);
  assert.deepEqual(await schema.values(propOf(schema, 'age'), '', abort), []);
  const reads = [];
  const name = df.columns.byName('name');
  Object.defineProperty(name, 'categories', {get: () => (reads.push(1), ['Ann', 'Bob', 'Cid'])});
  await schema.values(propOf(schema, 'name'), 'a', abort);
  assert.equal(reads.length, 0, 'categories are read once per column, not per keystroke');
  schema.refresh();
  assert.deepEqual(await schema.values(propOf(schema, 'name'), 'a', abort), [{value: 'Ann', label: 'Ann'}]);
  assert.equal(reads.length, 1, 'refresh re-reads them');
  const rows = Array.from({length: 25}, (_, i) => ({code: `c${i}`}));
  const many = FilterSchemas.forDataFrame(new DataFrame([{name: 'code', type: 'string'}], rows));
  assert.equal((await many.values(propOf(many, 'code'), 'c1', abort)).length, 11);
});

smoke('forDataFrame: refresh re-reads the columns — a semantic type detected after the snapshot changes the operators', () => {
  const df = new DataFrame([{name: 'smiles', type: 'string'}], [{smiles: 'CCO'}]);
  const schema = FilterSchemas.forDataFrame(df);
  const before = schema.properties;
  assert.deepEqual(before, [{name: 'smiles', type: 'string'}]);
  assert.equal(schema.renderer(propOf(schema, 'smiles')), undefined);
  df.columns.byName('smiles').dart.semType = 'Molecule';
  assert.equal(schema.properties, before, 'a snapshot until refresh');
  schema.refresh();
  assert.notEqual(schema.properties, before, 'replaced, not mutated');
  assert.equal(propOf(schema, 'smiles').semType, 'Molecule');
  assert.ok(schema.renderer(propOf(schema, 'smiles')));
});

// --- forDomainTable ---

function domainsDouble() {
  const calls = {rowProperties: [], facets: []};
  const tables = {
    'plates.plate': [
      new Property('name', 'string', {friendlyName: 'Name', nullable: false,
        get: (r) => r.name, set: (r, v) => r.name = v}),
      new Property('project_id', 'string', {friendlyName: 'Project', semType: 'Core.projects'}),
    ],
    'Core.projects': [new Property('name', 'string'), new Property('created_on', 'datetime')],
  };
  grok.dapi.domains = {
    registry: {
      rowProperties: async (address) => {
        calls.rowProperties.push(address);
        return tables[address];
      },
    },
    table: (address) => ({
      facets: async (spec) => {
        calls.facets.push([address, spec]);
        const column = spec.facets[0].column;
        return {facets: {v: {categories: column === 'project_id' ?
          [{value: 'p1', display: 'Alpha', total: 3, filtered: 1}] :
          [{value: 'Ann', total: 2, filtered: 2}]}}};
      },
    }),
  };
  return calls;
}

smoke('forDomainTable: registry properties are copied field by field, a dotted semType is the ref', async () => {
  const calls = domainsDouble();
  const schema = await FilterSchemas.forDomainTable('plates.plate');
  assert.deepEqual(calls.rowProperties, ['plates.plate']);

  const source = (await grok.dapi.domains.registry.rowProperties('plates.plate'))[0];
  assert.deepEqual(Object.keys({...source}), ['dart'], 'a spread of the platform property yields only its handle');
  const name = propOf(schema, 'name');
  assert.deepEqual([name.propertyType, name.friendlyName, name.nullable], ['string', 'Name', false]);
  assert.equal(typeof name.get, 'function');
  assert.equal(name.ref, undefined);

  const project = propOf(schema, 'project_id');
  assert.equal(project.ref, 'Core.projects');
  assert.equal(Filters.kindOf(project), 'ref');
  assert.equal(schema.renderer(project).caption({type: 'Core.projects', id: 'p1', name: 'Alpha'}), 'Alpha');
  assert.equal(schema.renderer(project).caption({type: 'Core.projects', id: 'p1'}), 'p1');
});

smoke('forDomainTable: values go through the categories facet, refs come back as FilterRef', async () => {
  const calls = domainsDouble();
  const schema = await FilterSchemas.forDomainTable('plates.plate');

  const refs = await schema.values(propOf(schema, 'project_id'), 'al', abort);
  assert.deepEqual(calls.facets, [['plates.plate', {facets: [
    {id: 'v', kind: 'categories', column: 'project_id', search: 'al', limit: 50}]}]]);
  assert.deepEqual(refs, [{value: {type: 'Core.projects', id: 'p1', name: 'Alpha'}, label: 'Alpha', count: 3}]);

  const names = await schema.values(propOf(schema, 'name'), '', abort);
  assert.deepEqual(names, [{value: 'Ann', label: 'Ann', count: 2}]);
});

smoke('forDomainTable: resolveRef loads the target table', async () => {
  const calls = domainsDouble();
  const schema = await FilterSchemas.forDomainTable('plates.plate');
  const project = propOf(schema, 'project_id');

  const target = await schema.resolveRef(project);
  assert.deepEqual(calls.rowProperties, ['plates.plate', 'Core.projects']);
  assert.deepEqual(target.properties.map((p) => p.name), ['name', 'created_on']);
  assert.equal(Filters.kindOf(propOf(target, 'created_on')), 'datetime');
});

// --- forEntityType ---

function metaDouble() {
  const calls = {propertiesOf: [], groups: []};
  const info = (name, type, extra = {}) => ({name, type, semType: null, friendlyName: name,
    description: null, refType: null, relationKind: null, ...extra});
  const types = {
    User: [info('login', 'string', {friendlyName: 'Login'}),
      info('group', 'object', {refType: 'Group', description: 'Personal group'})],
    Group: [info('friendlyName', 'string')],
  };
  grok.meta.propertiesOf = async (type, options) => {
    calls.propertiesOf.push([type, options]);
    return types[type] ?? null;
  };
  grok.dapi.groups = {
    list: async (options) => {
      calls.groups.push(options);
      return [new Entity('Developers', {id: 'g1', friendlyName: 'Developers'})];
    },
  };
  return calls;
}

smoke('forEntityType: filterable properties, refType as the ref, ref values from the dapi collection', async () => {
  const calls = metaDouble();
  const schema = await FilterSchemas.forEntityType('User');
  assert.deepEqual(calls.propertiesOf, [['User', {filterable: true}]]);
  assert.deepEqual(schema.properties.map((p) => p.name), ['login', 'group']);
  assert.deepEqual(propOf(schema, 'login'), {name: 'login', type: 'string', friendlyName: 'Login'});
  const group = propOf(schema, 'group');
  assert.deepEqual(group, {name: 'group', type: 'object', friendlyName: 'group', description: 'Personal group',
    ref: 'Group'});
  assert.equal(Filters.kindOf(group), 'ref');

  assert.deepEqual(await schema.values(group, 'dev', abort),
    [{value: {type: 'Group', id: 'g1', name: 'Developers'}, label: 'Developers'}]);
  assert.deepEqual(calls.groups, [{pageSize: 20, filter: 'dev'}]);
  assert.deepEqual(await schema.values(propOf(schema, 'login'), 'x', abort), [], 'no collection behind a plain property');

  const target = await schema.resolveRef(group);
  assert.deepEqual(calls.propertiesOf[1], ['Group', {filterable: true}]);
  assert.deepEqual(target.properties.map((p) => p.name), ['friendlyName']);

  const renderer = schema.renderer(group);
  assert.equal(schema.renderer(propOf(schema, 'login')), undefined);
  assert.equal(renderer.caption({type: 'Group', id: 'g1', name: 'Developers'}), 'Developers');
  assert.equal(renderer.caption({type: 'Group', id: 'g1'}), 'g1');
  assert.equal(renderer.listItem({type: 'Group', id: 'g1', name: 'Developers'}).textContent, 'Developers',
    'a ref row is its caption, never a handler');
  const {ObjectHandler} = await import('datagrok-api/dg');
  const handler = {isApplicable: (x) => x instanceof Entity, getCaption: (x) => `entity ${x.friendlyName}`};
  ObjectHandler.register(handler);
  try {
    assert.equal(renderer.caption(new Entity('Developers', {id: 'g1', friendlyName: 'Developers'})), 'entity Developers',
      'an entity value goes through its handler');
  } finally {
    ObjectHandler.registered.splice(ObjectHandler.registered.indexOf(handler), 1);
  }
  assert.equal(typeof renderer.icon, 'function');
  assert.equal(typeof renderer.tooltip, 'function');
});

smoke('forEntityType: a type without a catalog is an error, not an empty schema', async () => {
  metaDouble();
  await assert.rejects(FilterSchemas.forEntityType('Nope'), FilterError);
});

// --- the one editor resolver ---

/** The value editor the builder mounts for the one condition of `schema`'s first property. */
function rowEditor(schema, value) {
  const prop = schema.properties[0];
  const fb = new FilterBuilder({schema, value: Filters.group('and', [Filters.cond(prop.name, '=', value)])});
  document.body.append(fb.root);
  return {fb, editor: fb.root.querySelector('[data-u2-part="value"]')};
}

function testRule() {
  return Editors.register({
    match: (p) => p.semType === 'Test',
    create: (p, options) => {
      const input = new TextInput(options);
      input.root.dataset.test = 'registered';
      return input;
    },
  });
}

smoke('editors: a registered Editors rule is what the builder mounts', async () => {
  const unregister = testRule();
  try {
    const {fb, editor} = rowEditor(Filters.schema([{name: 'smiles', type: 'string', semType: 'Test'}]), 'x');
    assert.equal(editor.dataset.test, 'registered');
    fb.dispose();
  } finally {
    unregister();
  }
});

smoke('editors: a semType-ruled property gets the same editor in a form and in a filter row, hints notwithstanding',
  async () => {
    const unregister = testRule();
    try {
      const schema = Filters.schema([{name: 'smiles', type: 'string', semType: 'Test', inputType: 'TextArea'}]);
      const form = inputForProperty(schema.properties[0], {assumeWritable: true});
      assert.equal(form.root.dataset.test, 'registered', 'the rule wins over the inputType hint in the form');
      form.dispose();
      const {fb, editor} = rowEditor(schema, 'x');
      assert.equal(editor.dataset.test, 'registered', 'and in the row');
      fb.dispose();
    } finally {
      unregister();
    }
  });

smoke('editors: inputType/editor hints resolve after the rules; plain kinds return null', async () => {
  const schema = Filters.schema([{name: 'notes', type: 'string', inputType: 'TextArea'},
    {name: 'name', type: 'string'}, {name: 'age', type: 'int'}]);
  const notes = Editors.resolve(propOf(schema, 'notes'), {});
  assert.ok(notes instanceof TextArea);
  assert.equal(notes.enabled, true, 'a filter value is always writable');
  notes.dispose();
  assert.equal(Editors.resolve(propOf(schema, 'name'), {}), null);
  assert.equal(Editors.resolve(propOf(schema, 'age'), {}), null);
  assert.equal(FilterBuilder.defaultEditors(propOf(schema, 'age'), {}), null, 'the builder default is the resolver');
});

smoke('editors: a ref is a type-ahead over schema.values, whose value is the FilterRef', async () => {
  const items = [{value: {type: 'Group', id: 'g1', name: 'Developers'}, label: 'Developers'},
    {value: {type: 'Group', id: 'g2', name: 'Testers'}, label: 'Testers'}];
  const schema = {
    properties: [{name: 'group', type: 'object', ref: 'Group'}],
    values: async (prop, query) => items.filter((i) => i.label.toLowerCase().includes(query.toLowerCase())),
  };
  const input = new RefInput({prop: schema.properties[0], schema, debounceMs: 0});
  try {
    document.body.append(input.root);
    const box = input.root.querySelector('[data-u2="typeahead"] input');
    assert.ok(box);

    input.value.value = {type: 'Group', id: 'g2', name: 'Testers'};
    await flush();
    assert.equal(box.value, 'Testers', 'a written value shows its name');

    box.focus();
    box.value = 'dev';
    fire(box, 'input');
    await flush();
    await flush();
    fire(box, 'keydown', {key: 'Enter'});
    await flush();
    assert.deepEqual(input.value.value, {type: 'Group', id: 'g1', name: 'Developers'});
  } finally {
    input.dispose();
  }
  const withValues = rowEditor(schema, {type: 'Group', id: 'g2', name: 'Testers'});
  assert.equal(withValues.editor.dataset.u2, 'ref-input', 'the row mounts the type-ahead over the values');
  assert.equal(withValues.editor.querySelector('input').value, 'Testers');
  withValues.fb.dispose();
  const noValues = rowEditor({properties: schema.properties}, {type: 'Group', id: 'g2', name: 'Testers'});
  assert.equal(noValues.editor.dataset.u2, 'text-input', 'without a value source the core text editor applies');
  noValues.fb.dispose();
});

// --- toBitSet ---

smoke('toBitSet: one fromBytes over the mask bits, rowCount long', async () => {
  const df = people();
  const calls = [];
  const fromBytes = BitSet.fromBytes;
  BitSet.fromBytes = (buffer, length) => {
    calls.push([buffer, length]);
    return fromBytes(buffer, length);
  };
  try {
    const bitset = await toBitSet(df, Filters.group('and', [
      Filters.cond('age', '>', 30), Filters.cond('name', '!=', 'Cid')]));
    assert.equal(calls.length, 1);
    assert.ok(calls[0][0] instanceof ArrayBuffer);
    assert.equal(calls[0][1], 3);
    assert.equal(bitset.length, 3);
    assert.deepEqual([0, 1, 2].map((i) => bitset.get(i)), [true, false, false]);
  } finally {
    BitSet.fromBytes = fromBytes;
  }
  assert.equal(frameLike(df).column('nope'), null);
  assert.equal(frameLike(df).column('age').length, 3);
});
