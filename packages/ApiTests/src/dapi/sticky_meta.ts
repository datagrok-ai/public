import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

//@ts-ignore
import {after, assure, before, category, expect, test} from '@datagrok-libraries/test/src/test';

const idOf = (x: {dart: any}) => new DG.Entity(x.dart).id;

category('Dapi: sticky meta', () => {
  const name = `apitests-sm-${DG.Utils.randomString(6)}`;
  const semType = `${name}-semtype`;
  let schema: _DG.Schema;
  let entityType: _DG.EntityType;
  let property: _DG.EntityProperty;
  let handle = 'api-test';
  let value = 'val';
  let col: _DG.Column<String>;
  let valDf: _DG.DataFrame;

  before(async () => {
    schema = DG.Schema.create(name);
    entityType = DG.EntityType.create(`${name}-type`, `semtype=${semType}`);
    property = DG.EntityProperty.create('apiTest', 'string');
    schema.entityTypes = [entityType];
    schema.properties = [property];
    await grok.dapi.stickyMeta.saveSchema(schema);
    schema = (await grok.dapi.stickyMeta.getSchemas()).find((s) => s.name === name)!;
    assure.notNull(schema);
    col = DG.Column.string('key', 1);
    col.semType = semType;
    col.set(0, handle);
    valDf = DG.DataFrame.create(1);
    valDf.columns.addNewString(property.name);
    valDf.set(property.name, 0, value);
  });

  after(async () => {
    if (schema)
      await grok.dapi.stickyMeta.deleteSchema(idOf(schema));
    for (const t of (await grok.dapi.entityTypes.list()).filter((t) => t.name === entityType.name))
      await grok.dapi.entityTypes.delete(t);
  });

  test('save values', async () => {
    await grok.dapi.stickyMeta.setAllValues(schema, col, valDf);
  }, {stressTest: true});

  test('get values', async () => {
    const col = DG.Column.string('key', 1);
    col.set(0, handle);
    const df = await grok.dapi.stickyMeta.getAllValues(schema, col);
    expect(df.columns.length, schema.properties.length);  
    expect(df.col(property.name)?.name ?? "", property.name);
    expect(df.col(property.name)?.get(0), value);
  });

  test('list entity types', async () => {
    const types = await grok.dapi.entityTypes.list();
    expect(types.length > 0, true);
    for (const t of types)
      assure.notNull(t.name, 'entity type name');
    expect(types.some((t) => t.name.toLowerCase() == entityType.name.toLowerCase()), true);
  });
}, {
  owner: 'aparamonov@datagrok.ai'
});

// Sticky meta keyed by a canonical form: the platform runs the key through the canonicalizer registered
// for the entity type's semantic type, on write and on read, so the values follow the key, not its
// spelling or the table that holds it. Read the way the grid's "Sticky meta" columns read them
// (AddStickyPropertiesColumns). The semantic type, its detector (detectors.js) and its canonicalizer
// (ApiTests:canonicalizeApiTestsKey: trimmed, lower case) are this package's own, as Chem's are for molecules.
category('Dapi: sticky meta: canonical keys', () => {
  const name = `apitests-sm-${DG.Utils.randomString(6)}`;
  const typeName = `${name}-type`;
  const semType = 'apitests-key';
  const keyColumn = 'apitests_key';
  let schema: _DG.Schema;

  function keys(values: string[]): _DG.DataFrame {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(keyColumn, values)]);
    df.col(keyColumn)!.semType = semType;
    return df;
  }

  async function valuesOf(df: _DG.DataFrame): Promise<string> {
    const before = df.columns.names();
    await grok.functions.call('AddStickyPropertiesColumns', {table: df, column: df.col(keyColumn),
      entityTypeName: typeName, propertyIds: schema.properties.map(idOf)});
    const added = df.columns.toList().filter((c) => !before.includes(c.name));
    return added.map((c) => `${c.name}: ${c.toList().join(', ')}`).sort().join('; ');
  }

  const written = 'notes: excellent, good; rating: 5, 4';

  before(async () => {
    expect(await grok.functions.call('ApiTests:canonicalizeApiTestsKey', {key: ' Compound-A '}), 'compound-a');
    await grok.dapi.stickyMeta.createSchema(name, [{name: typeName, matchBy: `semtype=${semType}`}],
      [{name: 'rating', type: 'int'}, {name: 'notes', type: 'string'}]);
    schema = (await grok.dapi.stickyMeta.getSchemas()).find((s) => s.name === name)!;
    await grok.dapi.stickyMeta.setAllValues(schema, keys(['Compound-A', 'Compound-B']).col(keyColumn)!,
      DG.DataFrame.fromColumns([DG.Column.fromList('int', 'rating', [5, 4]),
        DG.Column.fromList('string', 'notes', ['excellent', 'good'])]));
  });

  after(async () => {
    if (schema)
      await grok.dapi.stickyMeta.deleteSchema(idOf(schema));
    for (const t of (await grok.dapi.entityTypes.list()).filter((t) => t.name === typeName))
      await grok.dapi.entityTypes.delete(t);
  });

  test('A clone of the table shows the values', async () => {
    expect(await valuesOf(keys(['Compound-A', 'Compound-B']).clone()), written);
  });

  test('A table exported to CSV and read back shows the values', async () => {
    const df = DG.DataFrame.fromCsv(keys(['Compound-A', 'Compound-B']).toCsv());
    await grok.data.detectSemanticTypes(df);
    expect(df.col(keyColumn)!.semType, semType);
    expect(await valuesOf(df), written);
  });

  test('Another spelling of the same key shows its values, another key none', async () => {
    expect(await valuesOf(keys([' COMPOUND-a', 'compound-b '])), written);
    expect(await valuesOf(keys(['Compound-C'])), 'notes: ; rating: ');
  });
}, {owner: 'aparamonov@datagrok.ai'});
