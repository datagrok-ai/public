import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

//@ts-ignore
import {after, assure, before, category, expect, test} from '@datagrok-libraries/test/src/test';

category('Dapi: sticky meta', () => { 
  let schema: _DG.Schema;
  let entityType: _DG.EntityType;
  let property: _DG.EntityProperty;
  let handle = 'api-test';
  let value = 'val';
  let col: _DG.Column<String>;
  let valDf: _DG.DataFrame;

  before(async () => {
    let schemas = await grok.dapi.stickyMeta.getSchemas();
    for (let i = 0; i < schemas.length; i++) {
      if (schemas[i].name == "default") {
        schema = schemas[i];
        break;
      }
    }
    if (schema == undefined) {
      schema = DG.Schema.create('default');
      await grok.dapi.stickyMeta.saveSchema(schema);
    }
    for (let i = 0; i < schema.entityTypes.length; i++) {
      if (schema.entityTypes[i].matching == 'semtype=test') {
        entityType = schema.entityTypes[i];
        break;
      }
    }
    if (entityType == undefined) {
      entityType = DG.EntityType.create('apiTest', 'semtype=test');
      const et = schema.entityTypes;
      et.push(entityType);
      schema.entityTypes = et;
      await grok.dapi.stickyMeta.saveSchema(schema);
    }
    for (let i = 0; i < schema.properties.length; i++) {
      if (schema.properties[i].name == 'apiTest') {
        property = schema.properties[i];
      }
    }
    if (property == undefined) {
      property = DG.EntityProperty.create('apiTest', 'string');
      const props = schema.properties;
      props.push(property);
      schema.properties = props;
      await grok.dapi.stickyMeta.saveSchema(schema);
    }
    assure.notNull(schema);
    assure.notNull(entityType);
    assure.notNull(property);
    col = DG.Column.string('key', 1);
    col.semType = 'test';
    col.set(0, handle);
    valDf = DG.DataFrame.create(1);
    valDf.columns.addNewString(property.name);
    valDf.set(property.name, 0, value);
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

// Sticky meta keyed by a molecule: the platform canonicalizes the SMILES (Chem's canonicalizer) on write
// and on read, so the values follow the molecule, not the string or the table that holds it. Read the
// way the grid's "Sticky meta" columns read them (AddStickyPropertiesColumns).
category('Dapi: sticky meta: molecules', () => {
  const name = `apitests-sm-${DG.Utils.randomString(6)}`;
  const typeName = `${name}-type`;
  let schema: _DG.Schema;
  const idOf = (x: {dart: any}) => new DG.Entity(x.dart).id;

  function molecules(smiles: string[]): _DG.DataFrame {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', smiles)]);
    df.col('smiles')!.semType = DG.SEMTYPE.MOLECULE;
    df.col('smiles')!.meta.units = DG.UNITS.Molecule.SMILES;
    return df;
  }

  async function valuesOf(df: _DG.DataFrame): Promise<string> {
    const before = df.columns.names();
    await grok.functions.call('AddStickyPropertiesColumns', {table: df, column: df.col('smiles'),
      entityTypeName: typeName, propertyIds: schema.properties.map(idOf)});
    const added = df.columns.toList().filter((c) => !before.includes(c.name));
    return added.map((c) => `${c.name}: ${c.toList().join(', ')}`).sort().join('; ');
  }

  const written = 'notes: excellent, good; rating: 5, 4';

  before(async () => {
    await grok.dapi.stickyMeta.createSchema(name, [{name: typeName, matchBy: 'semtype=molecule'}],
      [{name: 'rating', type: 'int'}, {name: 'notes', type: 'string'}]);
    schema = (await grok.dapi.stickyMeta.getSchemas()).find((s) => s.name === name)!;
    await grok.dapi.stickyMeta.setAllValues(schema, molecules(['CC(=O)Oc1ccccc1C(=O)O', 'Oc1ccccc1']).col('smiles')!,
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
    expect(await valuesOf(molecules(['CC(=O)Oc1ccccc1C(=O)O', 'Oc1ccccc1']).clone()), written);
  });

  test('A table exported to CSV and read back shows the values', async () => {
    const df = DG.DataFrame.fromCsv(molecules(['CC(=O)Oc1ccccc1C(=O)O', 'Oc1ccccc1']).toCsv());
    await grok.data.detectSemanticTypes(df);
    expect(df.col('smiles')!.semType, DG.SEMTYPE.MOLECULE);
    expect(await valuesOf(df), written);
  });

  test('Another SMILES of the same molecule shows its values, another molecule none', async () => {
    expect(await valuesOf(molecules(['OC(=O)c1ccccc1OC(C)=O', 'c1ccc(O)cc1'])), written);
    expect(await valuesOf(molecules(['CCCCCCCC'])), 'notes: ; rating: ');
  });
}, {owner: 'aparamonov@datagrok.ai'});
