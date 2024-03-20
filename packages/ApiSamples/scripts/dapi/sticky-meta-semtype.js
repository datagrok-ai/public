//name: sticky-meta-semtype
//tags: demo
//language: javascript
let t1 = grok.data.testData('molecules', 1000);
let type = DG.EntityType.create('apisamples-molecule', 'semtype=molecule');
let property = DG.EntityProperty.create('date', 'datetime');
let schema = DG.Schema.create('Date');
let schemas = await grok.dapi.stickyMeta.getSchemas();
for (var index in schemas) {
  var refSchema = schemas[index];

  if (schema.name === refSchema.name) {
    schema = refSchema;
  }
}
var hasProperty = false;
var hasType = false;
for (var index in schema.properties) {
  var refProperty =  schema.properties[index];
  if (refProperty.name === property.name) {
    hasProperty = true;
  }
}
for (var index in schema.entityTypes) {
  var refType =  schema.entityTypes[index];
  if (refType.name === type.name) {
    hasType = true;
  }
}
if (!hasProperty) {
	schema.properties = schema.properties.concat(property);
}
if (!hasType) {
 	schema.entityTypes = schema.entityTypes.concat(type);
}


await grok.dapi.stickyMeta.saveSchema(schema);
var startView = grok.shell.addTableView(t1);

// for fetching Chem's detectors and functions
setTimeout(async () => {
  var valuesColumn = DG.Column.dateTime('date', 1000);
  valuesColumn.set(0, '03/20/2024');
  
  await grok.dapi.stickyMeta.setAllValues(schema, t1.columns.byName('smiles'), DG.DataFrame.fromColumns([valuesColumn]));
}, 1000);

