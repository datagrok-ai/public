//name: sticky-meta-semtype
//tags: demo
//language: javascript
let t1 = grok.data.testData('molecules', 10000);
console.log(t1);
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



var valuesColumn = DG.Column.dateTime('date', 1000);
valuesColumn.set(0, '03/20/2024');


await grok.dapi.stickyMeta.setAllValues(schema, t1.columns.byName('smiles'), DG.DataFrame.fromColumns([valuesColumn]));

var autofilledColumn = t1.columns.addNew('sticky-meta-date', 'datetime');
autofilledColumn.setTag('dbPropertyName', property.name);
autofilledColumn.setTag('dbPropertySchema', schema.name);
autofilledColumn.setTag('dbPropertyReference', 'smiles');

