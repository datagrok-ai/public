//name: sticky-meta-tags
//tags: demo
//language: javascript
let t1 = grok.data.testData('demog', 10000);
t1.columns.byName('study').setTag('source', 'study');
let type = DG.EntityType.create('apisamples-study', 'source=study');
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

var newColumn = DG.Column.string('studies', 100);
newColumn.setTag('source', 'study');
var valuesColumn = DG.Column.dateTime('date', 100);
for (var i = 0; i < newColumn.length; i++) {
  newColumn.set(i, 'study ' + i);
  valuesColumn.set(i, '03/' + (i % 30 + 1) + '/2024');
}

var autofilledColumn = t1.columns.addNew('sticky-meta-date', 'datetime');
autofilledColumn.setTag('dbPropertyName', property.name);
autofilledColumn.setTag('dbPropertySchema', schema.name);
autofilledColumn.setTag('dbPropertyReference', 'study');

await grok.dapi.stickyMeta.setAllValues(schema, newColumn, DG.DataFrame.fromColumns([valuesColumn]));
grok.shell.addTableView(t1);
