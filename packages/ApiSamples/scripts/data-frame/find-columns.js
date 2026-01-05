// Gets columns by type or tags

let demog = grok.data.demo.demog();

for (let column of demog.columns.categorical)
  grok.shell.info(column.name);

for (let column of demog.columns.numerical)
  grok.shell.info(column.name);

demog.getCol('race')
  .setTag('tag1', 'value1')
  .setTag('tag2', 'value2')
  .setTag('tag', 'value3');

// searching for multiple tags at once
for (let column of demog.columns.byTags({'tag1': 'value1', 'tag2': 'value2'}))
  grok.shell.info(column.name);

// undefined or null means that any value passes
for (let column of demog.columns.byTags({'tag2': undefined}))
  grok.shell.info(column.name);