// Gets columns by type or tags

let demog = grok.data.demo.demog();

for (let column of demog.columns.categorical)
    grok.shell.info(column.name);

for (let column of demog.columns.numerical)
    grok.shell.info(column.name);


demog.getCol('race').setTag('tag1', '["Asian", "Black"]');
demog.getCol('age').setTag('tag2', 'something');

for (let column of demog.columns.byTags({'tag1': '["Asian", "Black"]'}))
    grok.shell.info(column.name);


for (let column of demog.columns.byTags({'tag2': undefined}))
    grok.shell.info(column.name);