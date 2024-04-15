//name: sticky-meta-tags
//tags: demo
//language: javascript

// uniqueName ensures that sample creates new configuration
var uniqueName = (prefix) => 'apisamles-' + prefix + '-' + (Math.random() + 1).toString(36).substring(7);

// Authors of experiments
var authors = DG.DataFrame.fromCsv(
  `name,email,age
   Joey,joey@datagrok.ai,25
   Ross,ross@datagrok.ai,26
   Rachel,rachel@datagrok.ai,24
   Chandler,chandler@datagrok.ai,26
   Monica,monica@datagrok.ai,24
   Phoebe,phoebe@datagrok.ai,26`
);

// Experiments: [experiment 1, experiment 2, ..., experiment 6]
var experiments = DG.Column.fromStrings('experiment', Array.from(Array(6).keys()).map((idx) => 'experiment ' + (idx + 1)));

// For datagrok, objects should be identifiable through tags or semtypes
var experimentTag = uniqueName('experiment-tag');
experiments.setTag('source', experimentTag);

// Create a schema named Author, connected with experiments, and having 3 parameters of different types
var schema = await grok.dapi.stickyMeta.createSchema(
  uniqueName('Author'),
  [{name: uniqueName('Experiment'), matchBy: 'source=' + experimentTag}],
  [{name: 'name', type: 'string'}, {name: 'email', type: 'string'}, {name: 'age', type: 'int'}]
);

// Save sticky meta values
await grok.dapi.stickyMeta.setAllValues(schema, experiments, authors);

// Create test data with tagged experiments.
var testDf = DG.DataFrame.fromCsv(
`experiment
experiment 1
experiment 3`
);
testDf.columns.byName('experiment').setTag('source', experimentTag);

// Enrich dataframe with sticky meta columns
var newData = await grok.dapi.stickyMeta.getAllValues(schema, testDf.columns.byName('experiment'));
for (var i = 0; i < newData.columns.length; i++) {
	testDf.columns.add(newData.columns.byIndex(i));
}

// Open enriched data.
grok.shell.addTableView(testDf);
