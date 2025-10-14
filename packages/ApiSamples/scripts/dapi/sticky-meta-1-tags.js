//name: sticky-meta-1-tags
//tags: demo
//language: javascript

// uniqueName ensures that created sticky meta configuration will be unique.
const uniqueName = (prefix) => 'apisamples-' + prefix + '-' + (Math.random() + 1).toString(36).substring(7);

// In this sample we will attach metadata to experiments
// They will be identified through strings: [experiment 1, experiment 2, ..., experiment 6].
const experiments = DG.Column.fromStrings('experiment', Array.from(Array(6).keys()).map((idx) => 'experiment ' + (idx + 1)));

// Objects that store sticky meta should be identifiable through tags or semtypes.
// We will use tag named "source".
const experimentTag = uniqueName('experiment-tag');
experiments.setTag('source', experimentTag);

// Authors of experiments.
// We will save them as metadata related to experiment.
const authors = DG.DataFrame.fromCsv(
  `name,email,age
   Joey,joey@datagrok.ai,25
   Ross,ross@datagrok.ai,26
   Rachel,rachel@datagrok.ai,24
   Chandler,chandler@datagrok.ai,26
   Monica,monica@datagrok.ai,24
   Phoebe,phoebe@datagrok.ai,26`
);


// Create a schema named Author.
// Author is connected with experiments and has 3 properties of different types
const schema = await grok.dapi.stickyMeta.createSchema(
  uniqueName('Author'),
  [{name: uniqueName('Experiment'), matchBy: 'source=' + experimentTag}],
  [{name: 'name', type: 'string'}, {name: 'email', type: 'string'}, {name: 'age', type: 'int'}]
);

// Save sticky meta values.
await grok.dapi.stickyMeta.setAllValues(schema, experiments, authors);

// Create test data with tagged column of experiment identifiers.
const testDf = DG.DataFrame.fromCsv(`experiment
experiment 1
experiment 3`
);
testDf.columns.byName('experiment').setTag('source', experimentTag);

// Fetch sticky meta related to the test dataframe.
const metaDataframe = await grok.dapi.stickyMeta.getAllValues(schema, testDf.columns.byName('experiment'));

// Show the result
for (var col of metaDataframe.columns)
  testDf.columns.add(col);
grok.shell.addTableView(testDf);