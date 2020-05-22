// An example of using missing values imputation.
//
// https://datagrok.ai/help/dialogs/missing-values-imputation

grok.data.loadDataFrame('https://public.datagrok.ai/demo/demog.csv')
    .then(t => grok.ml.missingValuesImputation(t, ['age', 'height', 'weight'], ['age', 'height', 'weight'], 5)
        .then(t => grok.shell.addTableView(t)));
