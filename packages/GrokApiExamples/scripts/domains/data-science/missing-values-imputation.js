// An example of using missing values imputation.
//
// https://datagrok.ai/help/dialogs/missing-values-imputation

grok.loadDataFrame('https://public.datagrok.ai/demo/demog.csv')
    .then(t => ml.missingValuesImputation(t, ['age', 'height', 'weight'], ['age', 'height', 'weight'], 5)
        .then(t => grok.addTableView(t)));
