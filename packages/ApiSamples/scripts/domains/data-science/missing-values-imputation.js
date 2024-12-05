// An example of using missing values imputation.
//
// https://datagrok.ai/help/transform/missing-values-imputation

let t1 = await grok.data.loadTable('https://public.datagrok.ai/demo/demog.csv');
let t2 = await grok.ml.missingValuesImputation(t1, ['age', 'height', 'weight'], ['age', 'height', 'weight'], 5);
grok.shell.addTableView(t2);
