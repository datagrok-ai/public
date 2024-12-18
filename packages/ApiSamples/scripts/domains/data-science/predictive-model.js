// An example of using predictive models.
// https://datagrok.ai/help/learn/predictive-modeling

let t1 = await grok.data.loadTable('https://public.datagrok.ai/demo/demog.csv');
let t2 = await grok.ml.applyModel('Samples:Models:PredictSexByBasicDemographics', t1);
grok.shell.addTableView(t2);