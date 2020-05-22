// An example of using predictive models.
//
// https://datagrok.ai/help/plugins/predictive-modeling

grok.data.loadDataFrame('https://public.datagrok.ai/demo/demog.csv')
    .then(t => grok.ml.applyModel('Demo:PredictSexByBasicDemographics', t)
        .then(t => grok.shell.addTableView(t)));
