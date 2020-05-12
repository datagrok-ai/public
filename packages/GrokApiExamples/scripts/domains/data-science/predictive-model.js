// An example of using predictive models.
//
// https://datagrok.ai/help/plugins/predictive-modeling

grok.loadDataFrame('https://public.datagrok.ai/demo/demog.csv')
    .then(t => ml.applyModel('Demo:PredictSexByBasicDemographics', t)
        .then(t => grok.addTableView(t)));
