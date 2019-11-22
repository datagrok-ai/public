// An example of using predictive models.
//
// https://datagrok.ai/help/plugins/predictive-modeling

gr.loadDataFrame('/demo/demog.csv')
    .then(t => ml.applyModel('Demo:PredictSexByBasicDemographics', t)
        .then(t => gr.addTableView(t)));
