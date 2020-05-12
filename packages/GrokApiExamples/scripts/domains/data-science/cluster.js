// An example of using clustering.
//
// https://datagrok.ai/help/dialogs/cluster-data

grok.loadDataFrame('https://public.datagrok.ai/demo/xclara.csv')
    .then(t => ml.cluster(t, ["V1", "V2"], 3)
        .then(function (t) {
            let view = grok.addTableView(t);
            view.scatterPlot({
                x: 'V1',
                y: 'V2',
                color: 'clusters',
            });
        }));
