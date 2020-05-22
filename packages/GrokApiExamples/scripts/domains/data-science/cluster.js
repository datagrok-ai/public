// An example of using clustering.
//
// https://datagrok.ai/help/dialogs/cluster-data

grok.data.loadDataFrame('https://public.datagrok.ai/demo/xclara.csv')
    .then(t => grok.ml.cluster(t, ["V1", "V2"], 3)
        .then((t) => {
            let view = grok.shell.addTableView(t);
            view.scatterPlot({
                x: 'V1',
                y: 'V2',
                color: 'clusters',
            });
        }));
