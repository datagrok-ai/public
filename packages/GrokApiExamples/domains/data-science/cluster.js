// An example of using clustering.
//
// https://datagrok.ai/help/dialogs/cluster-data

gr.loadDataFrame('/demo/xclara.csv')
    .then(t => ml.cluster(t, ["V1", "V2"], 3)
        .then(function (t) {
            let view = gr.addTableView(t);
            view.scatterPlot({
                x: 'V1',
                y: 'V2',
                color: 'clusters',
            });
        }));
