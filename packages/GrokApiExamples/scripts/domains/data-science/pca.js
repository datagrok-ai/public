// An example of using PCA (Principal Component Analysis).
//
// https://datagrok.ai/help/dialogs/pca

grok.data.loadTable('https://public.datagrok.ai/demo/cars.csv')
    .then(t => grok.ml.pca(t, ["wheel.base", "length", "width", "height", "city.mpg", "price"], 2, true, true)
        .then(t => grok.shell.addTableView(t)));
