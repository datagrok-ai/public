// An example of using PCA (Principal Component Analysis).
//
// https://datagrok.ai/help/explore/pca

let t1 = await grok.data.loadTable('https://public.datagrok.ai/demo/cars.csv');
let t2 = await  grok.ml.pca(t1, ["wheel.base", "length", "width", "height", "city.mpg", "price"], 2, true, true);
grok.shell.addTableView(t2);