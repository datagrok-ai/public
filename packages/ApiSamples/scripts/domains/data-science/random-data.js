// An example of using random data generator.
//
// https://datagrok.ai/help/transform/random-data

let t = grok.data.testData('demog', 1000);

async function generate() {
  await grok.ml.randomData(t, 'normal', {sd: 3.0, mean: 1.0}, 77);
  await grok.ml.randomData(t, 'uniform', {min: 0.0, max: 1.0}, 77);
  await grok.ml.randomData(t, 'binomial', {size: 100, prob: 0.7}, 77);
}

await generate();
let view = grok.shell.addTableView(t);
view.histogram({value: 'normal'});
view.histogram({value: 'uniform'});
view.histogram({value: 'binomial'});