// An example of using random data generator.
//
// https://datagrok.ai/help/dialogs/random-data

var t = grok.testData('demog', 1000);

async function generate() {
    await ml.randomData(t, 'normal', {sd: 3.0, mean: 1.0}, 77);
    await ml.randomData(t, 'uniform', {min: 0.0, max: 1.0}, 77);
    await ml.randomData(t, 'binomial', {size: 100, prob: 0.7}, 77);
}

generate().then(function () {
    let view = grok.addTableView(t);
    view.histogram({value: 'normal'});
    view.histogram({value: 'uniform'});
    view.histogram({value: 'binomial'});
});
