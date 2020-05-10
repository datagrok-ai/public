let table = grok.testData('random walk', 100000, 2);

measure = function (f) {
    let t0 = performance.now();
    f();
    let t1 = performance.now();
    return t1 - t0;
};

let bitsetIteration = measure(function() {
    for (let i = 0; i < table.rowCount; i++)
        table.filter.get(i);
});

let valuesIteration = measure(function() {
    let col = table.columns.byIndex(0);
    for (let i = 0; i < table.rowCount; i++)
        col.get(i);
});

ui.dialog()
    .add(`bitset: ${bitsetIteration}`)
    .add(`values: ${valuesIteration}`)
    .show();