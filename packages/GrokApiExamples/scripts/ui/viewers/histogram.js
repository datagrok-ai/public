// https://datagrok.ai/help/viewers/histogram

let view = grok.addTableView(grok.testData('demog', 5000));

view.histogram({
    value: 'age'
});