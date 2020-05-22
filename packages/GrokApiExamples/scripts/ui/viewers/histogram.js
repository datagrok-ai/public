// https://datagrok.ai/help/viewers/histogram

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));

view.histogram({
    value: 'age'
});
