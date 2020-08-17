// https://datagrok.ai/help/viewers/histogram

let view = grok.shell.addTableView(grok.data.demo.demog());

view.histogram({
    value: 'age'
});
