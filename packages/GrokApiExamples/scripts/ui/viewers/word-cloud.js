// https://datagrok.ai/help/viewers/word-cloud

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));

view.wordCloud();
