// https://datagrok.ai/help/visualize/viewers/word-cloud

let view = grok.shell.addTableView(grok.data.demo.demog());

view.addViewer(DG.VIEWER.WORD_CLOUD, {column:'site'});
