// Customize TableView toolbox elements

let t = grok.data.demo.demog();
let view = grok.shell.addTableView(t);

let acc = view.toolboxPage.accordion;
acc.addPane('Demo', () => ui.divText(`Cells count: ${t.rowCount * t.columns.length}`), true, acc.panes[0]);
