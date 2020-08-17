// https://datagrok.ai/help/viewers/grid

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.setOptions({
    colHeaderHeight: 80,
});
