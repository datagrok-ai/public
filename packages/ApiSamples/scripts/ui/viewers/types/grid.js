// https://datagrok.ai/help/visualize/viewers/grid

let view = grok.shell.addTableView(grok.data.demo.demog());

view.grid.setOptions({
    colHeaderHeight: 80,
});
