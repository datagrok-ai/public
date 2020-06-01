// Assigning custom category colors

let view = grok.shell.addTableView(grok.data.testData('demog', 5000));

view.grid.col('sex').categoryColors = {
    'M' : 0xFF0000FF,
    'F' : 0xFF800080
};
