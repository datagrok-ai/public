// Automatic discovery of renderers based on the semantic type

let wells = grok.data.demo.wells();
wells.col('pos').semType = 'demo_plate';
grok.shell.addTableView(wells);