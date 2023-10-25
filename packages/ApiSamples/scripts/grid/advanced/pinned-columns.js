const df = grok.data.demo.demog();
const tv = grok.shell.addTableView(df);

const pinnedColumn = tv.grid.columns.byName('subj');

// calling the function from the PowerGrid package
await grok.functions.call('PowerGrid:addPinnedColumn', {gridCol: pinnedColumn});