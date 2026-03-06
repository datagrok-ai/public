const t = await grok.data.files.openTable('System:AppData/ApiSamples/molblocks.csv');
grok.shell.addTableView(t);
console.log(await t.toCsvEx({moleculesAsSmiles: true}));