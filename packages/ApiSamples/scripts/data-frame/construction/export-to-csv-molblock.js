const t = await grok.data.files.openTable('System:AppData/Chem/molblocks.csv');
grok.shell.addTableView(t);
console.log(await t.toCsvEx({moleculesAsSmiles: true}));