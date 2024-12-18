// An example of getting R-Groups
// https://datagrok.ai/help/domains/chem/r-group-analysis

let table = await grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv');
let groupedTable = await grok.chem.rGroup(table, 'smiles', 'O=C1Nc2ccccc2C(C2CCCCC2)=NC1');
grok.shell.addTableView(groupedTable);