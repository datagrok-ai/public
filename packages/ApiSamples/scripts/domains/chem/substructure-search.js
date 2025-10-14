// An example of using substructure search

let t = await grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv');
let bs = await grok.chem.searchSubstructure(t.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');

t.filter.copyFrom(bs);
grok.shell.addTableView(t);