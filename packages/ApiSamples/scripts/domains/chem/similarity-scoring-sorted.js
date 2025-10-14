let molecules = await grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv');
let similar = await grok.chem.findSimilar(molecules.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1')    
grok.shell.addTableView(similar);
