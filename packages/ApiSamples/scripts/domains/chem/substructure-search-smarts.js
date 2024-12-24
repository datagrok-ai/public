// An example of using substructure search with SMARTS pattern

let t = await grok.data.getDemoTable('chem/smiles.csv');
let bs = await grok.chem.searchSubstructure(t.col('canonical_smiles'), '[!#6&!#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1')
t.filter.copyFrom(bs);
grok.shell.addTableView(t);