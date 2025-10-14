// An example of getting Maximum Common Substructure (MCS).

let t = await grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv');
let mcs = await grok.chem.mcs(t, 'smiles'); 
grok.shell.info(mcs)