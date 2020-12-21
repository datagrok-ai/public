// An example of getting Maximum Common Substructure (MCS).

grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv')
  .then(t => grok.chem.mcs(t.col('smiles'))
    .then(mcs => grok.shell.info(mcs)));

