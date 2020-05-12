// An example of getting Maximum Common Substructure (MCS).

grok.loadDataFrame('https://public.datagrok.ai/demo/sar_small.csv')
    .then(t => chem.mcs(t.col('smiles'))
        .then(mcs => grok.balloon.info(mcs)));

