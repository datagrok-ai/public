// An example of getting Maximum Common Substructure (MCS).

gr.loadDataFrame('/demo/sar_small.csv')
    .then(t => chem.mcs(t.col('smiles'))
        .then(mcs => gr.balloon.info(mcs)));

