// An example of getting R-Groups.
//
// https://datagrok.ai/help/domains/chem/r-group-analysis

gr.loadDataFrame('/demo/sar_small.csv')
    .then(t => chem.rGroup(t, 'smiles', 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1')
        .then(t => gr.addTableView(t)));

