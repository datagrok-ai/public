// An example of getting R-Groups
// https://datagrok.ai/help/domains/chem/r-group-analysis

grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv')
  .then(t => grok.chem.rGroup(t, 'smiles', 'O=C1Nc2ccccc2C(C2CCCCC2)=NC1')
    .then(t => grok.shell.addTableView(t)));