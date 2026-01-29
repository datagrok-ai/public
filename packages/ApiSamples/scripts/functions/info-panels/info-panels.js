/**
 * Register a function that provides an information panel for a molecule.
 * It gets natively integrated with the platform's data augmentation system.
 * See also: https://datagrok.ai/help/discover/info-panels
 */

grok.functions.register({
  signature: 'widget molWidget(String/Molecule smiles)',
  tags: 'panel',
  run: (smiles) => new DG.Widget(ui.divText(`mol panel: ${smiles}`))
});

// To see the panel, click on a cell with a molecule
let df = await grok.data.getDemoTable('chem/smiles.csv');
grok.shell.addTableView(df);
