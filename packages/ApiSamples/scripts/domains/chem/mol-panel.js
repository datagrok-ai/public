/**
 * Register a function that provides an information panel for a molecule.
 * It gets natively integrated with the platform's data augmentation system.
 * To see the panel in action, open a file containing molecules, and click on a molecule.
 * See also: https://datagrok.ai/help/discover/info-panels
 */

grok.functions.register({
  signature: 'widget molWidget(String/Molecule smiles)',
  tags: 'panel',
  run: (smiles) => {
    return new DG.Widget(ui.divText(`mol panel: ${smiles}`));
  }
});
