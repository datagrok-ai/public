// Register a function that provides an information panel for a molecule.
// It gets natively integrated with the platform's data augmentation system.
// To see the panel in action, open a file containing molecules, and click on a molecule.
//
// See also: // https://datagrok.ai/help/concepts/info-panels

gr.functions.register({
    signature: 'widget molWidget(String/Molecule smiles)',
    tags: 'panel',
    run: function(smiles) {
        return new Widget(ui.divText(`mol panel: ${smiles}`));
    }
});