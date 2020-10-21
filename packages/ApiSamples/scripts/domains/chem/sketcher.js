// Sketching a molecule

let molfileInput = ui.textArea();
let smilesInput = ui.textInput();

function onChanged(smiles, molfile) {
    smilesInput.value = smiles;
    molfileInput.value = molfile;
}

grok.shell.newView('sketcher', [
    grok.chem.sketcher(onChanged, 'CC(=O)Oc1ccccc1C(=O)O)'),
    smilesInput,
    molfileInput,
]);