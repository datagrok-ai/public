// Sketching a molecule

let molfileInput = ui.input.textArea('');
let smilesInput = ui.input.string('');

function onChanged(smiles, molfile) {
  smilesInput.value = smiles;
  molfileInput.value = molfile;
}

grok.shell.newView('sketcher', [
  grok.chem.sketcher(onChanged, 'CC(=O)Oc1ccccc1C(=O)O)'),
  smilesInput.root,
  molfileInput.root,
]);