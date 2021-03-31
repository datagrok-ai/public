// Sketching a molecule

let molfileInput = ui.textInput();
let smilesInput = ui.StringInput();

function onChanged(smiles, molfile) {
  smilesInput.value = smiles;
  molfileInput.value = molfile;
}

grok.shell.newView('sketcher', [
  grok.chem.sketcher(onChanged, 'CC(=O)Oc1ccccc1C(=O)O)'),
  smilesInput.root,
  molfileInput.root,
]);