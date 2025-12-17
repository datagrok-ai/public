
//name: #{NAME}
//description: Creates an info panel
//meta.role: panel
//input: string smiles {semType: Molecule}
//output: widget result
//condition: true
export function #{NAME}(smiles: string) {
  let mol = ui.div(grok.chem.svgMol(smiles));
  return DG.Widget.fromRoot(mol);
}
