
//name: #{NAME}
//description: Creates an info panel
//tags: panel
//input: string smiles {semType: Molecule}
//output: widget result
//condition: true
export function #{NAME}(smiles) {
    let panels = ui.div([
        this.createSearchPanel('Exact', smiles),
        this.createSearchPanel('Similar', smiles),
        this.createSearchPanel('Substructure', smiles)
    ]);
    return DG.Widget.fromRoot(panels);
}
