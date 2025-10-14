// Retrieve all molecules from chembl (~2 millions), and visualize in the grid

let t = await grok.data.query('Chembl:CbAllChemblStructures', {});
grok.shell.addTableView(t);
