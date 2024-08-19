// Retrieve all molecules from chembl (~2 millions), and visualize in the grid

grok.data.query('Chembl:CbAllChemblStructures', {}).then(t => grok.shell.addTableView(t));
