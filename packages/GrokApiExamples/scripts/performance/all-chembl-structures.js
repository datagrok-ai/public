// Retrieve all molecules from chembl (~2 millions), and visualize in the grid

grok.query('AllChemblStructures', {}).then(t => grok.addTableView(t));
