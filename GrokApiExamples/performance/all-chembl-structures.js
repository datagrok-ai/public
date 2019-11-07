// Retrieve all molecules from chembl (~2 millions), and visualize in the grid

gr.query('AllChemblStructures', {}).then(t => gr.addTableView(t));
