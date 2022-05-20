// adding the specified filter to the already existing filter group

let tv = grok.shell.addTableView(grok.data.demo.molecules());
let options = {
  type: "Chem:substructureFilter",
  column: "smiles",
  columnName: "smiles"
};
tv.filtersGroup.add(options);