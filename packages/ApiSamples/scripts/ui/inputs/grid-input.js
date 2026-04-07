const props = [
  //DG.Property.js("fragment", "string", {semType: "Molecule"}),
  DG.Property.js('chemist', 'string'),
  DG.Property.js('stage', 'string', {choices: ['Screening', 'Optimization']}),
  DG.Property.js('color', 'string', {inputType: 'Color'})
];

let molColors = [
  {chemist: 'Jon Snow', stage: 'Screening', fragment: 'CC(=O)OC1=CC=CC=C1C(=O)O', color: '#ff0000'},
  {chemist: 'John Wick', stage: 'Optimization', fragment: 'OC1=CC', color: '#00ff00'}
];

const grid = ui.input.grid(molColors, props);
grid.root;