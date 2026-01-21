// Create inputs dynamically for properties

const properties = [
  {
    'name': 'name',
    'type': 'string',
  },
  {
    'name': 'structure',
    'type': 'string',
    'semType': 'Molecule'
  }];

let props = properties.map((p) => DG.Property.fromOptions(p));

let objects = [{
  name: 'aspirin',
  structure: 'CC(=O)OC1=CC=CC=C1C(=O)O',
},
{
  name: 'coffee',
  structure: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
}];

ui.tableFromProperties(objects, props);