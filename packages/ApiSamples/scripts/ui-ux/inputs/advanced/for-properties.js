// Create inputs dynamically for properties

const properties = [
  {
    "name": "reactorType",
    "type": "string",
    "choices": ["Experimental", "Production"]
  },
  {
    "name": "structure",
    "type": "string",
    "semType": "Molecule"
  },
  {
    "name": "weight",
    "type": DG.TYPE.FLOAT,
    "valueValidators": [(w) => w > 10]
  },
  {
    "name": "started",
    "type": DG.TYPE.DATE_TIME,
  },
  {
    "name": "rating",
    "type": "int",
    "editor": "slider",
    "min": 0,
    "max": 10,
  }];

let props = properties.map((p) => DG.Property.fromOptions(p))

let object = {
  reactorType: 'Production',
  structure: 'CC(=O)OC1=CC=CC=C1C(=O)O',
  rating: 7,
  department: 'IT'
};

ui.input.form(object, props)