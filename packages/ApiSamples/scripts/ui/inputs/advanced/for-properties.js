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
    "valueValidators": [(w) => w > 10 ? null : 'Too small']
  },
  {
    "name": "started",
    "type": DG.TYPE.DATE_TIME,
  },
  {
    "name": "stage",
    "type": "int",
    "showPlusMinus": true
  },
  {
    "name": "slider",
    "type": DG.TYPE.FLOAT,
    "showSlider": true,
    "min": 2.2,
    "max": 5.5,
    "nullable": false,
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
  stage: 5,
  slider: 3.3,
  rating: 7,
  department: 'IT'
};

ui.input.form(object, props)