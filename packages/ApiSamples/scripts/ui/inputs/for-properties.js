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
    "name": "started",
    "type": DG.TYPE.DATE_TIME,
  },
  {
    "name": "rating",
    "type": "int",
    "editor": "slider",
    "min": 0,
    "max": 10,
  },
  {
    "name": "department",
    "type": "string",
    "choices": ["IT", "Business"]
  }];

let props = properties.map((p) => DG.Property.fromOptions(p))

let object = {};

ui.input.form(object, props)