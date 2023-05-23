const df = grok.data.demo.molecules();

const properties = [
  { "name": "int",         "inputType": "Int", "valueValidators": [(w) => w > 10] },
  { "name": "bigInt",      "inputType": "BigInt" },
  { "name": "float",       "inputType": "Float" },
  { "name": "qnum",        "inputType": "QNum" },
  { "name": "slider",      "inputType": "Slider", "min": 0, "max": 10 },
  { "name": "bool",        "inputType": "Bool" },
  { "name": "textArea",    "inputType": "TextArea" },
  { "name": "text",        "inputType": "Text"},
  { "name": "date",        "inputType": "Date"},
  { "name": "map",         "inputType": "Map"},
  { "name": "list",        "inputType": "List"},
  { "name": "color",       "inputType": "Color"},
  { "name": "column",      "inputType": "Column"},
  { "name": "radio",       "inputType": "Radio"},
  { "name": "choice",      "inputType": "Choice"},
  { "name": "multiChoice", "inputType": "MultiChoice"},
  { "name": "table",       "inputType": "Table"},
  ];
let props = properties.map((p) => DG.Property.fromOptions(p))
let object = {
  int: 12,
  float: 4.5,
  slider: 10,
  bool: true,
  textArea: 'What a beautiful day',
  text: 'Only text',
  list: ['1', '2', '4'],
  column: df.columns.byName('smiles'),
  choice: ['1'],
  table: df
};
ui.input.form(object, props);