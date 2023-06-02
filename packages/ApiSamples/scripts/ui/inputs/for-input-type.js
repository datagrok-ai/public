const df = grok.data.demo.molecules();

const properties = [
  { "name": "int",         "inputType": "Int", "valueValidators": [(w) => w > 10], "friendlyName": "IntInput"},
  { "name": "bigInt",      "inputType": "BigInt",                                  "friendlyName": "BigIntInput"},
  { "name": "float",       "inputType": "Float",                                   "friendlyName": "FloatInput"},
  { "name": "qnum",        "inputType": "QNum",                                    "friendlyName": "QNumInput" },
  { "name": "slider",      "inputType": "Slider", "min": 0, "max": 10,             "friendlyName": "SliderInput"},
  { "name": "bool",        "inputType": "Bool",                                    "friendlyName": "BoolInput"},
  { "name": "textArea",    "inputType": "TextArea",                                "friendlyName": "TextAreaInput"},
  { "name": "text",        "inputType": "Text",                                    "friendlyName": "TextInput"},
  { "name": "date",        "inputType": "Date",                                    "friendlyName": "DateInput"},
  { "name": "map",         "inputType": "Map",                                     "friendlyName": "MapInput"},
  { "name": "list",        "inputType": "List",                                    "friendlyName": "ListInput"},
  { "name": "color",       "inputType": "Color",                                   "friendlyName": "ColorInput"},
  { "name": "column",      "inputType": "Column",                                  "friendlyName": "ColumnInput"},
  { "name": "radio",       "inputType": "Radio",                                   "friendlyName": "RadioInput"},
  { "name": "choice",      "inputType": "Choice",                                  "friendlyName": "ChoiceInput"},
  { "name": "multiChoice", "inputType": "MultiChoice",                             "friendlyName": "MultiChoiceInput"},
  { "name": "table",       "inputType": "Table",                                   "friendlyName": "TableInput"},
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
