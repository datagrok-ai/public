const df = grok.data.demo.demog();
grok.shell.addTable(df);

const properties = [
  { "name": "int",         "inputType": "Int", "showPlusMinus": true },
  { "name": "float",       "inputType": "Float", min: 0, max: 10, "showSlider": true},
  //{ "name": "bigInt",      "inputType": "BigInt"},
  { "name": "qnum",        "inputType": "QNum"},
  { "name": "slider",      "inputType": "Slider", "min": 0, "max": 10},
  { "name": "bool",        "inputType": "Bool"},
  { "name": "textArea",    "inputType": "TextArea"},
  { "name": "text",        "inputType": "Text"},
  //{ "name": "date",        "inputType": "Date"},
  //{ "name": "map",         "inputType": "Map"},
  { "name": "list",        "inputType": "List"},
  { "name": "color",       "inputType": "Color"},
  { "name": "table",       "inputType": "Table"},
  { "name": "column",      "inputType": "Column"},
  { "name": "radio",       "inputType": "Radio", choices: ["Apple", "Banana"]},
  { "name": "choice",      "inputType": "Choice", choices: ["Apple", "Banana"]},
  { "name": "multiChoice", "inputType": "MultiChoice", choices: ["Apple", "Banana"]},
  { "name": "users",       "inputType": "UserGroups"},
];

let props = properties.map((p) => DG.Property.fromOptions(p))
let object = {
  int: 12,
  //bigInt: "12345678901234567890",
  float: 4.5,
  qnum: DG.Qnum.less(5),
  slider: 10,
  bool: true,
  textArea: "Bins rise and they fall,\nDistribution's gentle dance,\nHistogram unveils",
  text: 'Only text',
  //date: '2023/10/30', //dayjs(1977, 10, 2),
  //map: { key1: 'value1', key2: 'value2'},
  color: '#e66465',
  list: ['1', '2', '4'],
  table: df,
  column: df.columns.byName('race'),
  choice: 'Apple',
  multiChoice: ['Apple'],
  radio: 'Apple',
  users: [DG.User.current()]
};

const div = ui.divV([
  ui.input.form(object, props),
  ui.button('Show', () => grok.shell.info(JSON.stringify(object)))
])

grok.shell.newView('Inputs', [div]);