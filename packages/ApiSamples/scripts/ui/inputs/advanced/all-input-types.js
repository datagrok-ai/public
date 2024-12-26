const df = grok.data.demo.demog();
grok.shell.addTable(df);

const properties = [
  {'name': 'int',         'inputType': 'Int', 'showPlusMinus': true},
  {'name': 'bigInt',      'inputType': 'BigInt'},
  {'name': 'float',       'inputType': 'Float', min: 0, max: 10, 'showSlider': true},
  {'name': 'qnum',        'inputType': 'QNum'},
  {'name': 'slider',      'inputType': 'Slider', 'min': 0, 'max': 10},
  {'name': 'bool',        'inputType': 'Bool'},
  {'name': 'switch',      'inputType': 'Switch'},
  {'name': 'textArea',    'inputType': 'TextArea'},
  {'name': 'text',        'inputType': 'Text'},
  {'name': 'search',      'inputType': 'Search'},
  {'name': 'date',        'inputType': 'Date'},
  //{'name': 'map',         'inputType': 'Map'},
  {'name': 'list',        'inputType': 'List'},
  {'name': 'color',       'inputType': 'Color'},
  {'name': 'table',       'inputType': 'Table'},
  {'name': 'column',      'inputType': 'Column'},
  {'name': 'radio',       'inputType': 'Radio', choices: ['Apple', 'Banana']},
  {'name': 'choice',      'inputType': 'Choice', choices: ['Apple', 'Banana']},
  {'name': 'multiChoice', 'inputType': 'MultiChoice', choices: ['Apple', 'Banana']},
  //{'name': 'file',        'inputType': 'File'},
  {'name': 'user',       'inputType': 'User'},
  {'name': 'userGroups',  'inputType': 'UserGroups'},
];

let props = properties.map((p) => DG.Property.fromOptions(p))
let object = {
  int: 12,
  //bigInt: '12345678901234567890',
  float: 4.5,
  qnum: DG.Qnum.less(5),
  slider: 10,
  bool: true,
  switch: true,
  textArea: 'Bins rise and they fall,\nDistribution\'s gentle dance,\nHistogram unveils',
  text: 'Only text',
  search: '',
  //date: '2023/10/30', //dayjs(1977, 10, 2),
  //map: { key1: 'value1', key2: 'value2'},
  list: ['1', '2', '4'],
  color: '#e66465',
  table: df,
  column: df.columns.byName('race'),
  radio: 'Apple',
  choice: 'Apple',
  multiChoice: ['Apple'],
  //file: '',
  user: [DG.User.admin],
  userGroups: [DG.Group.developers],
};

const div = ui.divV([
  ui.input.form(object, props),
  ui.button('Show', () => grok.shell.info(JSON.stringify(object)))
]);

grok.shell.newView('Inputs', [div]);