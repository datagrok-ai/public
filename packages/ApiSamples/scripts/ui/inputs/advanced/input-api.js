const df = grok.data.demo.demog();
grok.shell.addTable(df);

const inputs = [
  ui.input.int('int', {value: 6, nullable: false, min: 2, max: 12, step: 2, showPlusMinus: true}),
  // ui.input.bigInt('bigInt', {value: '12345678901234567890'}), // doesn't work yet with value
  ui.input.float('float', {value: 6.5758, min: 0, max: 10, step: 0.5, format: '#0.000', showSlider: true}),
  ui.input.qNum('qNum', {value: 5.959}),
  ui.input.slider('slider', {value: 6, min: 0, max: 20, step: 2}),
  ui.input.bool('bool', {value: false, tooltipText: 'Bool input tooltip'}),
  ui.input.toggle('toggle', {value: true, tooltipText: 'Toggle input tooltip'}),
  ui.input.textArea('textArea', {value: 'Text area', size: {width: 150, height: 35}}),
  ui.input.string('string', {value: 'Text', placeholder: 'Text', clearIcon: true, escClears: true, icon: ui.iconFA('square')}),
  ui.input.search('search', {value: 'Search'}),
  ui.input.date('date', {value: dayjs('1970-5-10')}),
  // ui.input.map('map'), // doesn't work yet
  ui.input.file('file', {value: (await grok.dapi.files.list('System:AppData/ApiTests/datasets/'))[1]}),
  ui.input.list('list', {value: [1, 2, 3]}),
  ui.input.color('color', {value: '#e66465'}),
  ui.input.column('column', {value: df.col('height'), table: df, filter: (col) => col.type === DG.COLUMN_TYPE.FLOAT}),
  ui.input.columns('columns', {table: df, value: df.columns.byNames(['weight', 'age']), available: ['height', 'weight', 'age']}),
  // ui.input.columnsMap('columnsMap'), // doesn't work yet
  ui.input.radio('radio', {items: ['apple', 'orange'], value: 'apple'}),
  ui.input.choice('choice', {items: ['apple', 'orange'], value: 'apple'}),
  ui.input.multiChoice('multiChoice', {items: ['apple', 'orange', 'pear'], value: ['apple', 'orange']}),
  ui.input.table('table', {value: df}),
  ui.input.molecule('molecule', {value: 'CN1CCC(O)(CC1)c2ccccc2'}),
  ui.input.user('user', {value: [DG.User.current()]}),
  ui.input.userGroups('userGroups', {value: [DG.User.current()]}),
  ui.input.image('image', {value: 'https://datagrok.ai/img/logo.svg'}),
  ui.input.tags('tags', {tags: ['Tag 1', 'Tag2', 'Tag3'], showButton: true}),
  ui.input.code('code', {script: 'print(\'Hello world!\')', mode: 'python', placeholder: 'print(\'Hello world!\')'}),
  await ui.input.markdown('markdown'),
];

const form = DG.InputForm.forInputs(inputs);
form.onInputChanged.subscribe((ed) => grok.shell.info(`Input - ${ed.args.input.caption}, new value - ${ed.args.input.value}`));

grok.shell.newView('Inputs', [form.root]);
