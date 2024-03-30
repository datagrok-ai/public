const df = grok.data.demo.demog();
grok.shell.addTable(df);

const inputs = [
  ui.input.int('int'),
  // ui.input.bigInt('bigInt'), // doesn't work yet
  ui.input.float('float'),
  ui.input.qNum('qNum'),
  ui.input.slider('slider'),
  ui.input.bool('bool'),
  ui.input.toggle('toggle'),
  ui.input.textArea('textArea'),
  ui.input.string('string'),
  ui.input.search('search'),
  ui.input.date('date'),
  // ui.input.map('map'), // doesn't work yet
  ui.input.file('file'),
  // ui.input.list('list'), // doesn't work yet
  ui.input.color('color'),
  // ui.input.column('column'), // doesn't work yet
  // ui.input.columns('columns'), // doesn't work yet
  // ui.input.columnsMap('columnsMap'), // doesn't work yet
  ui.input.radio('radio', {items: ['apple', 'orange']}),
  ui.input.choice('choice', {items: ['apple', 'orange']}),
  ui.input.multiChoice('multiChoice', {items: ['apple', 'orange', 'pear']}),
  ui.input.table('table'),
  // ui.input.molecule('molecule'), // doesn't work yet
  // ui.input.users('users'), // doesn't work yet
  ui.input.userGroups('userGroups'),
  // ui.input.dynamic('dynamic'), // doesn't work yet
  // ui.input.image('image'), // doesn't work yet
  // ui.input.jsInput('jsInput'), // doesn't work yet,
  ui.input.tags('tags', {tags: ['Tag 1', 'Tag2', 'Tag3'], showButton: true}),
  ui.input.code('code', {script: 'print(\'Hello world!\')', mode: 'python', placeholder: 'print(\'Hello world!\')'}),
  await ui.input.markdown('markdown'),
];

grok.shell.newView('Inputs', [ui.form(inputs)]);