// Different form types

let v = grok.shell.newView('demo: form types');

//Basic form
//To set labels aligned to left add class .ui-form-left
let form = ui.form([
  ui.h1('Basic form'),
  ui.stringInput('First name', 'Arthur'),
  ui.stringInput('Last name', 'Dent'),
  ui.input.int('Age', {value: 30}),
]);

//Add class .ui-form-wide to make the form wide
let wideForm = ui.form([
  ui.h1('Wide form'),
  ui.stringInput('First name', 'Arthur'),
  ui.stringInput('Last name', 'Dent'),
  ui.input.int('Age', {value: 30}),
], 'ui-form-wide');

//Condensed form. You can also use .ui-form-condensed
let narrowForm = ui.narrowForm([
  ui.h1('Condensed form'),
  ui.stringInput('First name', 'Arthur'),
  ui.stringInput('Last name', 'Dent'),
  ui.input.int('Age', {value: 30}),
]);
v.append(ui.panel([narrowForm]));

v.append(ui.divV([
  form,
  wideForm,
  narrowForm,
]));
