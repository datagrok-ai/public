// Different form types

let v = grok.shell.newView('demo: form types');

//Basic form
//To set labels aligned to left add class .ui-form-left
let form = ui.form([
  ui.input.string('First name', {value: 'Arthur'}),
  ui.input.string('Last name', {value: 'Dent'}),
  ui.input.int('Age', {value: 30}),
]);

//Add class .ui-form-wide to make the form wide
let wideForm = ui.wideForm([
  ui.input.string('First name', {value: 'Arthur'}),
  ui.input.string('Last name', {value: 'Dent'}),
  ui.input.int('Age', {value: 30}),
]);

//Condensed form. You can also use .ui-form-condensed
let narrowForm = ui.narrowForm([
  ui.input.string('First name', {value: 'Arthur'}),
  ui.input.string('Last name', {value: 'Dent'}),
  ui.input.int('Age', {value: 30}),
]);

v.append(ui.divV([
  ui.h1('Basic form'),
  form,
  ui.h1('Wide form'),
  wideForm,
  ui.h1('Condensed form'),
  narrowForm,
]));
