// dynamically create properties and edit objects in the property grid

const person = {name: 'John', sex: 'M', weight: 67.2, age: 32, checked: true};

// generate properties
const name = DG.Property.string('name', (x) => x.name, (x, v) => x.name = v, '');
const sex = DG.Property.string('sex', (x) => x.sex, (x, v) => x.sex = v, '');

// shorter notation for JS properties
const weight = DG.Property.jsFloat('weight');
const age = DG.Property.jsInt('age');
const checked = DG.Property.jsBool('checked');

// change properties if you wish
sex.choices = ['M', 'F'];

const propGrid = new DG.PropertyGrid();
propGrid.update(person, [name, sex, weight, age, checked]);

//ui.dialog().add(propGrid).show({width: 400, height: 250});

ui.dialog()
  .add(ui.h1('Property grid'))
  .add(propGrid)
  .add(ui.h1('Form'))
  .add(ui.input.form(person, [name, sex, weight, age, checked]))
  .show({width: 350, height: 500});
