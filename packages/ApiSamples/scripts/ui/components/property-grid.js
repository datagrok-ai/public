// dynamically create properties and edit objects in the property grid

let person = {name: 'John', sex: 'M', weight: 67.2, age: 32, checked: true }

// generate properties
let name = DG.Property.string('name', x => x.name, (x, v) => x.name = v, '');
let sex = DG.Property.string('sex', x => x.sex, (x, v) => x.sex = v, '');

// shorter notation for JS properties
let width = DG.Property.jsFloat('weight');
let age = DG.Property.jsInt('age');
let checked = DG.Property.jsBool('checked');

// change properties if you wish
sex.choices = ['M', 'F'];

let propGrid = new DG.PropertyGrid();
propGrid.update(person, [name, sex, width, age, checked]);

ui.dialog().add(propGrid).show({width: 400, height: 250});