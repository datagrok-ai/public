// Different input types

let v = grok.shell.newView('demo: inputs');

let search = ui.input.search('Search');
let name = ui.input.string('Name', {value: 'Arthur Dent'});
let age = ui.input.int('Age', {value: 30});
let date = ui.input.date('Birthday', {value: dayjs('1970-5-10')});
let alien = ui.input.bool('Alien', {value: false});
let friends = ui.input.multiChoice('Friends', {items: ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast'],
  value: ['Ford', 'Fenchurch']});
let bestFriend = ui.radioInput('Best friend', 'Ford', ['Ford', 'Fenchurch']);
let active = ui.input.toggle('Active', {value: true});
let color = ui.colorInput('Favorite color', '#ff0000');
let country = ui.typeAhead('Country', {source: {
  local: ['USA', 'Ukraine', 'Antigua', 'United Kingdom', 'United Arab Emirates']}});
let data = ui.tableInput('Data', null);
let tags = ui.input.tags('Skills', {tags: ['HTML', 'CSS', 'JS'], showBtn: true});
let memo = ui.textInput('Memo', '');

name.addValidator(s => s.length < 15 ? 'Too short' : null);

let form = ui.form([search, data, name, age, date, alien, friends, bestFriend, active, color, country, tags, memo]);

form.append(ui.buttonsInput([ui.bigButton('Submit', ()=>{grok.shell.info('Submit');})]));
v.append(form);