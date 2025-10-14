// Different input types

let v = grok.shell.newView('demo: inputs');

let search = ui.input.search('Search');
let data = ui.input.table('Data');
let name = ui.input.string('Name', {value: 'Arthur Dent'});
let age = ui.input.int('Age', {value: 30});
let date = ui.input.date('Birthday', {value: dayjs('1970-5-10')});
let alien = ui.input.bool('Alien', {value: false});
let friends = ui.input.multiChoice('Friends', {items: ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast'],
  value: ['Ford', 'Fenchurch']});
let bestFriend = ui.input.radio('Best friend', {items: ['Ford', 'Fenchurch'], value: 'Ford'});
let active = ui.input.toggle('Active', {value: true});
let color = ui.input.color('Favorite color', {value: '#ff0000'});
let country = ui.typeAhead('Country', {source: {
  local: ['USA', 'Ukraine', 'Antigua', 'United Kingdom', 'United Arab Emirates']}});
let tags = ui.input.tags('Skills', {tags: ['HTML', 'CSS', 'JS'], showButton: true});
let memo = ui.input.textArea('Memo');

name.addValidator(s => s.length < 15 ? 'Too short' : null);

let form = ui.form([search, data, name, age, date, alien, friends, bestFriend, active, color, country, tags, memo]);

form.append(ui.buttonsInput([ui.bigButton('Submit', () => { grok.shell.info('Submit'); })]));
v.append(form);