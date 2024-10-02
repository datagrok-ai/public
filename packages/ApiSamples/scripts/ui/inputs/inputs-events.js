let v = grok.shell.newView('demo: inputs-events');

const msg = (v) => grok.shell.info(v);

let name = ui.input.string('Name', {value: 'Arthur Dent', onValueChanged: (v) => msg(v)});
let age = ui.input.int('Age', {value: 30, onValueChanged: (v) => msg(v)});
let sex = ui.input.choice('Sex', {items: ['Male', 'Female'], value: 'Male', onValueChanged: (v) => msg(v)});
let date = ui.input.date('Birthday', {value: dayjs('1970-5-10'), onValueChanged: (v) => msg(v)});
let alien = ui.input.bool('Alien', {value: false, onValueChanged: (v) => msg(v)});
let friends = ui.input.multiChoice('Friends', {items: ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast'],
  value: ['Ford', 'Fenchurch'], onValueChanged: (v) => msg(v)});
let active = ui.input.toggle('Active', {value: true, onValueChanged: (v) => msg(v)});

let inputs = [name, age, sex, date, alien, friends, active];
let container = ui.div();
v.append(container);
container.appendChild(ui.inputs(inputs));

container.appendChild(ui.bigButton('POST', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));
