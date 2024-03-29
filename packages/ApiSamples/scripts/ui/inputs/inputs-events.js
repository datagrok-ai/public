let v = grok.shell.newView('demo: inputs-events');

const msg = (v) => grok.shell.info(v);

let name = ui.stringInput('Name', 'Arthur Dent', v => msg(v));
let age = ui.input.int('Age', {value: 30, onValueChanged: v => msg(v)});
let sex = ui.input.choice('Sex', {items: ['Male', 'Female'], value: 'Male', onValueChanged: v => msg(v)});
let date = ui.dateInput('Birthday', dayjs('1970-5-10'), v => msg(v));
let alien = ui.boolInput('Alien', false, v => msg(v));
let friends = ui.input.multiChoice('Friends', {items: ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast'],
  value: ['Ford', 'Fenchurch'], onValueChanged: v => msg(v)});
let active = ui.switchInput('Active', true, (v)=>msg(v));

let inputs = [name, age, sex, date, alien, friends, active];
let container = ui.div();
v.append(container);
container.appendChild(ui.inputs(inputs));

container.appendChild(ui.bigButton('POST', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));
