let v = grok.shell.newView('demo: inputs-events');

const msg = (v) => grok.shell.info(v);

let name = ui.stringInput('Name', 'Arthur Dent', v => msg(v));
let age = ui.intInput('Age', 30, v => msg(v));
let sex = ui.choiceInput('Sex', 'Male', ['Male', 'Female'], v => msg(v));
let date = ui.dateInput('Birthday', DG.DateTime.fromDate(new Date(1970, 5, 10)), v => msg(v));
let alien = ui.boolInput('Alien', false, v => msg(v));
let friends = ui.multiChoiceInput('Friends',
  ['Ford', 'Fenchurch'], ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast'], v => msg(v));

let inputs = [name, age, sex, date, alien, friends];
let container = ui.div();
v.append(container);
container.appendChild(ui.inputs(inputs));

container.appendChild(ui.bigButton('POST', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));
