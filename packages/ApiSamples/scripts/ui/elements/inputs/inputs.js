let v = grok.shell.newView('demo: inputs');

let name = ui.stringInput('Name', 'Arthur Dent');
let age = ui.intInput('Age', 30);
let sex = ui.choiceInput('Sex', 'Male', ['Male', 'Female']);
let date = ui.dateInput('Birthday', DG.DateTime.fromDate(new Date(1970, 5, 10)));
let alien = ui.boolInput('Alien', false);
let friends = ui.multiChoiceInput('Friends', ['Ford', 'Fenchurch'], ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast']);
let active = ui.switchInput('Active', true);
let color = ui.colorInput('Favorite color', '#ff0000');

let inputs = [name, age, sex, date, alien, friends, active];
let container = ui.div();
v.append(container);
container.appendChild(ui.inputs(inputs));

container.appendChild(ui.bigButton('POST', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));
