let v = grok.newView('demo: inputs');

let name = ui.stringInput('Name', 'Arthur Dent');
let age = ui.intInput('Age', 30);
let sex = ui.choiceInput('Sex', 'Male', ['Male', 'Female']);
let date = ui.dateInput('Birthday', DateTime.fromDate(new Date(1970, 5, 10)));
let alien = ui.boolInput('Alien', false);

let inputs = [name, age, sex, date, alien];
v.append(ui.inputs(inputs));

v.append(ui.bigButton('POST', () => {
    grok.balloon.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));
