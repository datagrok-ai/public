// Different input types



let v = grok.shell.newView('demo: inputs');

let type = ui.choiceInput('Form type', 'Default', ['Default', 'Left', 'Wide', 'Condensed'], (x)=>{
  if (x == 'Default')
    form.className = 'ui-form';
  if (x == 'Left')
    form.className = 'ui-form ui-form-left';  
  if (x == 'Wide')
    form.className = 'ui-form ui-form-wide'; 
  if (x == 'Condensed')
    form.className = 'ui-form ui-form-condensed';   
});

let name = ui.stringInput('Name', 'Arthur Dent');
let age = ui.intInput('Age', 30);
let date = ui.dateInput('Birthday', dayjs('1970-5-10'));
let alien = ui.boolInput('Alien', false);
let friends = ui.multiChoiceInput('Friends', ['Ford', 'Fenchurch'], ['Ford', 'Fenchurch', 'Zaphod', 'Slartibartfast']);
let active = ui.switchInput('Active', true);
let color = ui.colorInput('Favorite color', '#ff0000');
let data = ui.tableInput('Details', null);
let memo = ui.textInput('Memo', '');

let form = ui.form([type, name, age, date, alien, friends, active, color, data, memo]);

form.append(ui.buttonsInput([ui.bigButton('Submit', ()=>{grok.shell.info('Form type: '+form.className)})]));
v.append(form);
