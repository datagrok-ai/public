// Different input types

// Adding an icon to each input that lets you inspect the value
// grok.events.onInputCreated.subscribe(
//  (i) => i.root.appendChild(ui.iconFA('eye', () => grok.shell.o = i.value)));

// Different input types

// Adding an icon to each input that lets you inspect the value
// grok.events.onInputCreated.subscribe(
//  (i) => i.root.appendChild(ui.iconFA('eye', () => grok.shell.o = i.value)));

let v = grok.shell.newView('demo: inputs');

let type = ui.choiceInput('Form type', 'Default', ['Default', 'Left', 'Wide', 'Condensed'], (x)=>{
  switch (x){
    case 'Default': 
      form.className = 'ui-form';
      break;  
    case 'Left': 
      form.className = 'ui-form ui-form-left';
      break;
    case 'Wide': 
      form.className = 'ui-form ui-form-wide';
      break;
    case 'Condensed': 
      form.className = 'ui-form ui-form-condensed';
      break;  
    default: form.className = '.ui-form';
  }
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

form.append(ui.buttonsInput([ui.bigButton('Submit', ()=>{grok.shell.info(form.className)})]));
v.append(form);
