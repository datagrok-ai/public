// Deprecated. Use UserSettingsStorage.
// User data storage demo: Saving values to storage

const STORAGE_NAME = 'user-data-storage-demo';

let v = grok.shell.newView('demo: user data storage');

let age = ui.input.int('Age', {value: 30});
let sex = ui.input.choice('Sex', {items: ['Male', 'Female'], value: 'Male'});
let music = ui.input.multiChoice('Music genres', {items: ['Classic', 'Rock', 'Pop', 'Jazz']});

let inputs = [age, sex, music];
v.append(ui.inputs(inputs));

let storageButton = ui.iconFA('database', () => {
  grok.dapi.userDataStorage.get(STORAGE_NAME).then((entities) => {
    if (entities !== null && Object.keys(entities).length === 0)
      grok.shell.info('Storage is empty. Try to post something to the storage');
    else {
      let menu = DG.Menu.popup();
      for (let time of Object.keys(entities)) {
        let values = JSON.parse(entities[time]);
        menu.item(values.map(v => `${v.caption}: ${v.value}`).join(', '), () => {
          for (let input of values)
            inputs.find(i => i.caption === input.caption).load(input.value);
        });
      }
      menu.show();
    }
  });
});
storageButton.style.margin = '8px 24px 0 24px';

let postButton = ui.button('Post to storage', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
  grok.dapi.userDataStorage.postValue(STORAGE_NAME, `${Math.floor(Math.random() * Math.floor(1000))}`,
    JSON.stringify(inputs.map(i => {
      return {caption: i.caption, value: i.save()};
    })));
});

let clearButton = ui.button('Clear storage', () => {
  grok.dapi.userDataStorage.remove(STORAGE_NAME, null);
});

v.append(ui.divH([
  storageButton,
  postButton,
  clearButton
]));
