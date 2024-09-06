---
title: "User settings storage"
---

It is often the case that certain settings need to be shared between different applications or
different instances of the same application. This functionality is implemented in the form of User Settings Storage 
that allows storage of key-value pairs.

Table of contents:

* [General Structure](#general-structure)
* [JavaScript API](#javascript-api)
    * [Saving](#saving)
    * [Loading](#loading)
    * [Exploring](#exploring)
    * [Deleting](#deleting)
* [Examples and Use Cases](#examples-and-use-cases)

## General structure

The main class `grok.userSettings` has the following 6 methods:

| Method        | Parameters                                                                                                           | Function                                                  |
|---------------|----------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| `.add()`      | <b>name</b>: <i>string</i>, <b>key</b>: <i>string</i>, <b>value</b>: <i>string</i>, <b>isPrivate</b>: <i>boolean</i> | Saves a single value to storage                           |
| `.addAll()`   | <b>name</b>: <i>string</i>, <b>data</b>: <i>Map</i>, <b>isPrivate</b>: <i>boolean</i>                                | Saves a map to storage, will be appended to existing data |
| `.put()`      | <b>name</b>: <i>string</i>, <b>data</b>: <i>Map</i>, <b>isPrivate</b>: <i>boolean</i>                                | Saves a map to storage, will replace existing data        |
| `.get()`      | <b>name</b>: <i>string</i>, <b>isPrivate</b>: <i>boolean</i>                                                         | Retrieves a map from storage                              |
| `.getValue()` | <b>name</b>: <i>string</i>, <b>key</b>: <i>string</i>, <b>isPrivate</b>: <i>boolean</i>                              | Retrieves a single value from storage                     |
| `.delete()`   | <b>name</b>: <i>string</i>, <b>key</b>: <i>string</i>, <b>isPrivate</b>: <i>boolean</i>                              | Removes a single value from storage                       |

## JavaScript API

### Saving

When using `.add()` specifically, the input can only be represented by a *string*. The best way to organize
entries in user storage is by creating *objects* which can then be converted to *
string* format with `JSON.stringify()`. For example:

```js
let inputObj = {x: 1, y: 2, z: 3};
function store(STORAGE_NAME,key,value) {
  grok.userSettings.add(STORAGE_NAME, key, JSON.stringify(value));
}

store('coordinate-storage', 'coordinates', inputObj)
```

:::note

Maximum allowed value length is 5000 symbols.

:::

### Loading

To return values, we execute the following code:

```js
let entry = grok.userSettings.getValue('coordinate-storage','coordinates');
let outputObj = JSON.parse(entry);
```

### Exploring

Unfortunately we can't explore the content of user storage without first retrieving it from the memory. Once retrieved,
we can use a standard `for` loop or a `forEach()` method to scroll through the entries:

 ```js
let entries = grok.userSettings.get('coordinate-storage');
Object.keys(entries).forEach((key) => {
  let entry = JSON.parse(entries[key]);
  //from here we can view the content of each record
});
```

### Deleting

We can clear a specific value from storage by executing:

```js
grok.userSettings.delete('coordinate-storage', 'coordinates');
```

or wipe the storage completely by running:

```js
grok.userSettings.delete('coordinate-storage', null);
```

## Examples and use cases

User storage is especially handy when we want to save sets of input parameters, e.g. if a user wishes to reproduce their
previous calculations in following sessions. All of the above snippets are summarized in the example below:

```js
const STORAGE_NAME = 'user-data-storage-demo';
let v = grok.shell.newView('demo: user data storage');
let age = ui.input.int('Age', {value: 30});
let sex = ui.input.choice('Sex', {value: 'Male', items: ['Male', 'Female']});
let music = ui.input.multiChoice('Music genres', {value: null, items: ['Classic', 'Rock', 'Pop', 'Jazz']});
let inputs = [age, sex, music];
v.append(ui.inputs(inputs));

let storageButton = ui.iconFA('database', () => {
  let entities = grok.userSettings.get(STORAGE_NAME);
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
storageButton.style.margin = '8px 24px 0 24px';

let postButton = ui.button('Post to storage', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
  grok.userSettings.add(STORAGE_NAME, `${Math.floor(Math.random() * Math.floor(1000))}`,
      JSON.stringify(inputs.map(i => {
        return {caption: i.caption, value: i.save()};
      })));
});

let clearButton = ui.button('Clear storage', () => {
  grok.userSettings.delete(STORAGE_NAME, null);
});

v.append(ui.divH([
  storageButton,
  postButton,
  clearButton
]));
```
