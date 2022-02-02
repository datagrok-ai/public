<!-- TITLE: User storage -->
<!-- SUBTITLE: -->

# User storage

It is often the case that certain settings or inputs/outputs need to be shared between different applications or
different instances of the same application. This functionality is implemented in the form of user data storage — a
virtual memory buffer which can be filled with new entries and from which these entries can later be retrieved.

Table of contents:

* [General Structure](#general-structure)
* [JavaScript API](#javascript-api)
    * [Saving](#saving)
    * [Loading](#loading)
    * [Exploring](#exploring)
    * [Deleting](#deleting)
* [Examples and Use Cases](#examples-and-use-cases)

## General structure

The main class `grok.dapi.userDataStorage` is extended by the following 6 methods:

| Method         | Parameters                                                                                                             | Function                                                            |
|----------------|------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| `.postValue()` | <b>name</b>: <i>string</i>, <b>key</b>: <i>string</i>, <b>value</b>: <i>string</i>, <b>currentUser</b>: <i>boolean</i> | Saves a single value to User Data Storage                           |
| `.post()`      | <b>name</b>: <i>string</i>, <b>data</b>: <i>Map</i>, <b>currentUser</b>: <i>boolean</i>                                | Saves a map to User Data Storage, will be appended to existing data |
| `.put()`       | <b>name</b>: <i>string</i>, <b>data</b>: <i>Map</i>, <b>currentUser</b>: <i>boolean</i>                                | Saves a map to User Data Storage, will replace existing data        |
| `.get()`       | <b>name</b>: <i>string</i>, <b>currentUser</b>: <i>boolean</i>                                                         | Retrieves a map from User Data Storage                              |
| `.getValue()`  | <b>name</b>: <i>string</i>, <b>key</b>: <i>string</i>, <b>currentUser</b>: <i>boolean</i>                              | Retrieves a single value from User Data Storage                     |
| `.remove()`    | <b>name</b>: <i>string</i>, <b>key</b>: <i>string</i>, <b>currentUser</b>: <i>boolean</i>                              | Removes a single value from User Data Storage                       |

Each of the above methods returns a **Promise** which has to be handled in one of the following ways:

1. By `.then()` method

    ```js
    grok.dapi.userDataStorage.method().then(() => {
        //continue your code
    })
    ```

2. By an `await` operator

    ```js
    async function fun() {
        let val = await grok.dapi.userDataStorage.method();
        return val;
    }
    ```

## JavaScript API

### Saving

When using `.postValue()` specifically, the input can only be represented by a *string*. The best way to organize
entries in user storage is by creating *objects* which can then be converted to *
string* format with `JSON.stringify()`. For example:

```js
let inputObj = {x: 1, y: 2, z: 3};
async function store(STORAGE_NAME,key,value) {
    await grok.dapi.userDataStorage.postValue(STORAGE_NAME, key, JSON.stringify(value));
}

store('coordinate-storage', 'coordinates', inputObj)
```

Please note that `JSON.stringify()` will only operate on simple *objects*, i.e. objects with `key:value` pairs where
`value = class DG.DataFrame` will return `undefined`. We suggest
using [`grok.dapi.tables.uploadDataFrame()`](https://dev.datagrok.ai/js/samples/data-access/save-and-load-df)
for a table storage.

### Loading

To return values, we execute the following code:

```js
grok.dapi.userDataStorage.getValue('coordinate-storage','coordinates').then((entry) => {
    let outputObj = JSON.parse(entry);
    //at this point outputObj === inputObj from the previous example
    //all follow up operations can be performed inside current then() statement
});
```

### Exploring

Unfortunately we can't explore the content of user storage without first retrieving it from the memory. Once retrieved,
we can use a standard `for` loop or a `forEach()` method to scroll through the entries:

 ```js
grok.dapi.userDataStorage.get('coordinate-storage').then((entries) => {
    Object.keys(entries).forEach((key) => {
        let entry = JSON.parse(entries[key]);
        //from here we can view the content of each record
    });
});
```

### Deleting

We can clear a specific value from storage by executing:

```js
await grok.dapi.userDataStorage.remove('coordinate-storage', 'coordinates');
```

or wipe the storage completely by running:

```js
await grok.dapi.userDataStorage.remove('coordinate-storage', null);
```

## Examples and use cases

User storage is especially handy when we want to save sets of input parameters, e.g. if a user wishes to reproduce their
previous calculations in following sessions. All of the above snippets are summarized in the example below:

```js
const STORAGE_NAME = 'user-data-storage-demo';
let v = grok.shell.newView('demo: user data storage');
let age = ui.intInput('Age', 30);
let sex = ui.choiceInput('Sex', 'Male', ['Male', 'Female']);
let music = ui.multiChoiceInput('Music genres', null, ['Classic', 'Rock', 'Pop', 'Jazz']);
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
```

Please note: here the `await` operator is not used with `remove()` and `postValue()` because both of these are attached
to buttons, i.e the order of execution is controlled manually.

See also:

* [JavaScript API Samples: User data storage](https://public.datagrok.ai/js/samples/misc/user-data-storage)
* [JavaScript API Samples: Dataframe upload](https://dev.datagrok.ai/js/samples/data-access/save-and-load-df)
