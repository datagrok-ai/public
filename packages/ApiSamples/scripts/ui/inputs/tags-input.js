let v = grok.shell.newView('Demo');

//standard tags input
const tagsInput = ui.input.tags('Tags', {tags: ['a', 'abb', 'aacc'], value: ['a']});

//user input with predefined list of items
const u1 = DG.User.create();
u1.firstName = 'test_name';
const u2 = DG.User.create();
u2.firstName = 'test_name_2';

const userInput = ui.input.user('User', {items: [u1, u2]});

const customTagsInput = ui.input.tags('Custom tags input', {
  itemToString: (item) => {return `${item.a}_${item.b}`;},
  getSuggestions: (inputText) => {return [{a: 'first', b: 'object'}, {a: 'second', b: 'object'}];},
  createTagLabel: (item) => {return ui.divH([ui.divText(`${item.a}_${item.b}`)]);},
  createNewItem: (text) => {return {a: 'created', b: 'object'};},
  allowNew: true,
});

v.append(ui.divV([
  tagsInput.root,
  userInput.root,
  customTagsInput.root,
  ui.button('show value', () => {
    console.log(tagsInput.value);
    console.log(userInput.value);
    console.log(customTagsInput.value);

  })]));