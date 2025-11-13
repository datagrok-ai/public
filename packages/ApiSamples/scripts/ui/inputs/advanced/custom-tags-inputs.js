let v = grok.shell.newView('Demo');

//standard tags input
const tsgsInput = ui.input.tags('Tags', {tags: ['a', 'abb', 'aacc'], value: ['a']});

//user input with predefined list of items
const u1 = DG.User.create();
u1.firstName = 'test_name';
const u2 = DG.User.create();
u2.firstName = 'test_name_2';

const userInput = ui.input.user('User', {items: [u1, u2]});

/* create custom tags input for any js object by overriding TagsInput methods
For instance, example object is of type
  type MyObj = {
    a: string,
    b: string,
  } */


class TagsInputCustom extends DG.TagsInput {

  constructor(dart) {
    super(dart);
  };

  itemToString(item) {return `${item.a}_${item.b}`};

  
  async getSuggestions(inputText) { 
    return [{a: 'first', b: 'object'}, {a: 'second', b: 'object'}];
  };

  createTag(item) {
    return ui.divH([ui.icons.delete(() => {
      const idxToRemove = this.selectedItems.findIndex((it) => it.a === item.a && it.b === item.b);
      if (idxToRemove !== -1)
        this.removeItemByIdx(idxToRemove);
    }), ui.divText(`${item.a}_${item.b}`)]);
  };

  async gatherItems(selectedItems, selectedItem, findItems) {
    if (selectedItem)
      selectedItems.push(selectedItem);
    return selectedItems;
  };
}

const customTagsInput = ui.input.tags('Custom tags input', {createCustomInputFunc: (dart) => new TagsInputCustom(dart)});

v.append(ui.divV([
  tsgsInput.root,
  userInput.root,
  customTagsInput.root,
  ui.button('show value', () => {
    console.log(tsgsInput.value);
    console.log(userInput.value);
    console.log(customTagsInput.value);

  })]));