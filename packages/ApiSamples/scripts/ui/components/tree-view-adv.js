let apple = {name: 'apple', type: 'fruit', color: 'red'};
let cherry = {name: 'cherry', type: 'fruit', color: 'red'};
let banana = {name: 'banana', type: 'fruit', color: 'orange'};
let ferrari = {name: 'ferrari', type: 'car', color: 'red'};
let items = [apple, cherry, banana, ferrari];

const tree = DG.TreeViewGroup.fromItemCategories(items,
  ['color', 'type'],
  {itemToString: (x) => x.name}).root;

grok.shell.newView('Tree', [tree]);