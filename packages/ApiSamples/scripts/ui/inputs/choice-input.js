// Modifying choice input list programmatically

// Method #1. Using grok's wrapper classes only
let v = grok.shell.newView('Demo');
let choices = ui.input.choice('Value', {items: ['A', 'B', 'C'], value: 'A'});
let container = ui.div();
v.append(container);
let inputs = ui.inputs([choices]);
container.appendChild(inputs);
let choicesNew = ui.input.choice('Value', {items: ['B', 'C', 'D'], value: 'B'});
inputs.replaceChild(choicesNew.root, choices.root);
choices = choicesNew;
choices.onChanged.subscribe((v) => {
  grok.shell.info('The selected value is changed');
});
choices.value = 'C';

// Method #2. Using raw DOM manipulation
let items = ['C', 'D', 'E'];
let choicesDOM = choices.input;
for (index = 0; index < items.length; ++index) {
  choicesDOM.options[index] = new Option(
    items[index],
    items[index]); // this is a name you select by
}
// Selecting an item with a name 'D'
// The event choices.onChanged WON'T be triggered
choicesDOM.value = 'D';

// Method #3. Using jQuery syntax
$(choicesDOM).append(new Option('F', 'F'));