// Modifying choice input list programmatically

// Method #1. Using grok's wrapper classes only
let v = grok.shell.newView('Demo');
let choices = ui.choiceInput('Value', 'A', ['A', 'B', 'C']);
let container = ui.div();
v.append(container);
let inputs = ui.inputs([choices]);
container.appendChild(inputs);
let choices_new = ui.choiceInput('Value', 'B', ['B', 'C', 'D']);
inputs.replaceChild(choices_new.root, choices.root);
choices = choices_new;
choices.onChanged((v) => {
  ui.dialog().show();
});

// Method #2. Using raw DOM manipulation
let items = ['C', 'D', 'E'];
let choicesDOM = choices.input;
for (index = 0; index < items.length; ++index) {
  choicesDOM.options[index] = new Option(items[index], index);
}

// Method #3. Using jQuery syntax
$(choicesDOM).append(new Option('F', 'F'));