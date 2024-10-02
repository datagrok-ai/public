let name = DG.Property.string('Name', null, null, '');
let age = DG.Property.int('Age', null, null, 30);
let height = DG.Property.float('Height', null, null, 175.5);
let properties = [name, age, height];

let v = grok.shell.newView('demo: inputs extended');

let inputHeight = DG.InputBase.forProperty(height);
inputHeight.format = 'two digits after comma';
let inputs = [
  DG.InputBase.forProperty(name),
  DG.InputBase.forProperty(age),
  inputHeight
];
let content = ui.div();
content.appendChild(ui.inputs(inputs));

for (let n = 0; n < inputs.length; n++) {
  let input = inputs[n];
  input.addCaption(properties[n].name);
  input.setTooltip(`Subject ${properties[n].name.toLowerCase()}`);
  input.nullable = false;
  input.fireChanged();
  let label = ui.divText('');
  label.style.margin = '10px';
  content.appendChild(label);
  input.onChanged.subscribe((value) => {
    label.innerText = `${input.caption}: ${value}`;
  });
}

function reset() {
  for (let input of inputs)
    input.value = properties.find((property) => property.name === input.caption).defaultValue;
}

reset();

content.appendChild(ui.bigButton('Post', function () {
  grok.shell.info(inputs.map((input) => `${input.caption}: ${input.stringValue}`).join('<br>'));
}));
content.appendChild(ui.bigButton('Toggle enabled', function () {
  for (let input of inputs)
    input.enabled = !input.enabled;
}));
content.appendChild(ui.bigButton('Reset', reset));

v.root.appendChild(content);
