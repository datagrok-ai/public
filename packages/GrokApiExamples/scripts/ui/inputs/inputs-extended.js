let name = grok.Property.string('Name', null, null, '');
let age = grok.Property.int('Age', null, null, 30);
let height = grok.Property.float('Height', null, null, 175.5);
var properties = [name, age, height];

let v = grok.shell.newView('demo: inputs extended');

var inputHeight = ui.InputBase.forProperty(height);
inputHeight.format = 'two digits after comma';
var inputs = [
    ui.InputBase.forProperty(name),
    ui.InputBase.forProperty(age),
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
    input.onChanged(function () {
        label.innerText = `${input.caption}: ${input.stringValue}`;
    });
}

var reset = function () {
    for (var input of inputs)
        input.value = properties.find((property) => property.name === input.caption).defaultValue;
}
reset();

content.appendChild(ui.bigButton('Post', function () {
    grok.shell.balloon.info(inputs.map((input) => `${input.caption}: ${input.stringValue}`).join('<br>'));
}));
content.appendChild(ui.bigButton('Toggle enabled', function () {
    for (let input of inputs)
        input.enabled = !input.enabled;
}));
content.appendChild(ui.bigButton('Reset', reset));

v.root.appendChild(content);
grok.addView(v);
