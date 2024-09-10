let v = grok.shell.newView('demo: column inputs');
let t = grok.data.testData('demog', 100);

let predict = ui.input.column('Predict', {value: t.col('age'), table: t});
let features = ui.input.columns('Features', {table: t});
features.value = [t.col('height'), t.col('weight')];

// events
features.onInput.subscribe(() => grok.shell.info(features.value.map((col) => col.name).join()));

let inputs = [predict, features];
let container = ui.div();
v.append(container);
container.appendChild(ui.inputs(inputs));

container.appendChild(ui.bigButton('Build', () => {
  grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));