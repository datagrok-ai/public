let v = grok.shell.newView('demo: column inputs');
let t = grok.data.testData('demog', 100);

let predict = ui.columnInput('Predict', t,  t.col('age'));
let features = ui.columnsInput('Features', t);

// events
predict.onInput(() => grok.shell.info(predict.value.map((col) => col.name)) );

let inputs = [predict, features];
let container = ui.div();
v.append(container);
container.appendChild(ui.inputs(inputs));

container.appendChild(ui.bigButton('Build', () => {
    grok.shell.info(inputs.map((i) => `${i.caption}: ${i.stringValue}`).join('<br>'));
}));
