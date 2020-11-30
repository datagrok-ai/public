// Panels

let v = grok.shell.newView('Tabs');

// Simple API
v.append(ui.div([
    ui.h1('First panel'),
    'Panel content'
]));

v.append(ui.div([
    ui.h1('Second panel'),
    ui.inputs([
        ui.stringInput('Name', 5),
        ui.choiceInput('Algorithm', 'SPE', ['SPE', 'PCA']),
        ui.intInput('Number of steps', 5)
    ]),
    ui.bigButton('RUN')
]));