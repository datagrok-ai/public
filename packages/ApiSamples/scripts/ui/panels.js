// TabControl

let v = grok.shell.newView('Tabs');

// Simple API
v.append(ui.div([
    ui.h1('First panel'),
    'Panel content'
]));

v.append(ui.div([
    ui.h1('Second panel'),
    ui.inputs([
        
    ])
]));