// Creating custom views

let v = grok.newView('list');

v.root.appendChild(ui.h1('List'));

v.root.appendChild(ui.list([
    'element 1',
    User.current
]));
