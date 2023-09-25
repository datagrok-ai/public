// Creating custom views

let v = grok.shell.newView('list');

v.root.appendChild(ui.h1('List'));

v.root.appendChild(ui.list([
  'element 1',
  grok.shell.user
]));
