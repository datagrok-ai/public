// Box

let d = ui.div();
for (let i = 0; i < 100; i++)
  d.append(ui.p('More Text'));
let box = ui.box(d);

box.style.border = '1px dashed red'; 

grok.shell.newView('Box', [box]);