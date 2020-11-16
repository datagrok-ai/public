// You can use custom HTML as the current object

let root = ui.div();
root.appendChild(ui.h2('Custom HTML header'));
root.appendChild(ui.span(['Here you can add whatever you want!']));

grok.shell.o = root;