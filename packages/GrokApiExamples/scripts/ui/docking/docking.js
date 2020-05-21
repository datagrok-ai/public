// Docking an arbitrary element in the platform

let e = document.createElement('DIV');
e.innerText = 'This element has been created in JavaScript';
grok.shell.dockElement(e, 'JS', 'left', 0.5);