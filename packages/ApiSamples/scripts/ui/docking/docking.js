// Different options to dock windows

// dock to the right of the root
let node1 = grok.shell.dockManager.dock(ui.divText('first'), 'right', null, 'First');

// to the bottom of the previously created node
let node2 = grok.shell.dockManager.dock(ui.divText('second'), 'down', node1, 'Second');

// inside node 2 (tab group will be created)
let node3 = grok.shell.dockManager.dock(ui.divText('third'), 'fill', node2, 'Third');

// to the bottom of the current document
grok.shell.dockManager.dock(ui.divText('bottom'), 'down', grok.shell.v.dockNode, 'Bottom');

// handling events
grok.shell.dockManager.onClosed.subscribe((e) => grok.shell.info(e));
