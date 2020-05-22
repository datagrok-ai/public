// Docking an arbitrary element in the platform

// dock to the right of the root
grok.shell.dockManager.dock(ui.divText('element 1'), 'right', null, 'Title');

// dock to the right of the document container
grok.shell.dockManager.dock(ui.divText('element 2'), 'right', grok.shell.dockManager.documentContainer, 'Title');

// floating window
grok.shell.dockManager.dock(ui.divText('Floating'), 'right')
    .container.undock();