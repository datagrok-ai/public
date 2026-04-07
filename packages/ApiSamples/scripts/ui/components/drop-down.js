// DropDown component examples

grok.shell.newView('Drop-Down', [
  ui.ribbonPanel([
    // Simple menu - object syntax
    ui.dropDown('Menu', {
      'Add': () => grok.shell.info('Add clicked'),
      'Edit': () => grok.shell.info('Edit clicked'),
      'Delete': () => grok.shell.info('Delete clicked'),
    }),

    // Menu with separators - builder function
    ui.dropDown('File', (menu) => {
      menu.item('New', () => grok.shell.info('New'));
      menu.item('Open', () => grok.shell.info('Open'));
      menu.item('Save', () => grok.shell.info('Save'));
      menu.separator();
      menu.item('Exit', () => grok.shell.info('Exit'));
    }),

    // Nested menu - object syntax
    ui.dropDown('Edit', {
      'Undo': () => grok.shell.info('Undo'),
      'Redo': () => grok.shell.info('Redo'),
      'Clipboard': {
        'Operations': {
          'Cut': () => grok.shell.info('Cut'),
          'Copy': () => grok.shell.info('Copy'),
          'Paste': () => grok.shell.info('Paste'),
        },
        'Clipboard operations': () => grok.shell.info('Clipboard operations clicked'),
      },
    }),

    // Icon dropdown
    ui.dropDown(ui.iconFA('plus'), {
      'Add Item': () => grok.shell.info('Add item'),
      'Add Group': () => grok.shell.info('Add group'),
    }),

    // List dropdown
    ui.dropDown('Select', ['Option 1', 'Option 2', 'Option 3'], {
      onItemClick: (item) => grok.shell.info(`Selected: ${item}`)
    }),

    // Custom content dropdown
    ui.dropDown('Custom', () => ui.divV([
      ui.h3('Actions'),
      ui.button('Run Analysis', () => grok.shell.info('Running analysis...')),
      ui.button('Export Data', () => grok.shell.info('Exporting...')),
      ui.element('hr'),
      ui.link('View Help', () => grok.shell.info('Opening help...')),
    ])),

    // Custom content that stays open
    ui.dropDown('Settings', () => ui.divV([
      ui.h3('Options'),
      ui.input.bool('Enable feature', {value: true}),
      ui.input.choice('Mode', {value: 'Auto', items: ['Auto', 'Manual', 'Off']}),
    ]), {closeOnClick: false}),
  ])
]);
