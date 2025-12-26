// Buttons

grok.shell.newView('Buttons', [
  ui.button('Regular button', ()=> grok.shell.info('clicked')),
  ui.button(ui.iconFA('info'), ()=> grok.shell.info('clicked')),
  ui.bigButton('Big button', ()=> grok.shell.info('clicked'))
]);  