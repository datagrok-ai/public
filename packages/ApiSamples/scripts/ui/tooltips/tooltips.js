let v = grok.shell.newView('Tooltips demo');

v.root.appendChild(ui.divV([
  ui.tooltip.bind(ui.span(['Static']), 'Simple string'),
  ui.tooltip.bind(ui.span(['Dynamic string']), () => `Function that returns a string. Time: ${new Date()}`),
  ui.tooltip.bind(ui.span(['ddd']), ui.button('Click me', () => grok.shell.info('Clicked'))),
  ui.tooltip.bind(ui.span(['Element']), ui.span([ui.iconFA('search'), 'DOM Element'])),
  ui.tooltip.bind(ui.span(['Markup']), () => ui.markdown('Learn more [here](https://wikipedia.org)')),
]));


//ui.tooltip.bind(ui.span(['Element']), () => ui.button('Click me', () => grok.shell.info('Clicked'))),
