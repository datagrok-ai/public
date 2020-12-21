let v = grok.shell.newView('Tooltips demo');

v.root.appendChild(ui.divV([
  ui.tooltip.bind(ui.span(['Static']), 'Parameter is a string parameter'),
  ui.tooltip.bind(ui.span(['Dynamic string']), () => `Function that returns a string. Time: ${new Date()}`),
  ui.tooltip.bind(ui.span(['Element']), ui.span([ui.iconFA('search'), 'DOM Element'])),
  ui.tooltip.bind(ui.span(['Dynamic string']), () => ui.span([ui.iconFA('search'), `Function that returns DOM element. Time: ${new Date()}`])),
]));
