let v = gr.newView('Tooltips demo');

v.root.appendChild(ui.divV([
    ui.tooltip(ui.span(['Static']), 'Parameter is a string parameter'),
    ui.tooltip(ui.span(['Dynamic string']), () => `Function that returns a string. Time: ${new Date()}`),
    ui.tooltip(ui.span(['Element']), ui.span([ui.iconFA('search'), 'DOM Element'])),
    ui.tooltip(ui.span(['Dynamic string']), () => ui.span([ui.iconFA('search'), `Function that returns DOM element. Time: ${new Date()}`])),
]));
// String tooltip
