// Resolve a DOM element to the u2 control that owns it (DG.U2.Control.forElement) —
// the door for inspectors and automation. Every DG.Widget is a Control, so widgets resolve too.
const widget = DG.Widget.fromRoot(ui.divV([ui.divText('Click me')]));
grok.shell.newView('forElement demo', [widget.root]);

widget.root.addEventListener('click', (e) => {
  const control = DG.U2.Control.forElement(e.target);
  grok.shell.info(`owned by: ${control === widget ? 'the widget' : control}`);
});
