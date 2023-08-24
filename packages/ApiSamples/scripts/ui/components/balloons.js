/// Different way of showing balloons

grok.shell.info('Just a string');
grok.shell.error(ui.div([
  ui.divText('Actionable hint'),
  ui.button('RUN', () => {})]));