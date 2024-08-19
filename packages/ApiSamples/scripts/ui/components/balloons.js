/// Different way of showing balloons

grok.shell.info('Just a string');
grok.shell.warning('Warning ballon');
grok.shell.error(ui.div([
  ui.divText('Actionable hint'),
  ui.bigButton('RUN', () => {}),
  ui.button('Cancel', ()=>{})
]));