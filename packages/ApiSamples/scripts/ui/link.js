// Text links

grok.shell.newView('Links',[
	ui.divV([
	  ui.link('datagrok','https://datagrok.ai','tooltip message'),
	  ui.link('Hello',()=>{grok.shell.info('hello')},'tooltip message')
	])
  ]);