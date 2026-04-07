// Typography

grok.shell.newView('Typography', [
  ui.h1('Header h1'),
  ui.h2('Header h2'),
  ui.h3('Header h3'),
  ui.p('paragraph'),
  ui.span(['span 1', 'span 2', ui.divText('span 3')]),
  ui.label('label'),
  ui.divV([
    ui.link('link string', 'https://datagrok.ai', 'tooltip message'),
    ui.link('Link func', ()=>{grok.shell.info('hello');}, 'tooltip message'),
  ]),
  ui.divText('div text'),
  ui.inlineText(['Inline ', ui.link('text', ()=>{grok.shell.info('');}, 'click me', '')])   
]);