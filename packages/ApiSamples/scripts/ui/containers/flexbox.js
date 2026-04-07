// Flexbox

grok.shell.newView('Flexbox', [
  ui.divH([
    ui.span(['item 1']),
    ui.span(['item 2']),
    ui.span(['item 3'])
  ], 'myclass'), //rows
  ui.divV([
    ui.span(['item 1']),
    ui.span(['item 2']),
    ui.span(['item 3'])
  ], 'myclass'), //collumns
  ui.divV([
    ui.span(['item1']),
    ui.divH([
      ui.span(['item 2']),
      ui.span(['item 3'])
    ])
  ], 'myclass'),
]);

$('.myclass').css('border', '1px dashed red');