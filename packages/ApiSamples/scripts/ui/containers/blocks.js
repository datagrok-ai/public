// Blocks

let view = grok.shell.newView('Block', [
  ui.block([ui.panel('100% block width')], 'myclass'),
  //two blocks next to each other
  ui.block75([ui.panel('75% block width')], 'myclass'),  
  ui.block25([ui.panel('25% block width')], 'myclass'),
  //three blocks next ot each other
  ui.block50([ui.panel('50% block width')], 'myclass'),
  ui.block25([ui.panel('25% block width')], 'myclass'),
  ui.block25([ui.panel('25% block width')], 'myclass'),
  //four blocks next ot each other
  ui.block25([ui.panel('25% block width')], 'myclass'),
  ui.block25([ui.panel('25% block width')], 'myclass'),
  ui.block25([ui.panel('25% block width')], 'myclass'),
  ui.block25([ui.panel('25% block width')], 'myclass')
]);
  
$('.myclass').css('border', '1px dashed red');