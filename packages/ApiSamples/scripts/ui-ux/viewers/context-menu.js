// Add an item to the viewer's context menu
// See also ../components/menu.js

grok.shell
  .addTableView(grok.data.demo.demog())
  .barChart()
  .onContextMenu.subscribe(menu => menu.item('Foo', () => grok.shell.info('Foo!')));