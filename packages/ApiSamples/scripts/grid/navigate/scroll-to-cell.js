// Scroll to a particular cell.

let grid = grok.shell.addTableView(grok.data.demo.randomWalk(100, 100)).grid;

ui.dialog('Scroll')
  .add(ui.button('Jump to #100x80', () => grid.scrollToCell('#55', 99)))
  .show();