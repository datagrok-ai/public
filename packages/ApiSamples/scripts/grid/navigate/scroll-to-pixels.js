// Scroll to pixels.

const grid = grok.shell.addTableView(grok.data.demo.randomWalk(100, 100)).grid;
let cntr = 0;

ui.dialog('Scroll')
  .add(ui.button('Scroll by 10 pixels', () => {
    cntr += 10;
    grid.scrollToPixels(cntr, 0);
  }))
  .show();