let grid = grok.shell.addTableView(grok.data.demo.randomWalk(100, 100)).grid;

let vert = grid.vertScroll;

function info() {
  grok.shell.info(`Max: ${vert.max}, min: ${vert.min}, min range: ${vert.minRange}, max range: ${vert.maxRange}`)
}

ui.dialog('Scroll bars')
  .add(ui.button('Get', info))
  .show();