// Categorical and continuous color palettes

let v = grok.shell.newView('palettes');

function getBlock(c) {
  let block = ui.div([DG.Color.toRgb(c), '  hex:', DG.Color.toHtml(c)]);
  block.style.backgroundColor = DG.Color.toRgb(c);
  block.style.color = DG.Color.toRgb(DG.Color.getContrastColor(c));
  return block;
}

v.appendAll([
  ui.h1('Categorical palette with contrast text color'),
  ui.div(DG.Color.categoricalPalette.map(getBlock)),
  ui.h1('Get contrast color from white'),
  ui.div(getBlock(DG.Color.getContrastColor(DG.Color.white)))
]);