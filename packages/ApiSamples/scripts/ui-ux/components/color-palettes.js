// Categorical and continuous color palettes

let v = grok.shell.newView('palettes');

function getBlock(c) {
  let block = ui.divText(DG.Color.toRgb(c));
  block.style.backgroundColor = DG.Color.toRgb(c);
  block.style.color = DG.Color.toRgb(DG.Color.getContrastColor(c));
  return block;
}

v.appendAll([
  ui.h1('Categorical palette with contrast text color'),
  ui.div(DG.Color.categoricalPalette.map(getBlock)),
  ui.h1('Category colors (looping over the palette)'),
  ui.div(DG.utils.identity(30).map(DG.Color.getCategoricalColor).map(getBlock))
]);
