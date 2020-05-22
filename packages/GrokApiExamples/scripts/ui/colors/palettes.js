// Categorical and continuous color palettes

let v = grok.shell.newView('palettes');

function getBlock(c) {
    let block = ui.divText(ui.Color.toRgb(c));
    block.style.backgroundColor = ui.Color.toRgb(c);
    block.style.color = ui.Color.toRgb(ui.Color.getContrastColor(c));
    return block;
}

v.appendAll([
    ui.h1('Categorical palette with contrast text color'),
    ui.div(ui.Color.categoricalPalette.map(getBlock)),
    ui.h1('Category colors (looping over the palette)'),
    ui.div(grok.Utils.identity(30).map(ui.Color.getCategoricalColor).map(getBlock))
]);