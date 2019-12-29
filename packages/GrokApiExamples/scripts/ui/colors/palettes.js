// Categorical and continuous color palettes

let v = grok.newView('palettes');

function getBlock(c) {
    let block = ui.divText(Color.toRgb(c));
    block.style.backgroundColor = Color.toRgb(c);
    block.style.color = Color.toRgb(Color.getContrastColor(c));
    return block;
}

v.appendAll([
    ui.h1('Categorical palette with contrast text color'),
    ui.div(Color.categoricalPalette.map(getBlock)),
    ui.h1('Category colors (looping over the palette)'),
    ui.div(Utils.identity(30).map(Color.getCategoricalColor).map(getBlock))
]);