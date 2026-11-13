//api: ui.colorPicker, ui.showColorPicker, DG.Color.opaqueIfNoAlpha
// Color pickers with an alpha channel: the picker gets an alpha strip next to the hue slider,
// and the reported color keeps its alpha (ARGB, so DG.Color.a reads the opacity).
// A color written without an alpha (0xRRGGBB) is drawn opaque: DG.Color.opaqueIfNoAlpha says so.

let v = grok.shell.newView('color picker');
let swatch = ui.div([], {style: {width: '40px', height: '40px', border: '1px solid var(--grey-3)'}});
let readout = ui.divText('');
let report = (c) => readout.textContent = `${DG.Color.toHtml(c)}, alpha ${DG.Color.a(c)}`;

v.appendAll([
  ui.h1('Click the swatch'),
  ui.colorPicker(DG.Color.fromHtml('#ff000080'), report, swatch, null, null, true),
  readout,
  ui.divText(`0xFF0000 is drawn with alpha ${DG.Color.a(DG.Color.opaqueIfNoAlpha(0xFF0000))}`),
  ui.button('Open a standalone picker', () => ui.showColorPicker(DG.Color.fromHtml('#0000ff40'), report, null, null, true)),
]);
