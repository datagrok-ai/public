//api: DG.AnnotationRegionsHelper.getFillColor
// The color a viewer fills an annotation region with: `fill` is the color with its alpha as the opacity;
// regions saved before it keep `fillColor` and `opacity` (30% when missing).

const regions = [
  {fill: 0x80FF0000},
  {fill: '#ff0000'},
  {fillColor: 0xFFFF0000},
  {fillColor: 0xFFFF0000, opacity: 60},
];

for (const region of regions) {
  const fill = DG.AnnotationRegionsHelper.getFillColor(region);
  grok.shell.info(`${JSON.stringify(region)}: ${DG.Color.toHtml(fill)}, alpha ${DG.Color.a(fill)}`);
}
