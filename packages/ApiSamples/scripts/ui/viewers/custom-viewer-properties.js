const view = grok.shell.addTableView(grok.data.demo.geo(1000));
const map = view.addViewer('Map');
map.setOptions({
  latitudeColumnName: 'lat',
  longitudeColumnName: 'lng',
});

grok.shell.info(map.props.getProperties().map((p) => p.name).join(', '));