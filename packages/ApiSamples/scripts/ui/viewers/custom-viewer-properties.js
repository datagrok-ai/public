let view = grok.shell.addTableView(grok.data.demo.geo(10000));
let leaflet = view.addViewer('Leaflet');
leaflet.setOptions({
  latitudeColumnName: 'lat',
  longitudeColumnName: 'lng'
});

grok.shell.info(leaflet.props.getProperties().map((p) => p.name).join(', '));