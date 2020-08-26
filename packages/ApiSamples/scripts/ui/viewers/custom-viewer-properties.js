let view = grok.shell.addTableView(grok.data.demo.geo(10000));
view.addViewer('Leaflet').setOptions({
    latitudeColumnName: 'lat',
    longitudeColumnName: 'lon'
});
