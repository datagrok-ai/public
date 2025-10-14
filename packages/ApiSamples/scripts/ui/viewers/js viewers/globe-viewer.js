// https://datagrok.ai/help/visualize/viewers/globe

let view = grok.shell.addTableView(grok.data.demo.geo());

view.addViewer('globe', {
  'latitudeColumnName': 'lat',
  'longitudeColumnName': 'lng',
  'magnitudeColumnName': 'lat',
  'colorByColumnName': 'value',
  'pointRadius': 15,
  'pointAltitude': 50,
  'autorotation': false
});