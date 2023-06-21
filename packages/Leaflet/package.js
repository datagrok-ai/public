class LeafletPackage extends DG.Package {

  //name: Leaflet
  //description: Leaflet map
  //tags: viewer
  //output: viewer result
  //meta.icon: icons/leaflet-viewer.svg
  leafletViewer() {
    return new LeafletViewer();
  }

  //name: Leaflet
  //input: column long { semType: longitude }
  //input: column lat { semType: latitude }
  showOnMap(long, lat) {
    grok.shell.info('map it!');
  }
}