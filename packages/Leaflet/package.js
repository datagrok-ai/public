class LeafletPackage extends grok.Package {

    //name: Leaflet
    //description: Leaflet map
    //tags: viewer
    //output: viewer result
    leafletViewer() {
        return new LeafletViewer();
    }

    //name: Leaflet
    //input: column long { semType: longitude }
    //input: column lat { semType: latitude }
    showOnMap(long, lat) {
        grok.balloon.info('map it!');
    }
}