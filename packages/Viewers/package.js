class ViewersPackage extends DG.Package {

    //name: Sankey
    //description: creates a sankey viewer
    //tags: viewer
    //output: viewer result
    sankeyViewer() {
        return new SankeyViewer();
    }

    //name: Sunburst
    //description: creates a sunburst viewer
    //tags: viewer
    //output: viewer result
    sunburstViewer() {
        return new SunburstViewer();
    }

    //name: Globe
    //description: creates a globe viewer
    //tags: viewer
    //output: viewer result
    globeViewer() {
        return new GlobeViewer(this.webRoot + '/globe/');
    }
}
