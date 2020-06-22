class ViewersPackage extends DG.Package {

    //name: Sankey
    //description: Creates a sankey viewer
    //tags: viewer
    //output: viewer result
    sankeyViewer() {
        return new SankeyViewer();
    }

    //name: Sunburst
    //description: Creates a sunburst viewer
    //tags: viewer
    //output: viewer result
    sunburstViewer() {
        return new SunburstViewer();
    }

    //name: Globe
    //description: Creates a globe viewer
    //tags: viewer
    //output: viewer result
    globeViewer() {
        return new GlobeViewer(this.webRoot + '/globe/');
    }
}
