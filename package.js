class SunburstPackage extends DG.Package {

    //name: Sunburst
    //description: Sunburst map
    //tags: viewer
    //output: viewer result
    sunburstViewer() {
        return new SunburstViewer();
    }

    //name: Sunburst
    //input: column long { semType: longitude }
    //input: column lat { semType: latitude }
    showSunburst(level1column, level2column, level3column) {
        grok.shell.info('Sunburst');
    }
}
