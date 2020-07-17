class Sunburst2Package extends DG.Package {

    //name: Sunburst2
    //description: Sunburst chart
    //tags: viewer
    //output: viewer result
    sunburst2Viewer() {
        return new Sunburst2Viewer();
    }

    // //name: Sunburst2
    // //input: column level1 { semType: text }
    // //input: column level2 { semType: text }
    // //input: column level3 { semType: text }
    // showSunburst2(level1, level2, level3) {
    //     grok.shell.info('Sunburst2');
    // }
}
