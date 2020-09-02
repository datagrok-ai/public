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

    //name: flagCellRenderer
    //tags: cellRenderer, cellRenderer-flag
    //meta-cell-renderer-sem-type: flag
    //output: grid_cell_renderer result
    flagCellRenderer() { return new FlagCellRenderer(); }
}


class FlagCellRenderer extends DG.GridCellRenderer {

    get name() { return 'RDKit cell renderer'; }

    get cellType() { return 'flag'; }

    render(g, x, y, w, h, gridCell, cellStyle) {
        g.fillStyle = 'black';
        g.fillText('flag', x, y);
    }
}