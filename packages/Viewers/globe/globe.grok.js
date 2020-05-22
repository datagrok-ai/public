class GlobeViewer extends grok.JsViewer {

    constructor(webRoot) {
        super();
        this.webRoot = webRoot;
    }

    init() {
        var globe = new DAT.Globe(this.root, {imgDir: this.webRoot});
        var latCol = this.table.getCol(this.options.lat);
        var lonCol = this.table.getCol(this.options.lon);
        var magCol = this.table.getCol(this.options.mag);
        var points = [];
        for (var i = 0; i < this.table.rowCount; i++) {
            points.push(latCol.get(i));
            points.push(lonCol.get(i));
            points.push(magCol.get(i));
        }
        TWEEN.start();
        globe.addData(points, {format: 'magnitude', name: 'Series A', animated: true});
        globe.createPoints();
        new TWEEN.Tween(globe).to({time: 0},500).easing(TWEEN.Easing.Cubic.EaseOut).start();
        globe.animate();
    }

    onFrameAttached(dataFrameHandle) {
        this.options = {lat: 'Latitude', lon: 'Longitude', mag: 'Magnitude'};

        this.table = new DG.DataFrame(dataFrameHandle);
        this.init();

        this.subs.push(this.table.selection.onChanged.subscribe((_) => this.render()));
        this.subs.push(this.table.filter.onChanged.subscribe((_) => this.render()));
        this.render();
    }

    render() {
    }
}
