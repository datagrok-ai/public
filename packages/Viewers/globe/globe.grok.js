class GlobeViewer extends DG.JsViewer {

    constructor() {
        super();

        // properties
        this.latitude = this.string('latitudeColumnName');
        this.longitude = this.string('longitudeColumnName');
        this.magnitude = this.float('magnitude');
    }

    init() {
        let globe = new DAT.Globe(this.root, {imgDir: this.webRoot});

        let latCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.LATITUDE);
        let lonCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.LONGITUDE);
        let magCol = this.dataFrame.getCol(this.options.magnitude);

        let points = [];
        for (let i = 0; i < this.dataFrame.rowCount; i++) {
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

    onTableAttached() {
        this.init();
        this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
        this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));
        this.render();
    }

    render() {}
}
