class LeafletViewer extends JsViewer {
    constructor() {
        super();

        // properties
        this.latitude = this.string('latitudeColumnName');
        this.longitude = this.string('longitudeColumnName');
        this.renderType = this.string('renderType', 'heat map');
        this.markerSize = this.int('markerSize', 10);
        this.markerOpacity = this.float('markerOpacity', 0.8);

        this.getProperty('renderType').choices = ['markers', 'heat map'];

        this.layers = [];
        this.coordinates = [];
    }

    init() {
        let mapDiv = ui.div([], 'd4-viewer-host');
        this.root.appendChild(mapDiv);
        this.map = L.map(mapDiv).setView([51.505, -0.09], 13);

        this.tileLayer = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }).addTo(this.map);
    }

    onFrameAttached(dataFrameHandle) {
        this.table = new DataFrame(dataFrameHandle);
        console.log('foo');
        this.init();

        this.latitude = this.table.columns.bySemType(SEMTYPE_LATITUDE);
        this.longitude = this.table.columns.bySemType(SEMTYPE_LONGITUDE);

        this.subs.push(debounce(this.table.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(debounce(this.table.filter.onChanged, 50).subscribe((_) => this.render()));

        this.render(true);
    }

    onPropertyChanged(prop) {
        this.render();
    }

    onSizeChanged(w, h) { this.map.invalidateSize(); }

    render(fit = false) {
        for (const layer of this.layers)
            this.map.removeLayer(layer);
        this.layers.length = 0;

        this.getCoordinates();

        if (this.renderType === 'heat map')
            this.renderHeat();
        else if (this.renderType === 'markers')
            this.renderMarkers();

        if (fit)
            this.map.fitBounds(this.coordinates);
    }

    getCoordinates() {
        this.coordinates.length = 0;
        let indexes = this.table.filter.getSelectedIndexes();
        let lat = this.latitude.getRawData();
        let lon = this.longitude.getRawData();

        for (let i = 0; i < indexes.length; i++)
            this.coordinates.push([lat[indexes[i]], lon[indexes[i]]]);
    }

    renderHeat() {
        this.layers.push(L.heatLayer(this.coordinates, {radius: this.markerSize}).addTo(this.map));
    }

    renderMarkers() {
        let markerOptions = {
            radius: this.markerSize,
            fillColor: "#ff7800",
            color: "#000",
            weight: 1,
            opacity: 1,
            fillOpacity: this.markerOpacity
        };

        let markers = [];
        for (let i = 0; i < this.coordinates.length; i++)
            markers.push(L.circleMarker(this.coordinates[i], markerOptions));
        this.layers.push(L.featureGroup(markers).addTo(this.map));
    }
}