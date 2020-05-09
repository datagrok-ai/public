class LeafletViewer extends JsViewer {

    constructor() {
        super();

        // properties
        this.latitude = this.string('latitudeColumnName');
        this.longitude = this.string('longitudeColumnName');
        this.markerSize = this.int('markerSize', 10);

        this.layers = [];
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

        this.table.selection.onChanged(() => this.render());
        this.table.filter.onChanged(() => this.render());
        this.render();
    }

    onPropertyChanged(prop) {
        this.render();
    }

    onSizeChanged(w, h) { this.map.invalidateSize(); }

    render() {
        for (const layer of this.layers)
            this.map.removeLayer(layer);

        this.renderHeat();
    }

    renderHeat() {
        let coordinates = [];
        for (let i = 0; i < this.table.rowCount; i++)
            coordinates.push([this.latitude.get(i), this.longitude.get(i)]);

        this.layers.push(L.heatLayer(coordinates, {radius: this.markerSize}).addTo(this.map));
    }

    renderMarkers() {
        var markerOptions = {
            radius: this.markerSize,
            fillColor: "#ff7800",
            color: "#000",
            weight: 1,
            opacity: 1,
            fillOpacity: 0.8
        };

        let markers = [];
        for (let i = 0; i < this.table.rowCount; i++)
            markers.push(L.circleMarker([this.latitude.get(i), this.longitude.get(i)], markerOptions));

        const group = L.featureGroup(markers).addTo(this.map);
        this.map.fitBounds(group.getBounds());
    }
}