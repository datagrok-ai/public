class LeafletViewer extends JsViewer {

    constructor() {
        super();

        // properties
        this.latitude = this.string('latitudeColumnName');
        this.longitude = this.string('longitudeColumnName');
    }

    init() {
        let mapDiv = ui.div([], 'd4-viewer-host');
        this.root.appendChild(mapDiv);
        this.map = L.map(mapDiv).setView([51.505, -0.09], 13);

        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }).addTo(this.map);

        // L.marker([51.5, -0.09]).addTo(this.map)
        //     .bindPopup('A pretty CSS3 popup.<br> Easily customizable.')
        //     .openPopup();
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

    render() {
        let markers = [];
        for (let i = 0; i < this.table.rowCount; i++)
            markers.push(L.marker([this.latitude.get(i), this.longitude.get(i)]));

        const group = L.featureGroup(markers).addTo(this.map);
        this.map.fitBounds(group.getBounds());
    }
}