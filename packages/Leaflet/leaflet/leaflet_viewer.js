class LeafletViewer extends JsViewer {

    init() {
        this.root.style.height = '100px';
        this.map = L.map(this.root).setView([51.505, -0.09], 13);

        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }).addTo(map);

        L.marker([51.5, -0.09]).addTo(map)
            .bindPopup('A pretty CSS3 popup.<br> Easily customizable.')
            .openPopup();
    }

    onFrameAttached(dataFrameHandle) {
        this.table = new DataFrame(dataFrameHandle);
        this.init();

        this.table.selection.onChanged(() => this.render());
        this.table.filter.onChanged(() => this.render());
        this.render();
    }

    render() {

    }
}