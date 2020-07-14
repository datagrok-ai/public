class MultiFormViewer extends DG.JsViewer {
    constructor() {
        super();

        // properties
        this.fields = this.stringList('fieldsColumnNames');

        // init
        let dis = this;
        this.root.classList.add('d4-multi-form');
        this.columnHeadersDiv = ui.divV([]);
        this.virtualView = ui.virtualView(0, function(i) { return dis.renderForm(i); });
        this.root.appendChild(ui.divH([this.columnHeadersDiv, this.virtualView.root]));
    }

    onTableAttached() {
        if (this.fields == null)
            this.fields = this.dataFrame.columns.names();

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 20).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

        this.render();
    }

    onPropertyChanged(prop) {
        this.render();
    }

    renderHeader() {
        return ui.divV(this.fields.map((name) => ui.divText(name, 'd4-multi-form-column-name')));
    }

    renderForm(row) {
        console.log('' + row);
        console.log('' + this.fields);
        return ui.divV(
            this.fields.map((name) =>
                ui.divText(this.dataFrame.get(name, row), 'd4-multi-form-value')));
    }

    render() {
        let dis = this;
        let indexes = this.dataFrame.selection.getSelectedIndexes();
        ui.empty(this.columnHeadersDiv);
        this.columnHeadersDiv.appendChild(this.renderHeader());
        this.virtualView.setData(indexes.length, function(i) { return dis.renderForm(i); });
    }
}
