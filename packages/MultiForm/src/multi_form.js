class MultiFormViewer extends DG.JsViewer {
    constructor() {
        super();

        // properties
        this.fields = this.stringList('fieldsColumnNames');

        // init
        this.root.classList.add('d4-multi-form');
        this.columnHeadersDiv = ui.divV([]);
        this.virtualView = ui.virtualView(0, renderItem);
        this.root.append(ui.divH([this.columnHeadersDiv, this.virtualView.root]));
    }

    onTableAttached() {
        this.init();
        this.fields ??= this.dataFrame.columns.names();

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 20).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

        this.render(true);
    }

    onPropertyChanged(prop) {
        this.render();
    }

    renderHeader() {
        return ui.divV(this.fields.map((name) => ui.divText(name, 'd4-multi-form-column-name')));
    }

    renderForm(row) {
        return ui.divV(
            this.fields.map((name) =>
                ui.divText(this.dataFrame.get(name, row), 'd4-multi-form-value')));
    }

    render() {
        let indexes = this.dataFrame.selection.getSelectedIndexes();
        ui.empty(this.columnHeadersDiv);
        this.columnHeadersDiv.appendChild(this.renderHeader());
        this.virtualView.setData(indexes.length, this.renderForm)
    }
}
