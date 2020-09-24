/**
 * Single-option categorical filter
 * */
class RadioButtonFilter extends DG.Filter {

    constructor() {
        super();
        this.root = ui.divV(null, 'd4-radio-button-filter');
    }

    attach(dFrame) {
        this.dataFrame = DG.toJs(dFrame);
        this.column = this.dataFrame.columns.categorical[0];
        this.dataFrame.onRowsFiltering.subscribe((_) => applyFilter());
        this.render();
    }

    applyFilter() {
        let indexes = this.column.getRawData();
        const checks = this.root.querySelectorAll('input');
        for (let checkIdx = 0; checkIdx < checks.length; checkIdx++) {
            if (!checks[checkIdx].checked) {
                for (let i = 0; i < this.dataFrame.rowCount; i++)
                    if (indexes[i] === checkIdx)
                        this.dataFrame.filter.set(i, false, false);
            }
        }

        this.dataFrame.filter.fireChanged();
    }

    render() {
        $(this.root).empty();
        this.root.appendChild(ui.divText('radio'));

        for (let category of this.column.categories)
            $(`<input>${category}</input>`)
                .on('change', () => this.applyFilter())
                .prop({type: 'radio' })
                .prop({name: category})
                .prop({value: 'yes'})
                .appendTo(this.root);
    }

    onTableAttached() {}
    onPropertyChanged(prop) { }
}
