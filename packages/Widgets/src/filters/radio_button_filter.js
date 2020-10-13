let groupId = 0;
let buttonId = 0;

/**
 * Single-option categorical filter that demonstrates the concept of collaborative filtering:
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */
class RadioButtonFilter extends DG.Filter {

    constructor() {
        super();
        this.root = ui.divV(null, 'd4-radio-button-filter');
    }

    attach(dFrame) {
        this.dataFrame = DG.toJs(dFrame);
        this.column = this.dataFrame.columns.categorical[0];
        this.dataFrame.onRowsFiltering.subscribe((_) => this.applyFilter());
        this.render();
    }

    applyFilter() {
        let indexes = this.column.getRawData();
        let checked = this.root.querySelector("input[type='radio']:checked");
        let categoryIdx = parseInt(checked.getAttribute('data-category-id'));

        for (let i = 0; i < this.dataFrame.rowCount; i++)
            if (indexes[i] !== categoryIdx)
                this.dataFrame.filter.set(i, false, false);

        this.dataFrame.filter.fireChanged();
    }

    render() {
        let name = `radio_${groupId++}`;
        $(this.root).empty();

        for (let i = 0; i < this.column.categories.length; i++) {
            let category = this.column.categories[i];
            let id = `rb_${buttonId++}`;
            let radioButon = $(`<input type="radio" id="${id}" name="${name}" data-category-id="${i}">`)
                .on('change', () => this.dataFrame.rows.requestFilter())
            let label = $(`<label for="${id}">${category}</label>`)
            this.root.appendChild(ui.div([radioButon[0], label[0]]));
        }
    }
}
