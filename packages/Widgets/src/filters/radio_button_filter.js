/**
 * Single-option categorical filter
 * */
class RadioButtonFilter extends DG.Filter {
    constructor() {
        super();
        this.root = ui.div(null, 'd4-radio-button-filter');
    }

    onTableAttached() {
        this.render();
        this.column = this.dataFrame.columns.categorical[0];
    }

    render() {
        $(this.root).empty();
        for (let category of this.column.categories)
            this.root.appendChild($('<input>').prop({
                type: 'radio',
                value: category
            }));
    }

    onPropertyChanged(prop) { }
}
