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
        this.render();
    }

    onTableAttached() {
    }

    render() {
        $(this.root).empty();
        this.root.appendChild(ui.divText('radio'));

        for (let category of this.column.categories)
            $(`<input>${category}</input>`)
                .prop({type: 'radio' })
                .appendTo(this.root);
    }

    onPropertyChanged(prop) { }
}
