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
    this.subs = [];
  }

  get isFiltering() { return true; }

  get filterSummary() { return this.column.getCategory(this.checkedCategoryId); }

  get checkedCategoryId() {
    let checkedInput = this.root.querySelector("input[type='radio']:checked");
    return parseInt(checkedInput.getAttribute('data-category-id'));
  }

  attach(dataFrame) {
    super.attach(dataFrame);
    this.column = DG.Utils.firstOrNull(this.dataFrame.columns.categorical);
    this.columnName = this.column.name;

    this.render();
  }

  applyState(state) {
    super.applyState(state);
    this.render();
  }

  detach() {
    super.detach();
    console.log('radio button filter detached');
  }

  applyFilter() {
    let indexes = this.column.getRawData();
    let categoryIdx = this.checkedCategoryId;
    const filter = this.dataFrame.filter;
    const rowCount = this.dataFrame.rowCount;

    for (let i = 0; i < rowCount; i++)
      if (indexes[i] !== categoryIdx)
        filter.set(i, false, false);

    this.dataFrame.filter.fireChanged();
  }

  render() {
    let name = `radio_${groupId++}`;
    $(this.root).empty();

    for (let i = 0; i < Math.min(20, this.column.categories.length); i++) {
      let category = this.column.categories[i];
      let id = `rb_${buttonId++}`;
      let radioButon = $(`<input type="radio" id="${id}" name="${name}" data-category-id="${i}">`)
        .on('change', () => this.dataFrame.rows.requestFilter())
      let label = $(`<label for="${id}">${category}</label>`)
      this.root.appendChild(ui.div([radioButon[0], label[0]]));
    }
  }
}
