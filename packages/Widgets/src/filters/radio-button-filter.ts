import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import '../styles/widgets.css';

/**
 * Single-option categorical filter that demonstrates the concept of collaborative filtering:
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */
export class RadioButtonFilter extends DG.Filter {
  static groupId = 0;
  static buttonId = 0;

  constructor() {
    super();
    this.root = ui.divV([], 'd4-radio-button-filter');
    this.subs = [];
  }

  get isFiltering() {return super.isFiltering && this.checkedCategoryId != null;}

  get filterSummary() {
    const categoryIdx = this.checkedCategoryId;
    return categoryIdx == null ? '' : this.column!.getCategory(categoryIdx);
  }

  get checkedCategoryId(): number | null {
    const checkedInput = this.root.querySelector('input[type=\'radio\']:checked');
    return checkedInput == null ? null : parseInt(checkedInput.getAttribute('data-category-id')!);
  }

  attach(dataFrame: DG.DataFrame) {
    super.attach(dataFrame);
    this.column = DG.Utils.firstOrNull(this.dataFrame!.columns.categorical);
    this.columnName = this.column!.name;

    this.render();
  }

  applyState(state: string) {
    super.applyState(state);
    this.render();
  }

  onPropertyChanged(property: DG.Property | null): void {
    console.log(property?.get(this));
  }

  detach() {
    super.detach();
    console.log('radio button filter detached');
  }

  applyFilter() {
    const categoryIdx = this.checkedCategoryId;
    if (categoryIdx == null)
      return;
    const indexes = this.column!.getRawData();
    const filter = this.dataFrame!.filter;
    const rowCount = this.dataFrame!.rowCount;

    for (let i = 0; i < rowCount; i++)
      filter.set(i, indexes[i] === categoryIdx, false);

    this.dataFrame!.filter.fireChanged();
  }

  render() {
    const name = `radio_${RadioButtonFilter.groupId++}`;
    $(this.root).empty();

    for (let i = 0; i < Math.min(20, this.column!.categories.length); i++) {
      const category = this.column!.categories[i];
      const id = `rb_${RadioButtonFilter.buttonId++}`;
      const radioButon = $(`<input type="radio" id="${id}" name="${name}" data-category-id="${i}">`)
        .on('change', () => this.dataFrame!.rows.requestFilter());
      const label = $(`<label for="${id}">${category}</label>`);
      this.root.appendChild(ui.div([radioButon[0], label[0]]));
    }
  }
}
