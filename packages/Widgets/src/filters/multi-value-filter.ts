import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

export class MultiValueFilter extends DG.Filter {
  valuesHost: HTMLDivElement;
  buttonId = 0;

  constructor() {
    super();
    this.valuesHost = ui.divV([]);
    //this.typeSelect = ui.choiceInput('Type', 'any', ['any', 'all']);
    this.root = ui.divV([this.valuesHost]);
    this.subs = [];
  }

  //@ts-ignore: it doesn't return boolean
  get isFiltering() {return this.root.querySelectorAll('input[type=\'checkbox\']:not(:checked)');}

  get filterSummary() {return `${this.getSelectedInputs().length}`;}

  attach(dataFrame: DG.DataFrame) {
    super.attach(dataFrame);
    this.column = this.dataFrame!.columns.byIndex(0);

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
    console.log('multi filter detached');
  }

  /** @return {NodeListOf<HTMLInputElement>} */
  getSelectedInputs() {
    return this.root.querySelectorAll('input[type=\'checkbox\']:checked');
  }

  applyFilter() {
    const separator = ' | ';
    const checkedValues = new Set();
    const checkedNodes = this.getSelectedInputs();
    checkedNodes.forEach((check) => checkedValues.add(check.getAttribute('name')));

    const filter = this.dataFrame!.filter;
    const rowCount = this.dataFrame!.rowCount;

    for (let i = 0; i < rowCount; i++) {
      const vv: string[] = this.column!.get(i).split(separator);
      const contains = vv.some((v) => checkedValues.has(v));

      if (!contains)
        filter.set(i, false, false);
    }

    this.dataFrame!.filter.fireChanged();
  }

  render() {
    const values = new Set();
    const separator = ' | '; //;this.column.tags[DG.TAGS.MULTI_VALUE_SEPARATOR];

    for (let i = 0; i < this.column!.length; i++) {
      const vv = this.column!.get(i).split(separator);
      //let vv = s.map((s) => s.trim());
      for (const v of vv) {
        if (v !== '')
          values.add(v);
      }
    }

    $(this.valuesHost).empty();

    for (const v of values) {
      const id = `rb_${this.buttonId++}`;
      const check = $(`<input type="checkbox" id="${id}" name="${v}">`)
        .on('change', () => this.dataFrame!.rows.requestFilter());
      const label = $(`<label for="${id}">${v}</label>`);
      this.valuesHost.appendChild(ui.div([check[0], label[0]]));
    }
  }
}
