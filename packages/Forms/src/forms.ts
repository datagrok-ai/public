import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Observable} from 'rxjs';
import {filter} from 'rxjs/operators';


import '../css/forms.css';


export class FormsViewer extends DG.JsViewer {
  fieldsColumnNames: string[];
  colorCode: boolean;
  showCurrentRow: boolean;
  showMouseOverRow: boolean;
  indexes: number[];
  columnHeadersDiv: HTMLDivElement;
  virtualView: DG.VirtualView;

  constructor() {
    super();

    // properties
    this.fieldsColumnNames = this.addProperty('fieldsColumnNames', DG.TYPE.COLUMN_LIST);
    this.colorCode = this.bool('colorCode', true);
    this.showCurrentRow = this.bool('showCurrentRow', true);
    this.showMouseOverRow = this.bool('showMouseOverRow', true);

    //fields
    this.indexes = [];

    // init
    this.root.classList.add('d4-multi-form');
    this.columnHeadersDiv = ui.div([], 'd4-multi-form-header');
    this.virtualView = ui.virtualView(0, (i: number) => this.renderForm(i), false, 1);
    this.root.appendChild(ui.divH([this.columnHeadersDiv, this.virtualView.root]));
  }

  onTableAttached() {
    if (this.fieldsColumnNames === null)
      this.fieldsColumnNames = this.dataFrame.columns.names();

    const sub = (stream: Observable<unknown>, action: Function) => {
      this.subs.push(DG.debounce(stream, 50).subscribe((_) => action()));
    };

    sub(this.dataFrame.selection.onChanged, () => this.render());
    sub(this.dataFrame.filter.onChanged, () => this.render());

    sub(this.dataFrame.onCurrentRowChanged.pipe(filter((_) => this.showCurrentRow)),
      () => this.virtualView.refreshItem(this.currentRowPos!));

    sub(this.dataFrame.onMouseOverRowChanged.pipe(filter((_) => this.showMouseOverRow)),
      () => this.virtualView.refreshItem(this.mouseOverPos!));

    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    this.render();
  }

  renderHeader() {
    ui.empty(this.columnHeadersDiv);
    const form = this.renderForm(0);
    form.classList.add('temp');
    document.body.appendChild(form);

    for (const name of this.fieldsColumnNames) {
      const columnLabel = ui.bind(this.dataFrame.col(name), ui.divText(name, 'd4-multi-form-column-name'));
      const formField = form.querySelector('[column="' + name + '"]');
      columnLabel.style.top = `${(formField as HTMLElement).offsetTop + 10}px`;
      this.columnHeadersDiv.appendChild(columnLabel);
    }

    form.remove();
  }

  get currentRowPos() {return this.showCurrentRow ? 0 : null;}
  get mouseOverPos() {return this.showMouseOverRow ? (this.showCurrentRow ? 1 : 0) : null;}

  renderForm(row: number) {
    if (this.showCurrentRow && row === this.currentRowPos)
      row = this.dataFrame.currentRowIdx;
    else if (this.showMouseOverRow && row === this.mouseOverPos)
      row = this.dataFrame.mouseOverRowIdx;
    else
      row = this.indexes[row - (this.showCurrentRow ? 1 : 0) - (this.showMouseOverRow ? 1 : 0)];

    const form = ui.divV(
      this.fieldsColumnNames.map((name) => {
        if (row === -1)
          return ui.div();

        const input = DG.InputBase.forColumn(this.dataFrame.col(name)!);
        input.input.setAttribute('column', name);
        input.value = this.dataFrame.get(name, row);

        if (this.colorCode) {
          const grid = this.view.grid;
          const color = grid.cell(name, row).color;
          input.input.style.color = DG.Color.toHtml(DG.Color.getContrastColor(color));
          input.input.style.backgroundColor = DG.Color.toHtml(color);
        }

        return input.input;
      },
      ), 'd4-multi-form-form');

    form.onclick = () => this.dataFrame.currentRowIdx = row;
    form.onmouseenter = () => this.dataFrame.mouseOverRowIdx = row;
    form.onmouseleave = () => this.dataFrame.mouseOverRowIdx = -1;
    return form;
  }

  render() {
    this.indexes = this.dataFrame.selection.trueCount > 0 ?
      Array.from(this.dataFrame.selection.getSelectedIndexes()) : [];

    ui.empty(this.columnHeadersDiv);

    this.renderHeader();

    this.virtualView.setData(
      this.indexes.length + (this.showCurrentRow ? 1 : 0) + (this.showMouseOverRow ? 1 : 0),
      (i: number) => this.renderForm(i));
  }
}
