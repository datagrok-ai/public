import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Observable} from 'rxjs';
import {filter} from 'rxjs/operators';
import '../../css/forms.css';


export class FormsViewer extends DG.JsViewer {
  get type(): string { return 'FormsViewer'; }

  fieldsColumnNames: string[];
  colorCode: boolean;
  showCurrentRow: boolean;
  showMouseOverRow: boolean;
  showSelectedRows: boolean;
  showFixedRows = false;
  indexes: number[];
  columnHeadersDiv: HTMLDivElement;
  virtualView: DG.VirtualView;


  set dataframe(df: DG.DataFrame) {
    this.dataFrame = df;
    this.onTableAttached();
  }

  set fixedRowNumbers(idxs: number[]) {
    this.indexes = idxs;
    this.showCurrentRow = false;
    this.showMouseOverRow = false;
    this.showSelectedRows = false;
    this.showFixedRows = true;
    this.render();
  }

  set columns(colNames: string[]) {
    this.fieldsColumnNames = colNames;
    this.render();
  }

  constructor() {
    super();

    // properties
    this.fieldsColumnNames = this.addProperty('fieldsColumnNames', DG.TYPE.COLUMN_LIST);
    this.colorCode = this.bool('colorCode', true);
    this.showCurrentRow = this.bool('showCurrentRow', true);
    this.showMouseOverRow = this.bool('showMouseOverRow', true);
    this.showSelectedRows = this.bool('showSelectedRows', true);

    //fields
    this.indexes = [];

    // init
    this.root.classList.add('d4-multi-form');
    this.columnHeadersDiv = ui.div([], 'd4-multi-form-header');
    this.virtualView = ui.virtualView(0, (i: number) => this.renderForm(i), false, 1);
    const formWithHeaderDiv = ui.divH([this.columnHeadersDiv, this.virtualView.root]);
    formWithHeaderDiv.setAttribute('style', 'overflow:visible!important');
    this.root.appendChild(formWithHeaderDiv);
  }

  onTableAttached() {
    if (this.fieldsColumnNames === null)
      this.fieldsColumnNames = this.dataFrame.columns.names();

    const sub = (stream: Observable<unknown>, action: Function) => {
      this.subs.push(DG.debounce(stream, 50).subscribe((_) => action()));
    };

    sub(this.dataFrame.selection.onChanged, () => this.render());
    sub(this.dataFrame.filter.onChanged, () => this.render());

    sub(this.dataFrame.onMetadataChanged, () => this.render());

    sub(this.dataFrame.onColumnsRemoved, () => {
      this.updatefieldsColumnNames();
      this.render();
    });

    sub(this.dataFrame.onCurrentRowChanged.pipe(filter((_) => this.showCurrentRow)),
      () => this.virtualView.refreshItem(this.currentRowPos!));

    sub(this.dataFrame.onMouseOverRowChanged.pipe(filter((_) => this.showMouseOverRow)),
      () => this.virtualView.refreshItem(this.mouseOverPos!));

    this.render();
  }

  updatefieldsColumnNames() {
    const newColumns = this.dataFrame.columns.names();
    let counter = this.fieldsColumnNames.length;
    for (let i = 0; i < counter; i++) {
      const colIdx = newColumns.indexOf(this.fieldsColumnNames[i]);
      if (colIdx === -1) {
        this.fieldsColumnNames.splice(i, 1);
        counter--;
        i--;
      }
    }
  }

  onPropertyChanged(property: DG.Property): void {
    if (property.name === 'showCurrentRow' || property.name === 'showMouseOverRow' ||
      property.name === 'showSelectedRows' && property.get(this) === true && this.showFixedRows)
      grok.shell.warning(`Cannot set ${property.name} to true since fixed rows are set`);
    this.render();
  }

  renderHeader() {
    ui.empty(this.columnHeadersDiv);
    const form = this.renderForm(0);
    form.classList.add('temp');
    document.body.appendChild(form);

    for (const name of this.fieldsColumnNames) {
      const formField = form.querySelector('[column="' + name + '"]');
      if (formField) {
        const columnLabel = ui.bind(this.dataFrame.col(name), ui.divText(name, 'd4-multi-form-column-name'));
        columnLabel.style.top = `${(formField as HTMLElement).offsetTop + 10}px`;
        this.columnHeadersDiv.appendChild(columnLabel);
      }
    }

    form.remove();
  }

  get currentRowPos() { return this.showCurrentRow ? 0 : null; }
  get mouseOverPos() { return this.showMouseOverRow ? (this.showCurrentRow ? 1 : 0) : null; }

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
        if (input) {
          input.input.setAttribute('column', name);
          input.value = this.dataFrame.get(name, row);

          if (this.colorCode) {
            const grid = ((this.view ?? grok.shell.tv) as DG.TableView).grid;
            const color = grid.cell(name, row).color;
            input.input.style.color = DG.Color.toHtml(DG.Color.getContrastColor(color));
            input.input.style.backgroundColor = DG.Color.toHtml(color);
          }

          return input.input;
        }
        return ui.div();
      }
      ), 'd4-multi-form-form');

    form.onclick = () => this.dataFrame.currentRowIdx = row;
    form.onmouseenter = () => this.dataFrame.mouseOverRowIdx = row;
    form.onmouseleave = () => this.dataFrame.mouseOverRowIdx = -1;
    return form;
  }

  render() {
    if (!this.showFixedRows) {
      if (this.showSelectedRows) {
        this.indexes = this.dataFrame.selection.trueCount > 0 ?
          Array.from(this.dataFrame.selection.getSelectedIndexes()) : [];
      }
    }

    ui.empty(this.columnHeadersDiv);
    if (this.dataFrame.currentRowIdx === -1)
      this.dataFrame.currentRowIdx = 0;

    this.renderHeader();

    this.virtualView.setData(
      this.indexes.length + (this.showCurrentRow ? 1 : 0) + (this.showMouseOverRow ? 1 : 0),
      (i: number) => this.renderForm(i));
  }
}
