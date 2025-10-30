import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {Observable, Subject} from 'rxjs';
import {filter} from 'rxjs/operators';
import '../../css/forms.css';

const COLS_LIMIT_EXCEEDED_WARNING = `Number of columns is more than 20. First 20 columns are shown`;
const BOOLEAN_INPUT_TOP_MARGIN = 15;

export class FormsViewer extends DG.JsViewer {
  get type(): string { return 'FormsViewer'; }
  moleculeSize: 'small' | 'normal' | 'large';
  fieldsColumnNames: string[];
  colorCode: boolean;
  showCurrentRow: boolean;
  showMouseOverRow: boolean;
  showSelectedRows: boolean;
  showFixedRows = false;
  indexes: number[];
  columnHeadersDiv: HTMLDivElement;
  virtualView: DG.VirtualView;
  columnLabelWidth: number = 0;
  currentRowIndicator = ui.div('', 'd4-multi-form-form-indicator d4-multi-form-form-indicator-current-row');
  mouseOverRowIndicator = ui.div('', 'd4-multi-form-form-indicator d4-multi-form-form-indicator-mouse-over-row');
  splitColLeft: HTMLElement;
  splitColRight: HTMLElement;
  inputClicked: Subject<string> = new Subject();

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

  protected getRendererSize(renderer: DG.GridCellRenderer): DG.Point {
    let width = renderer.defaultWidth;
    let height = renderer.defaultHeight;
    if (!width || !height)
      return this.getMoleculeSize();

    switch (this.moleculeSize) {
      case 'normal': return new DG.Point(width, height);
      case 'large': return new DG.Point(Math.floor(width * 1.5), Math.floor(height * 1.5));
      case 'small':
      default: return new DG.Point(Math.floor(width * 0.66), Math.floor(height * 0.66));
    }
  }

  protected getMoleculeSize(): DG.Point {
    switch (this.moleculeSize) {
      case 'normal': return new DG.Point(200, 100);
      case 'large': return new DG.Point(300, 150);
      case 'small':
      default: return new DG.Point(120, 60);
    }
  }

  constructor() {
    super();

    this.helpUrl = 'https://datagrok.ai/help/visualize/viewers/forms';

    // properties
    this.fieldsColumnNames = this.addProperty('fieldsColumnNames', DG.TYPE.COLUMN_LIST);
    this.colorCode = this.bool('colorCode', true);
    this.showCurrentRow = this.bool('showCurrentRow', true);
    this.showMouseOverRow = this.bool('showMouseOverRow', true);
    this.showSelectedRows = this.bool('showSelectedRows', true);
    this.moleculeSize = this.string('moleculeSize', 'small', {choices: ['small', 'normal', 'large']}) as 'small' | 'normal' | 'large';

    //fields
    this.indexes = [];

    // init
    this.root.classList.add('d4-multi-form');
    this.columnHeadersDiv = ui.div([], 'd4-multi-form-header');
    this.virtualView = ui.virtualView(0, (i: number) => this.renderForm(i), false, 1);
    const columnHeadersBox = ui.div(this.columnHeadersDiv);
    const formWithHeaderDiv = ui.splitH([columnHeadersBox, this.virtualView.root], null, true);
    this.root.appendChild(formWithHeaderDiv);

    this.splitColLeft = formWithHeaderDiv.firstElementChild as HTMLElement;
    this.splitColRight = formWithHeaderDiv.lastElementChild as HTMLElement;

    ui.tooltip.bind(this.currentRowIndicator, 'Current row');
    ui.tooltip.bind(this.mouseOverRowIndicator, 'Mouse over row');

    ui.tools.waitForElementInDom(this.root).then((_) => {
      this.columnHeadersDiv.style.cssText = `
        oveflow:hidden!important;
        min-width: 150px;
        flex-shrink: 0;
      `;
      this.splitColLeft.style.cssText = `
        overflow: hidden!important;
      `;
      columnHeadersBox.style.cssText = `
        box-sizing: content-box;
        width:100%;
        position:relative;
        display:flex;
        padding-right:17px;
        overflow: scroll!important;
      `;

      columnHeadersBox.addEventListener('scroll', (e: Event) => {
        this.virtualView.root.scrollTop = columnHeadersBox.scrollTop;
      });

      this.virtualView.root.addEventListener('scroll', (e: Event) => {
        columnHeadersBox.scrollTop = this.virtualView.root.scrollTop;
      });
    });

    ui.tools.waitForElementInDom(this.virtualView.root).then((_) => {
      this.fitHeaderToLabelWidth();
    });
  }

  fitHeaderToLabelWidth(width?: number) {
    const w = width ?? this.columnLabelWidth;
    const rootWidth = this.root.getBoundingClientRect().width;
    this.splitColLeft.style.width = `${w + 30}px`;
    this.splitColRight.style.width = `${rootWidth - w - 30 - 4}px`;
  }

  getColumnWidth(el:HTMLElement) {
    return parseInt(el.style.width);
  }

  onTableAttached() {
    if (this.fieldsColumnNames === null)
      this.setFieldsColumnNames(this.dataFrame.columns.names());

    const sub = (stream: Observable<unknown>, action: Function) => {
      this.subs.push(DG.debounce(stream, 50).subscribe((_) => action()));
    };

    sub(this.dataFrame.selection.onChanged, () => this.render());
    sub(this.dataFrame.filter.onChanged, () => this.render());

    sub(this.dataFrame.onMetadataChanged, () => this.render());

    sub(this.dataFrame.onColumnsRemoved, () => {
      this.updateFieldsColumnNames();
      this.render();
    });

    sub(this.dataFrame.onCurrentRowChanged.pipe(filter((_) => this.showCurrentRow)),
      () => this.virtualView.refreshItem(this.currentRowPos!));

    sub(this.dataFrame.onMouseOverRowChanged.pipe(filter((_) => this.showMouseOverRow)),
      () => this.virtualView.refreshItem(this.mouseOverPos!));

    this.render();
  }

  setFieldsColumnNames(dfColumns: string[]) {
    this.fieldsColumnNames = dfColumns.length > 20 ? dfColumns.slice(0, 20) : dfColumns;
  }

  updateFieldsColumnNames() {
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
    if ((property.name === 'showCurrentRow' || property.name === 'showMouseOverRow' ||
      property.name === 'showSelectedRows' && property.get(this) === true) && this.showFixedRows)
      grok.shell.warning(`Cannot set ${property.name} to true since fixed rows are set`);
    this.render();
  }

  renderHeader() {
    ui.empty(this.columnHeadersDiv);
    const form = this.renderForm(0, true);
    form.classList.add('temp');
    document.body.appendChild(form);

    for (const name of this.fieldsColumnNames) {
      const formField = form.querySelector('[column="' + name + '"]') as HTMLInputElement;
      if (formField) {
        const columnLabel = ui.bind(this.dataFrame.col(name), ui.divText(name));
        const closeIcon = ui.iconFA('times', () => {
          this.fieldsColumnNames = this.fieldsColumnNames.filter((value) => value !== name);
          this.render();
        });
        const columnLabelContainer = ui.div([columnLabel, closeIcon], 'd4-multi-form-column-name d4-flex-row');
        const offsetTop = formField.type === 'checkbox' ?
          formField.offsetTop - BOOLEAN_INPUT_TOP_MARGIN : formField.offsetTop;
        columnLabelContainer.style.top = `${offsetTop + 10}px`;
        this.columnHeadersDiv.appendChild(columnLabelContainer);
        columnLabel.onclick = (e: MouseEvent) => {
          this.dataFrame.currentCol = this.dataFrame.col(name)!;
        };
        if (columnLabel.getBoundingClientRect().width > this.columnLabelWidth)
          this.columnLabelWidth = columnLabel.getBoundingClientRect().width;
        columnLabelContainer.onmouseenter = (e: MouseEvent) => closeIcon.style.visibility = 'visible';
        columnLabelContainer.onmouseleave = (e: MouseEvent) => closeIcon.style.visibility = 'hidden';
      }
    }

    form.remove();
  }

  get currentRowPos() { return this.showCurrentRow ? 0 : null; }
  get mouseOverPos() { return this.showMouseOverRow ? (this.showCurrentRow ? 1 : 0) : null; }

  getGrid() {
    if (this.view && (this.view as DG.TableView).grid)
      return (this.view as DG.TableView).grid;
    if (this.dataFrame && grok.shell.tv && grok.shell.tv.dataFrame === this.dataFrame && grok.shell.tv.grid)
      return grok.shell.tv.grid;
    const tv = Array.from(grok.shell.tableViews).find((tv) => tv.dataFrame === this.dataFrame);
    if (tv && tv.grid)
      return tv.grid;
    return null;
  }

  renderForm(row: number, header?: boolean) {
    const savedIdx = row;
    if (header)
      row = 0;
    else {
      if (this.showCurrentRow && row === this.currentRowPos)
        row = this.dataFrame.currentRowIdx;
      else if (this.showMouseOverRow && row === this.mouseOverPos)
        row = this.dataFrame.mouseOverRowIdx;
      else
        row = this.indexes[row - (this.showCurrentRow ? 1 : 0) - (this.showMouseOverRow ? 1 : 0)];
    }

    const form = ui.divV(
      this.fieldsColumnNames.map((name) => {
        let resDiv: HTMLElement = ui.div();
        if (row === -1)
          return resDiv;

        try {
          const grid = this.getGrid();
          // for molecules, use the gridCol renderer instead of inputBase
          const col = this.dataFrame.col(name)!;
          if (col.semType && grid?.col(name)?.renderer) {
            const renderer = grid.col(name)!.renderer!;
            const rendererSize = this.getRendererSize(renderer);
            const gridCell = DG.GridCell.fromColumnRow(grid, name, grid.tableRowToGrid(row));
            const canvas = ui.canvas(rendererSize.x, rendererSize.y);
            canvas.width = rendererSize.x * window.devicePixelRatio;
            canvas.height = rendererSize.y * window.devicePixelRatio;
            const ctx = canvas.getContext('2d')!;
            // ctx.fillRect(0, 0, canvas.width, canvas.height);
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
            
            canvas.setAttribute('column', name);

            gridCell.render({context: ctx, bounds: new DG.Rect(0, 0, rendererSize.x, rendererSize.y)});
            resDiv = canvas;
            resDiv.style.background = 'white'; // fixes before changes of viertual view not removing old items
            resDiv.onclick = (e: MouseEvent) => {
              this.dataFrame.currentCell = this.dataFrame.cell(row, name);
              this.inputClicked.next(name);
            };
            ui.tooltip.bind(resDiv, name);

          } else {
            const input = DG.InputBase.forColumn(col);
            if (input) {
              if (this.dataFrame.col(name)!.semType === DG.SEMTYPE.MOLECULE)
                input.input.classList.add(`d4-multi-form-molecule-input-${this.moleculeSize}`);
              input.input.setAttribute('column', name);
              input.value = this.dataFrame.get(name, row);
              input.readOnly = true;

              if (this.colorCode) {
                if (grid) {
                  const color = this.dataFrame.col(name)?.meta.colors.getColor(row);
                  if (color) {
                    if (grid.col(name)?.isTextColorCoded)
                      input.input.setAttribute('style', `color:${DG.Color.toHtml(color)}!important;`);
                    else {
                      input.input.setAttribute('style',
                        `color:${DG.Color.toHtml(DG.Color.getContrastColor(color))}!important;`);
                      input.input.style.backgroundColor = DG.Color.toHtml(color);
                    }
                  }
                }
              }
              input.input.onclick = (e: MouseEvent) => {
                this.dataFrame.currentCell = this.dataFrame.cell(row, name);
                this.inputClicked.next(name);
              };
              ui.tooltip.bind(input.input, name);

              resDiv = input.input;
            }
          }
        } catch (e) {
          console.error(e);
         }
        return resDiv;
      }
      ), 'd4-multi-form-form');

    if (this.showMouseOverRow && savedIdx == this.mouseOverPos)
      form.append(this.mouseOverRowIndicator);
    if (this.showCurrentRow && savedIdx === this.currentRowPos)
      form.append(this.currentRowIndicator);

    form.onclick = (event: MouseEvent) => {
      if (event.ctrlKey && event.shiftKey) {
        for (let i = 0; i <= row; i++)
          this.dataFrame.selection.set(i, false);
      } else {
        if (event.ctrlKey) {
          const currentSelection = this.dataFrame.selection.get(row);
          this.dataFrame.selection.set(row, !currentSelection);
        } else if (event.shiftKey) {
          for (let i = 0; i <= this.dataFrame.rowCount; i++)
            this.dataFrame.selection.set(i, i <= row);
        } else {
          if (!this.showCurrentRow || savedIdx !== this.currentRowPos)
            this.dataFrame.currentRowIdx = row;
        }
      }
    };
    if (!this.showMouseOverRow || savedIdx !== this.mouseOverPos) {
      form.onmouseenter = () => this.dataFrame.mouseOverRowIdx = row;
      form.onmouseleave = () => this.dataFrame.mouseOverRowIdx = -1;
    }
    return form;
  }

  render() {
    if (!this.showFixedRows) {
      if (this.showSelectedRows) {
        this.indexes = this.dataFrame.selection.trueCount > 0 ?
          Array.from(this.dataFrame.selection.getSelectedIndexes()) : [];
      }
      else
        this.indexes = [];
    }

    ui.empty(this.columnHeadersDiv);

    this.renderHeader();

    this.virtualView.setData(
      this.indexes.length + (this.showCurrentRow ? 1 : 0) + (this.showMouseOverRow ? 1 : 0),
      (i: number) => this.renderForm(i));
  }
}
