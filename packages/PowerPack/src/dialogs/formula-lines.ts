import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Observable } from 'rxjs';

// Formula Line types:
const enum ITEM_TYPE {
  PARALLEL_LINE = 'Parallel line',
  UNIVERSAL_LINE = 'Universal line',
  BAND = 'Band'
}

const FOR_TABLE = 'For table';

// Options for Formula Line types:
const ITEM_OPTIONS = {

  // Editor options for Parallel Line:
  [ITEM_TYPE.PARALLEL_LINE]: [
    { name: 'title', type: 'string', caption: 'Title' },
    { name: 'description', type: 'string', caption: 'Description' },
    { name: 'column', type: 'column', caption: 'Column' },
    { name: 'value', type: 'double', caption: 'Value' },
    { name: 'width', type: 'int', caption: 'Width' },
    { name: 'visible', type: 'bool', caption: 'Show' },
    { name: 'opacity', type: 'int', caption: 'Opacity', editor: 'slider', min: 0, max: 100 },
    { name: 'zIndex', type: 'int', caption: 'zIndex' },
    { name: 'color', type: 'string', caption: 'Color' },
    { name: 'style', type: 'string', caption: 'Style', choices: ['solid', 'dotted', 'dashed', 'longdash', 'dotdash'] },
    { name: 'min', type: 'int', caption: 'Min' },
    { name: 'max', type: 'int', caption: 'Max' }
  ],

  // Editor options for Universal Line:
  [ITEM_TYPE.UNIVERSAL_LINE]: [
    { name: 'title', type: 'string', caption: 'Title' },
    { name: 'description', type: 'string', caption: 'Description' },
    { name: 'formula', type: 'string', editor: 'textarea', caption: 'Formula' },
    { name: 'width', type: 'int', caption: 'Width' },
    { name: 'visible', type: 'bool', caption: 'Show' },
    { name: 'opacity', type: 'int', caption: 'Opacity', editor: 'slider', min: 0, max: 100 },
    { name: 'zIndex', type: 'int', caption: 'zIndex' },
    { name: 'color', type: 'string', caption: 'Color' },
    { name: 'style', type: 'string', caption: 'Style', choices: ['solid', 'dotted', 'dashed', 'longdash', 'dotdash'] },
    { name: 'min', type: 'int', caption: 'Min' },
    { name: 'max', type: 'int', caption: 'Max' }
  ],

  // Editor options for Band:
  [ITEM_TYPE.BAND]: [
    { name: 'title', type: 'string', caption: 'Title' },
    { name: 'description', type: 'string', caption: 'Description' },
    { name: 'column', type: 'string', caption: 'Column 1' },
    { name: 'formula', type: 'string', editor: 'textarea', caption: 'Formula' },
    { name: 'column2', type: 'string', caption: 'Column 2' },
    { name: 'visible', type: 'bool', caption: 'Show' },
    { name: 'opacity', type: 'int', caption: 'Opacity', editor: 'slider', min: 0, max: 100 },
    { name: 'zIndex', type: 'int', caption: 'zIndex' },
    { name: 'color', type: 'string', caption: 'Color' },
    { name: 'min', type: 'int', caption: 'Min' },
    { name: 'max', type: 'int', caption: 'Max' }
  ],

  // Options for displaying Formula Lines in the Table:
  [FOR_TABLE]: [
    { name: 'title', type: 'string', caption: 'Title' },
    { name: 'formula', type: 'string', editor: 'textarea', caption: 'Formula' },
    { name: 'visible', type: 'bool', caption: 'Show' },
    { name: 'new', type: 'string', caption: 'new' }
  ]
};

// Creates properties from options:
function propsFromOptions(options: any[]): DG.Property[] {
  return options.map((p: any) => DG.Property.fromOptions(p));
}

/**
 * Formula Lines Host Helper.
 * Reads, storages and saves Formula Lines to the host (DataFrame or Viewer).
 */
class Host {
  _src: DG.DataFrame | DG.Viewer;
  items: DG.FormulaLine[];

  constructor(src: DG.DataFrame | DG.Viewer) {
    if (!src)
      throw 'Host table/viewer not found.';

    this._src = src;
    this.items = src.meta.getFormulaLines();
  }

  save(): void {
    this._src.meta.removeFormulaLines();
    this._src.meta.addFormulaLines(this.items);
  }
}

/**
 * Table Helper for displaying and navigating Formula Lines list.
 */
class Table {
  _dataFrame?: DG.DataFrame;
  _grid: DG.Grid;
  _items: DG.FormulaLine[];
  _props: DG.Property[] = propsFromOptions(ITEM_OPTIONS[FOR_TABLE]);

  get _onCurrentItemChanged(): Observable<any> { return this._dataFrame!.onCurrentRowChanged; }
  get _onValuesChanged(): Observable<any> { return this._dataFrame!.onValuesChanged; }

  get _currentItemIdx(): number { return this._dataFrame!.currentRowIdx;}
  set _currentItemIdx(rowIdx: number) { this._dataFrame!.currentRowIdx = rowIdx; }

  currentItem?: DG.FormulaLine;
  set height(h: number) { this._grid.root.style.height = `${h}px`; }
  get root(): HTMLElement { return this._grid.root; }

  constructor(items: DG.FormulaLine[], onChangedAction: Function) {
    this._items = items;
    this._dataFrame = DG.DataFrame.fromObjects(items);
    this._grid = DG.Grid.create(this._dataFrame!);

    this._grid.setOptions({
      showRowHeader: false,
      showCellTooltip: false,
      showColumnTooltip: false,
      showCurrentCellOutline: false,
      showContextMenu: false
    });

    this._grid.columns.setVisible(['title', 'formula', 'visible']);

    this._grid.columns.byName('title')!.width = 120;
    this._grid.columns.byName('formula')!.width = 220;
    this._grid.columns.byName('visible')!.width = 40;

    (this._dataFrame!.columns as DG.ColumnList).addNew('New', DG.TYPE.STRING);
    let ctrlCol = this._grid.columns.byName('New')!;
    ctrlCol.width = 35;
    ctrlCol.cellType = 'html';

    // Shows "Delete" and "Add new" buttons:
    this._grid.onCellPrepare((cell) => {
      if (cell.isColHeader) {
        if (cell.gridColumn.name == 'title')
          cell.customText = 'Title';
        else if (cell.gridColumn.name == 'formula')
          cell.customText = 'Formula';
        else if (cell.gridColumn.name == 'visible')
          cell.customText = 'Show';
      }

      if (cell.gridColumn.name == 'New')
        cell.style.element = cell.isTableCell ? deleteBtn(cell.gridRow) : addBtn();
    });

    // Set current Formula Line and execute callback "onCurrentItemChangedAction":
    this._onCurrentItemChanged.subscribe((_) => {
      if (this._currentItemIdx < 0)
        return;
      let item = items[this._currentItemIdx];
      this.currentItem = item;
      onChangedAction(this._currentItemIdx);
    });

    this._onValuesChanged.subscribe((_) => {
      let item = this._items[this._currentItemIdx];
      item.title = this._dataFrame!.getCol('title').get(this._currentItemIdx);
      item.formula = this._dataFrame!.getCol('formula').get(this._currentItemIdx);
      item.visible = this._dataFrame!.getCol('visible').get(this._currentItemIdx);
      onChangedAction(this._currentItemIdx);
    });

    // Set first Formula Line as current:
    if (items.length > 0)
      this._currentItemIdx = 0;

    // Link to the Table instance to use in callbacks:
    let $this = this;

    // Creates "Add new" button:
    function addBtn(): HTMLElement {
      let btn = ui.button(ui.iconFA('plus'), () => {

      });
      return btn;
    }

    // Creates "Delete" button:
    function deleteBtn(itemIdx: number): HTMLElement {
      let btn = ui.button(ui.iconFA('trash-alt'), () => {
        items.splice(itemIdx, 1);
        $this._dataFrame!.rows.removeAt(itemIdx);
        if ($this._currentItemIdx > itemIdx)
          $this._currentItemIdx--;
        if ($this._currentItemIdx >= 0) {
          let item = items[$this._currentItemIdx];
          $this.currentItem = item;
          onChangedAction($this._currentItemIdx);
        }
      });
      btn.style.textAlign = 'center';
      btn.style.height = '20px';
      return btn;
    }
  }

  update(itemIdx: number): void {
    let item = this._items[itemIdx];
    this._dataFrame!.rows.setValues(itemIdx, [item.title, item.formula]);
    this._grid.invalidate();
  }
}

/**
 * Preview Helper for Formula Lines.
 * Scatter Plot viewer by default.
 */
class Preview {
  scatterPlot: DG.ScatterPlotViewer;
  dataFrame?: DG.DataFrame;
  _src: DG.DataFrame | DG.Viewer;
  _items: DG.FormulaLine[];

  set height(h: number) { this.scatterPlot.root.style.height = `${h}px`; }
  get root(): HTMLElement { return this.scatterPlot.root; }

  // Sets the corresponding axes:
  _setAxes(item: DG.FormulaLine): void {
    let [itemY, itemX] = this.scatterPlot.meta.getFormulaInfo(item, this.dataFrame!);
    let [previewY, previewX] = [itemY, itemX];

    // If the source is a Scatter Plot, then we try to set similar axes:
    if (this._src instanceof DG.ScatterPlotViewer) {
      let [srcY, srcX] = [this._src.props.yColumnName, this._src.props.xColumnName];
      [previewY, previewX] = [previewY ?? srcY, previewX ?? srcX];

      if (previewX == srcY || previewY == srcX)
        [previewY, previewX] = [previewX, previewY];

      if (previewX == previewY)
        previewY = srcY;
    }

    if (previewY)
      if (this.dataFrame!.getCol(previewY))
        this.scatterPlot.setOptions({y: previewY});

    if (previewX)
      if (this.dataFrame!.getCol(previewX))
        this.scatterPlot.setOptions({x: previewX});
  }

  constructor(items: DG.FormulaLine[], src: DG.DataFrame | DG.Viewer) {
    this._items = items;
    // Extract data for Formula Lines:
    this.dataFrame =
          src instanceof DG.DataFrame ? src
        : src instanceof DG.Viewer ? src.dataFrame
        : undefined;

    if (!this.dataFrame)
      throw 'Host is not DataFrame or Viewer.';

    this._src = src;

    this.scatterPlot = DG.Viewer.scatterPlot(this.dataFrame, {
      showDataframeFormulaLines: false,
      showViewerFormulaLines: true,
      showSizeSelector: false,
      showColorSelector: false,
      showContextMenu: false,
      axesFollowFilter: false,
      showMinMaxTickmarks: false,
      showMouseOverPoint: false,
      showCurrentPoint: false,
      zoomAndFilter: 'no action',
      axisFont: '11px Arial',
      legendVisibility: 'Never',
      xAxisHeight: 25
    });
  }

  update(itemIdx: number): boolean {
    if (this.scatterPlot == undefined)
      return false;

    // Duplicate the original item to display it even if it's hidden:
    let item = this._items[itemIdx];
    let previewItem = Object.assign({}, item);
    previewItem.visible = true;

    // Show the item:
    this.scatterPlot.meta.removeFormulaLines();
    try {
      this.scatterPlot.meta.addFormulaItem(previewItem);
      this._setAxes(previewItem);
      return true;
    } catch {
      this.scatterPlot.meta.removeFormulaLines();
      return false;
    }
  }
}

/**
 * Editor Helper for Formula Lines (form with corresponding inputs).
 */
class Editor {
  _form: HTMLElement;
  _preview: Preview;
  _items: DG.FormulaLine[];
  _onChangedAction: Function;

  // Set of properties for different Formula Line types:
  _props = {
    [ITEM_TYPE.PARALLEL_LINE]: propsFromOptions(ITEM_OPTIONS[ITEM_TYPE.PARALLEL_LINE]),
    [ITEM_TYPE.UNIVERSAL_LINE]: propsFromOptions(ITEM_OPTIONS[ITEM_TYPE.UNIVERSAL_LINE]),
    [ITEM_TYPE.BAND]: propsFromOptions(ITEM_OPTIONS[ITEM_TYPE.BAND])
  }

  // Returns corresponding properties for given Formula Line:
  _getItemProps(item?: DG.FormulaLine): DG.Property[] {
    if (item?.type == 'band')
      return this._props[ITEM_TYPE.BAND];
    return this._props[ITEM_TYPE.UNIVERSAL_LINE];
  }

  get root(): HTMLElement { return this._form; }

  constructor(items: DG.FormulaLine[], preview: Preview, onChangedAction: Function) {
    this._form = ui.form();
    this._items = items;
    this._preview = preview;
    this._onChangedAction = onChangedAction;
  }

  // Creates and fills editor for given Formula Line:
  update(itemIdx: number): void {
    if (this._form == undefined)
      return;

    let newForm = this._createForm(itemIdx);

    this._form.replaceWith(newForm);
    this._form = newForm;
  }

  _inputFormula(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibFormula = ui.textInput('', item.formula ?? '', (value: string) => {
      item.formula = value;
      let resultOk = this._onChangedAction(itemIdx);
      if (resultOk)
        (ibFormula.input as HTMLInputElement).classList.remove('d4-forced-invalid')
      else
        (ibFormula.input as HTMLInputElement).classList.add('d4-forced-invalid')
    });

    let elFormula = ibFormula.input as HTMLInputElement;
        elFormula.placeholder = 'Formula';
        elFormula.style.width = '360px';
        elFormula.style.height = '137px';
        elFormula.style.marginRight = '-6px';

    ui.tools.initFormulaAccelerators(ibFormula, this._preview.dataFrame!);

    return ibFormula.root;
  }

  _inputColor(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibColor = ui.colorInput('Color', item.color ?? '#000000', (value: string) => {
      item.color = value;
      this._onChangedAction(itemIdx);
    });

    let elColor = ibColor.input as HTMLInputElement;
        elColor.style.maxWidth = 'none';
        elColor.style.width = '204px';

    return ui.divH([ibColor.root]);
  }

  _inputOpacity(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibOpacity = ui.element('input');
        ibOpacity.type = 'range';
        ibOpacity.min = 0;
        ibOpacity.max = 100;
        ibOpacity.value = item.opacity ?? 100;
        ibOpacity.style.width = '204px';
        ibOpacity.style.marginLeft = '0px';
        ibOpacity.style.marginTop = '6px';
        ibOpacity.addEventListener('input', () => {
          item.opacity = parseInt(ibOpacity.value);
          this._onChangedAction(itemIdx);
        });

    let label = ui.label('Opacity', 'ui-label ui-input-label');

    return ui.divH([ui.div([label, ibOpacity], 'ui-input-root')]);
  }

  _inputStyle(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibStyle = ui.choiceInput('Style', item.style ?? 'solid', ['solid', 'dotted', 'dashed', 'longdash', 'dotdash'], (value: string) => {
      item.style = value;
      this._onChangedAction(itemIdx);
    });

    let elStyle = ibStyle.input as HTMLInputElement;
        elStyle.style.width = '135px';

    let ibWidth = ui.intInput('', item.width ?? 1, (value: number) => {
      item.width = value;
      this._onChangedAction(itemIdx);
    });

    let elWidth = ibWidth.input as HTMLInputElement;
        elWidth.placeholder = '1';
        elWidth.style.width = '61px';
        elWidth.style.paddingRight = '24px';

    let unit = ui.divText('px', {style: {marginTop: '10px', marginLeft: '-24px', zIndex: '1000'} });

    return ui.divH([ibStyle.root, ibWidth.root, unit]);
  }

  _inputRange(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibMin = ui.stringInput('Range', `${item.min ?? ''}`, (value: string) => {
      item.min = value.length == 0 ? undefined : Number(value);
      this._onChangedAction(itemIdx);
    });

    let elMin = ibMin.input as HTMLInputElement;
        elMin.placeholder = 'min';
        elMin.style.width = '98px';

    let ibMax = ui.stringInput('', `${item.max ?? ''}`, (value: string) => {
      item.max = value.length == 0 ? undefined : Number(value);
      this._onChangedAction(itemIdx);
    });

    let elMax = ibMax.input as HTMLInputElement;
        elMax.placeholder = 'max';
        elMax.style.width = '98px';

    return ui.divH([ibMin.root, ibMax.root]);
  }

  _inputArrange(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibArrange = ui.choiceInput('Arrange', item.zIndex && item.zIndex > 0 ? 'above markers' : 'below markers', ['above markers', 'below markers'], (value: string) => {
      item.zIndex = value == 'above markers' ? 100 : -100;
      this._onChangedAction(itemIdx);
    });

    let elArrange = ibArrange.input as HTMLInputElement;
        elArrange.style.maxWidth = 'none';
        elArrange.style.width = '204px';

    return ui.divH([ibArrange.root]);
  }

  _inputTitle(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibTitle = ui.stringInput('Title', item.title ?? '', (value: string) => {
      item.title = value;
      this._onChangedAction(itemIdx);
    });

    let elTitle = ibTitle.input as HTMLInputElement;
        elTitle.style.maxWidth = 'none';
        elTitle.style.width = '204px';

    return ui.divH([ibTitle.root]);
  }

  _inputDescription(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibDescription = ui.textInput('Description', item.description ?? '', (value: string) => {
      item!.description = value;
      this._onChangedAction(itemIdx);
    });

    let elDescription = ibDescription.input as HTMLInputElement;
        elDescription.style.width = '194px';
        elDescription.style.height = '40px';
        elDescription.style.paddingLeft = '6px';
        elDescription.style.marginRight = '-8px';
        elDescription.style.fontFamily = 'inherit';
        elDescription.style.fontSize = 'inherit';

    return ibDescription.root;
  }

  _inputColumn2(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];

    let ibColumn2 = ui.columnInput('Adjacent column', this._preview.dataFrame!, item.column2 ? this._preview.dataFrame!.col(item.column2) : null, (value: DG.Column) => {
      item!.column2 = value.name;
      this._onChangedAction(itemIdx);
    });

    let elColumn2 = ibColumn2.input as HTMLInputElement;
        elColumn2.style.maxWidth = 'none';
        elColumn2.style.width = '204px';

    return ui.divH([ibColumn2.root], {style: {marginLeft: '-8px'}});
  }

  _inputConstant(itemIdx: number, colName: string, value: string): HTMLElement {
    let item = this._items[itemIdx];

    let ibColumn = ui.columnInput('Column', this._preview.dataFrame!, colName ? this._preview.dataFrame!.col(colName) : null, (value: DG.Column) => {
      item!.formula = '${' + value + '} = ' + ibValue.value;
      this._onChangedAction(itemIdx);
    });

    let elColumn = ibColumn.input as HTMLInputElement;
        elColumn.style.maxWidth = 'none';
        elColumn.style.width = '204px';

    let ibValue = ui.stringInput('Value', value, (value: string) => {
      item!.formula = '${' + ibColumn.value + '} = ' + value;
      this._onChangedAction(itemIdx);
    });
    ibValue.nullable = false;

    let elValue = ibValue.input as HTMLInputElement;
        elValue.style.maxWidth = 'none';
        elValue.style.width = '204px';

    return ui.div([ibColumn.root, ibValue.root]);
  }

  _createForm(itemIdx: number): HTMLElement {
    let item = this._items[itemIdx];
    let type = 'Band';
    let itemY = '';
    let itemX = '';
    let expression = '';

    if (item.type != 'band') {
      [itemY, itemX, expression] = this._preview.scatterPlot.meta.getFormulaInfo(item, this._preview.dataFrame!);
      type = itemX ? 'Line' : 'Constant Line';
    }

    let mainPane = ui.div([], 'ui-form');
        if (type == 'Constant Line')
          mainPane.append(this._inputConstant(itemIdx, itemY, expression));
        else
          mainPane.append(this._inputFormula(itemIdx));
        if (type == 'Band')
          mainPane.append(this._inputColumn2(itemIdx));

    let formatPane = ui.div([], 'ui-form');
        formatPane.style.marginLeft = '-20px';
        formatPane.append(this._inputColor(itemIdx));
        formatPane.append(this._inputOpacity(itemIdx));
        if (type == 'Line')
          formatPane.append(this._inputStyle(itemIdx));
        formatPane.append(this._inputRange(itemIdx));
        formatPane.append(this._inputArrange(itemIdx));

    let tooltipPane = ui.div([], 'ui-form');
        tooltipPane.style.marginLeft = '-20px';
        tooltipPane.append(this._inputTitle(itemIdx));
        tooltipPane.append(this._inputDescription(itemIdx));

    let combinedPanels = ui.accordion();
        combinedPanels.addPane(type, () => mainPane, true);
        combinedPanels.addPane('Format', () => formatPane, true);
        combinedPanels.addPane('Tooltip', () => tooltipPane, true);

    return ui.div([combinedPanels.root]);
  }
}

/**
 * A Dialog window with Formula Lines list, preview and editor.
 */
export class FormulaLinesDialog {
  title: string = 'Formula Lines';
  helpUrl: string = '/help/develop/how-to/show-formula-lines.md';

  // Helpers:
  host: Host;
  preview: Preview;
  editor: Editor;
  table: Table;
  dialog: DG.Dialog;

  /**
   * Initializes all parameters and opens the Dialog window.
   * @param {DG.DataFrame | DG.Viewer} src - the source where the Formula Lines will be read from.
   */
  constructor(src: DG.DataFrame | DG.Viewer) {
    this.host = new Host(src);
    this.preview = new Preview(this.host.items, src);
    this.editor = new Editor(this.host.items, this.preview, (itemIdx: number): boolean => {
      this.table.update(itemIdx);
      return this.preview.update(itemIdx);
    });
    this.table = new Table(this.host.items, (itemIdx: number): boolean => {
      this.editor.update(itemIdx);
      return this.preview.update(itemIdx);
    });
    this.dialog = ui.dialog({ title: this.title, helpUrl: this.helpUrl });

    // Create Dialog window layout:
    let layout = ui.div([
      ui.block([
        this.table.root,
        this.preview.root
      ], { style: { width: '55%', paddingRight: '20px' } }),
      ui.block([
        this.editor.root
      ], { style: { width: '45%' } })
    ]);

    this.table.height = 230;
    this.preview.height = 300;

    // Setup and open the Dialog window:
    this.dialog
      .add(layout)
      .onOK(this.onOKAction.bind(this))
      .show({width: 850, height: 650});
  }

  onOKAction() {
    this.host.save();
  }
}
