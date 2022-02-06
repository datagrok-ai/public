import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Formula Line types */
const enum ITEM_TYPE {
  LINE = 'line',
  BAND = 'band'
}

/** Formula Line captions */
const enum ITEM_CAPTION {
  LINE = 'Line',
  BAND = 'Band',
  CONST_LINE = 'Constant line',
  VERT_LINE = 'Line - Vertical',
  HORZ_LINE = 'Line - Horizontal',
  HORZ_BAND = 'Band - Horizontal',
  VERT_BAND = 'Band - Vertical'
}

/** Formula Line sources */
const enum ITEM_SOURCE {
  VIEWER = 'Viewer',
  DATAFRAME = 'DataFrame'
}

/**
 * Returns Formula Line type by its user-friendly [caption].
 * @param {ITEM_CAPTION} caption
 * @returns {ITEM_TYPE}
 */
function getItemTypeByCaption(caption: string): string {
  switch (caption) {
    case ITEM_CAPTION.LINE:
    case ITEM_CAPTION.HORZ_LINE:
    case ITEM_CAPTION.VERT_LINE: return ITEM_TYPE.LINE;

    case ITEM_CAPTION.BAND:
    case ITEM_CAPTION.HORZ_BAND:
    case ITEM_CAPTION.VERT_BAND: return ITEM_TYPE.BAND;

    default: throw 'Unknown item caption.';
  }
}

/**
 * Formula Lines Host Helper.
 * Reads, storages and saves Formula Lines to the host (DataFrame and Viewer).
 */
class Host {
  dframeItems?: DG.FormulaLine[];
  viewerItems?: DG.FormulaLine[];

  _dframeHelper?: DG.FormulaLinesHelper;
  _viewerHelper?: DG.FormulaLinesHelper;

  constructor(src: DG.DataFrame | DG.Viewer) {
    if (!src)
      throw 'Source table/viewer not found.';

    if (src instanceof DG.DataFrame) {
      this._dframeHelper = src.meta.formulaLines;
      this.dframeItems = this._dframeHelper.items;
    } else if (src instanceof DG.Viewer) {
      this._viewerHelper = src.meta.formulaLines;
      this.viewerItems = this._viewerHelper.items;
      if (src.dataFrame) {
        this._dframeHelper = src.dataFrame.meta.formulaLines;
        this.dframeItems = this._dframeHelper.items;
      } else
        throw 'Viewer not attached to table.';
    }
    else
      throw 'Unknown source of Formula Lines.';
  }

  save() {
    if (this._dframeHelper)
      this._dframeHelper.clear().addAll(this.dframeItems!);
    if (this._viewerHelper)
      this._viewerHelper.clear().addAll(this.viewerItems!);
  }
}

/**
 * Table Helper for displaying and navigating Formula Lines list.
 */
class Table {
  _grid: DG.Grid;
  _onChangedAction: Function
  items: DG.FormulaLine[];

  get _dataFrame(): DG.DataFrame { return this._grid.dataFrame!; }

  get _currentItemIdx(): number { return this._dataFrame.currentRowIdx;}
  set _currentItemIdx(rowIdx: number) { this._dataFrame.currentRowIdx = rowIdx; }

  /** Creates "Delete" button */
  _deleteBtn(itemIdx: number): HTMLElement {
    let btn = ui.button(ui.iconFA('trash-alt'), () => {
      this.items.splice(itemIdx, 1);
      this._dataFrame.rows.removeAt(itemIdx);
      if (this._currentItemIdx > itemIdx)
        this._currentItemIdx--;
      if (this._currentItemIdx >= 0) {
        this._onChangedAction(this._currentItemIdx);
      }
    });
    btn.style.textAlign = 'center';
    btn.style.height = '20px';
    return btn;
  }

  /** Used to prevent onValuesChanged event when the grid changes itself */
  _notify: boolean = true;

  get root(): HTMLElement { return this._grid.root; }

  constructor(items: DG.FormulaLine[], onChangedAction: Function) {
    this.items = items;
    this._onChangedAction = onChangedAction;

    let dataFrame = this.items.length > 0
      ? DG.DataFrame.fromObjects(this.items)!
      : DG.DataFrame.fromColumns([
          DG.Column.fromList('string', 'title', []),
          DG.Column.fromList('string', 'formula', []),
          DG.Column.fromList('bool', 'visible', [])
        ]);

    /** Column for "trash" buttons */
    const BTN_COL_NAME = 'deleteBtn';
    dataFrame.columns.addNewString(BTN_COL_NAME);

    this._grid = DG.Grid.create(dataFrame);
    this._grid.setOptions({
      showRowHeader: false,
      showCellTooltip: false,
      showColumnTooltip: false,
      showCurrentCellOutline: false,
      showContextMenu: false
    });

    this._grid.columns.setVisible(['title', 'formula', 'visible', BTN_COL_NAME]);
    this._grid.columns.setOrder(['title', 'formula', 'visible', BTN_COL_NAME]);

    this._grid.col('title')!.width = 120;
    this._grid.col('formula')!.width = 220;
    this._grid.col('visible')!.width = 40;

    let deleteBtnCol = this._grid.col(BTN_COL_NAME)!;
    deleteBtnCol.width = 35;
    deleteBtnCol.cellType = 'html';

    this._grid.onCellPrepare((cell) => {
      if (cell.isColHeader)
        switch (cell.gridColumn.name) {
          case 'title': cell.customText = 'Title'; break;
          case 'formula': cell.customText = 'Formula'; break;
          case 'visible': cell.customText = 'Show'; break;
          case BTN_COL_NAME: cell.style.textColor = 0; break;
        }
      else if (cell.isTableCell && cell.gridColumn.name == BTN_COL_NAME)
        cell.style.element = this._deleteBtn(cell.gridRow);
    });

    this._dataFrame.onCurrentRowChanged.subscribe((_) => {
      this._onChangedAction(this._currentItemIdx);
    });

    this._dataFrame.onValuesChanged.subscribe((_) => {
      if (!this._notify)
        return;
      let item = this.items[this._currentItemIdx];
      item.title = this._dataFrame.get('title', this._currentItemIdx);
      item.formula = this._dataFrame.get('formula', this._currentItemIdx);
      item.visible = this._dataFrame.get('visible', this._currentItemIdx);
      this._onChangedAction(this._currentItemIdx);
    });

    this.setFirstItemAsCurrent();
  }

  setFirstItemAsCurrent() {
    if (this.items.length > 0) {
      this._currentItemIdx = 0;
      this._onChangedAction(0);
    }
    else
      this._onChangedAction(-1);
  }

  update(itemIdx: number) {
    let item = this.items[itemIdx];
    this._notify = false;
      this._dataFrame.set('title', itemIdx, item.title);
      this._dataFrame.set('formula', itemIdx, item.formula);
    this._notify = true;
  }

  add(item: DG.FormulaLine) {
    this.items.unshift(item);
    this._dataFrame.rows.insertAt(0);
    this._notify = false;
      this._dataFrame.set('title', 0, item.title);
      this._dataFrame.set('formula', 0, item.formula);
      this._dataFrame.set('visible', 0, item.visible);
      this._dataFrame.set('deleteBtn', 0, '');
    this._notify = true;
    this._currentItemIdx = 0;
  }
}

interface AxisNames {
  y?: string,
  x?: string
}

interface AxisColumns {
  y: DG.Column,
  x: DG.Column
}

/**
 * Preview Helper for Formula Lines.
 * Scatter Plot viewer by default.
 */
class Preview {
  scatterPlot: DG.ScatterPlotViewer;
  dataFrame: DG.DataFrame;
  items: DG.FormulaLine[];

  /** Source Scatter Plot axes */
  _scrAxes?: AxisNames;

  set height(h: number) { this.scatterPlot.root.style.height = `${h}px`; }
  get root(): HTMLElement { return this.scatterPlot.root; }

  /** Returns the current columns pair of the preview Scatter Plot */
  get axisCols(): AxisColumns {
    return {
      y: this.dataFrame.getCol(this.scatterPlot.props.yColumnName),
      x: this.dataFrame.getCol(this.scatterPlot.props.xColumnName)
    }
  }

  /** Sets the current axes of the preview Scatter Plot by column names */
  set _axes(names: AxisNames) {
    if (names && names.y && this.dataFrame.getCol(names.y))
      this.scatterPlot.setOptions({y: names.y});
    if (names && names.x && this.dataFrame.getCol(names.x))
      this.scatterPlot.setOptions({x: names.x});
  }

  /**
   * Extracts the axes names from the formula. If possible, adjusts the axes
   * of the formula to the axes of the original scatterplot.
   */
  _getItemAxes(item: DG.FormulaLine): AxisNames {
    let itemMeta = DG.FormulaLinesHelper.getMeta(item);
    let [previewY, previewX] = [itemMeta.funcName, itemMeta.argName];

    /** If the source axes exist, then we try to set similar axes */
    if (this._scrAxes) {
      [previewY, previewX] = [previewY ?? this._scrAxes.y, previewX ?? this._scrAxes.x];

      if (previewX == this._scrAxes.y || previewY == this._scrAxes.x)
        [previewY, previewX] = [previewX, previewY];

      if (previewX == previewY)
        previewY = this._scrAxes.y;
    }

    return { y: previewY, x: previewX }
  }

  constructor(items: DG.FormulaLine[], src: DG.DataFrame | DG.Viewer, onContextMenu: Function) {
    this.items = items;

    if (src instanceof DG.DataFrame)
      this.dataFrame = src;
    else if (src instanceof DG.Viewer) {
      this.dataFrame = src.dataFrame!;
      this._scrAxes = { y: src.props.yColumnName, x: src.props.xColumnName };
    }
    else
      throw 'Host is not DataFrame or Viewer.';

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

    /**
     * Creates special context menu for preview Scatter Plot.
     * Before opening the menu, it calculates the world coordinates of the click point.
     */
    this.scatterPlot.root.addEventListener('contextmenu', (event: MouseEvent) => {
      event.preventDefault();
      let worldPoint = this.scatterPlot.screenToWorld(event.offsetX, event.offsetY);
      onContextMenu(worldPoint.y, worldPoint.x);
    });
  }

  /**
   * Shows a line with [itemIdx] index on the Scatter Plot.
   * Returns true if the rendering was successful, false otherwise.
   */
  update(itemIdx: number): boolean {
    /** If there are no lines, try to set the axes as in the original Scatter Plot. */
    if (itemIdx < 0 && this._scrAxes)
      this._axes = this._scrAxes;

    /** Duplicate the original item to display it even if it's hidden */
    let item = this.items[itemIdx];
    let previewItem = Object.assign({}, item);
    previewItem.visible = true;

    /** Trying to show the item */
    this.scatterPlot.meta.formulaLines.clear();
    try {
      this.scatterPlot.meta.formulaLines.add(previewItem);
      this._axes = this._getItemAxes(previewItem);
      return true;
    } catch {
      this.scatterPlot.meta.formulaLines.clear();
      return false;
    }
  }
}

/**
 * Editor Helper for Formula Lines (form with corresponding inputs).
 */
class Editor {
  _form: HTMLElement;
  _dataFrame: DG.DataFrame;
  _onChangedAction: Function;

  items: DG.FormulaLine[];
  get root(): HTMLElement { return this._form; }

  constructor(items: DG.FormulaLine[], dataFrame: DG.DataFrame, onChangedAction: Function) {
    this._form = ui.form();
    this.items = items;
    this._dataFrame = dataFrame;
    this._onChangedAction = onChangedAction;
  }

  /** Creates and fills editor for given Formula Line */
  update(itemIdx: number) {
    let newForm = this._createForm(itemIdx);
    this._form.replaceWith(newForm);
    this._form = newForm;
  }

  _createForm(itemIdx: number): HTMLElement {
    let item = itemIdx >= 0 ? this.items[itemIdx] : { type: ITEM_TYPE.BAND };
    let caption = ITEM_CAPTION.BAND;
    let [itemY, itemX, expression] = ['', '', ''];

    if (itemIdx >= 0 && item.type != ITEM_TYPE.BAND) {
      let itemMeta = DG.FormulaLinesHelper.getMeta(item);
      [itemY, itemX, expression] = [itemMeta.funcName!, itemMeta.argName!, itemMeta.expression!];
      caption = itemX ? ITEM_CAPTION.LINE : ITEM_CAPTION.CONST_LINE;
    }

    let mainPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px'}});
    let formatPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px'}});
    let tooltipPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px'}});

    if (itemIdx >= 0) {
      /** Preparing the "Main" panel */
      mainPane.append(caption == ITEM_CAPTION.CONST_LINE
        ? this._inputConstant(itemIdx, itemY, expression)
        : this._inputFormula(itemIdx));
      if (caption == ITEM_CAPTION.BAND)
        mainPane.append(this._inputColumn2(itemIdx));

      /** Preparing the "Format" panel */
      formatPane.append(this._inputColor(itemIdx));
      formatPane.append(this._inputOpacity(itemIdx));
      if (caption != ITEM_CAPTION.BAND)
        formatPane.append(this._inputStyle(itemIdx));
      formatPane.append(this._inputRange(itemIdx));
      formatPane.append(this._inputArrange(itemIdx));

      /** Preparing the "Tooltip" panel */
      tooltipPane.append(this._inputTitle(itemIdx));
      tooltipPane.append(this._inputDescription(itemIdx));
    }

    /** Creating the accordion */
    let combinedPanels = ui.accordion();
    combinedPanels.addPane(itemIdx >= 0 ? caption : ITEM_CAPTION.LINE, () => mainPane, true);
    combinedPanels.addPane('Format', () => formatPane, true);
    combinedPanels.addPane('Tooltip', () => tooltipPane, true);

    return ui.div([combinedPanels.root]);
  }

  /** Creates textarea for item formula */
  _inputFormula(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibFormula = ui.textInput('', item.formula ?? '', (value: string) => {
      item.formula = value;
      let resultOk = this._onChangedAction(itemIdx);
      (ibFormula.input as HTMLInputElement).classList.toggle('d4-forced-invalid', !resultOk)
    });

    let elFormula = ibFormula.input as HTMLInputElement;
    elFormula.placeholder = 'Formula';
    elFormula.setAttribute('style', 'width: 360px; height: 137px; margin-right: -6px;');

    ui.tools.initFormulaAccelerators(ibFormula, this._dataFrame);

    return ibFormula.root;
  }

  /** Creates color picker for item color */
  _inputColor(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibColor = ui.colorInput('Color', item.color ?? '#000000', (value: string) => {
      item.color = value;
      this._onChangedAction(itemIdx);
    });

    let elColor = ibColor.input as HTMLInputElement;
    elColor.placeholder = '#000000';
    elColor.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibColor.root]);
  }

  /** Creates range slider for item opacity */
  _inputOpacity(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let elOpacity = ui.element('input');
    elOpacity.type = 'range';
    elOpacity.min = 0;
    elOpacity.max = 100;
    elOpacity.value = item.opacity ?? 100;
    elOpacity.addEventListener('input', () => {
      item.opacity = parseInt(elOpacity.value);
      this._onChangedAction(itemIdx);
    });
    elOpacity.setAttribute('style', 'width: 204px; margin-top: 6px; margin-left: 0px;');

    let label = ui.label('Opacity', 'ui-label ui-input-label');

    return ui.divH([ui.div([label, elOpacity], 'ui-input-root')]);
  }

  /** Creates combobox for item line style and text input for item width */
  _inputStyle(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

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
    elWidth.setAttribute('style', 'width: 61px; padding-right: 24px;');

    let unit = ui.divText('px', {style: {marginTop: '10px', marginLeft: '-24px', zIndex: '1000'} });

    return ui.divH([ibStyle.root, ibWidth.root, unit]);
  }

  /** Creates text inputs for min-max values of item */
  _inputRange(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibMin = ui.stringInput('Range', `${item.min ?? ''}`, (value: string) => {
      item.min = value.length == 0 ? undefined : Number(value);
      this._onChangedAction(itemIdx);
    });

    let elMin = ibMin.input as HTMLInputElement;
    elMin.placeholder = 'min';
    elMin.setAttribute('style', 'width: 98px;');

    let ibMax = ui.stringInput('', `${item.max ?? ''}`, (value: string) => {
      item.max = value.length == 0 ? undefined : Number(value);
      this._onChangedAction(itemIdx);
    });

    let elMax = ibMax.input as HTMLInputElement;
    elMax.placeholder = 'max';
    elMax.setAttribute('style', 'width: 98px;');

    return ui.divH([ibMin.root, ibMax.root]);
  }

  /** Creates combobox for item position (z-index) */
  _inputArrange(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibArrange = ui.choiceInput('Arrange', item.zIndex && item.zIndex > 0 ? 'above markers' : 'below markers', ['above markers', 'below markers'], (value: string) => {
      item.zIndex = value == 'above markers' ? 100 : -100;
      this._onChangedAction(itemIdx);
    });

    let elArrange = ibArrange.input as HTMLInputElement;
    elArrange.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibArrange.root]);
  }

  /** Creates text input for item title */
  _inputTitle(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibTitle = ui.stringInput('Title', item.title ?? '', (value: string) => {
      item.title = value;
      this._onChangedAction(itemIdx);
    });

    let elTitle = ibTitle.input as HTMLInputElement;
    elTitle.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibTitle.root]);
  }

  /** Creates textarea for item description */
  _inputDescription(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibDescription = ui.textInput('Description', item.description ?? '', (value: string) => {
      item!.description = value;
      this._onChangedAction(itemIdx);
    });

    let elDescription = ibDescription.input as HTMLInputElement;
    elDescription.setAttribute('style', 'width: 194px; height: 40px; padding-left: 6px; margin-right: -8px; font-family: inherit; font-size: inherit;');

    return ibDescription.root;
  }

  /** Creates column input for band second column */
  _inputColumn2(itemIdx: number): HTMLElement {
    let item = this.items[itemIdx];

    let ibColumn2 = ui.columnInput('Adjacent column', this._dataFrame, item.column2 ? this._dataFrame.col(item.column2) : null, (value: DG.Column) => {
      item!.column2 = value.name;
      this._onChangedAction(itemIdx);
    });

    let elColumn2 = ibColumn2.input as HTMLInputElement;
    elColumn2.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibColumn2.root]);
  }

  /** Creates column input and text input for constant item */
  _inputConstant(itemIdx: number, colName: string, value: string): HTMLElement {
    let item = this.items[itemIdx];

    let ibColumn = ui.columnInput('Column', this._dataFrame, colName ? this._dataFrame.col(colName) : null, (value: DG.Column) => {
      item!.formula = '${' + value + '} = ' + ibValue.value;
      this._onChangedAction(itemIdx);
    });

    let elColumn = ibColumn.input as HTMLInputElement;
    elColumn.setAttribute('style', 'width: 204px; max-width: none; margin-right: -10px;');

    let ibValue = ui.stringInput('Value', value, (value: string) => {
      item!.formula = '${' + ibColumn.value + '} = ' + value;
      this._onChangedAction(itemIdx);
    });
    ibValue.nullable = false;

    let elValue = ibValue.input as HTMLInputElement;
    elValue.setAttribute('style', 'width: 204px; max-width: none; margin-right: -10px;');

    return ui.div([ibColumn.root, ibValue.root]);
  }
}

/**
 * Helper that implements the logic of creating a Formula Line item of a given type.
 */
class CreationControl {
  /** Opens a popup menu with predefined new Formula Line item types */
  popupMenu: Function;

  /** Used to create constant lines passing through the mouse click point on the Scatter Plot */
  _getCols: Function;

  /** Creates a button and binds an item creation menu to it */
  get button(): HTMLElement {
    let btn = ui.bigButton('Add new', this.popupMenu);
    return ui.div([btn], {style: {width: '100%', textAlign: 'right'}});
  }

  constructor(getCols: Function, onCreatedAction: Function) {
    this._getCols = getCols;

    this.popupMenu = (valY?: number, valX?: number) => {
      let onClickAction = (itemCaption: string) => {
        let cols: AxisColumns = this._getCols();
        let colY = cols.y;
        let colX = cols.x;

        /** Create blank item with a given type */
        let item: DG.FormulaLine = DG.FormulaLinesHelper.setDefaults({
          type: getItemTypeByCaption(itemCaption),
          title: 'New title'
        });

        /** Fill the item with the necessary data */
        switch (itemCaption) {
          case ITEM_CAPTION.LINE:
            item.formula = '${' + colY.name + '} = ${' + colX.name + '}';
            break;

          case ITEM_CAPTION.VERT_LINE:
            let vertVal = (valX ?? colX.stats.q2).toFixed(1);
            item.formula = '${' + colX.name + '} = ' + vertVal;
            break;

          case ITEM_CAPTION.HORZ_LINE:
            let horzVal = (valY ?? colY.stats.q2).toFixed(1);
            item.formula = '${' + colY.name + '} = ' + horzVal;
            break;

          case ITEM_CAPTION.VERT_BAND:
            let left = colX.stats.q1.toFixed(1);
            let right = colX.stats.q3.toFixed(1);
            item.formula = '${' + colX.name + '} in(' + left + ', ' + right + ')';
            item.column2 = colY.name;
            break;

          case ITEM_CAPTION.HORZ_BAND:
            let bottom = colY.stats.q1.toFixed(1);
            let top = colY.stats.q3.toFixed(1);
            item.formula = '${' + colY.name + '} in(' + bottom + ', ' + top + ')';
            item.column2 = colX.name;
            break;
        }

        /** Used to update the Table, Preview and Editor states */
        onCreatedAction(item);
      }

      /** Construct popup menu */
      DG.Menu.popup()
        .items([ITEM_CAPTION.LINE,
                ITEM_CAPTION.VERT_LINE,
                ITEM_CAPTION.HORZ_LINE], onClickAction)
        .separator()
        .items([ITEM_CAPTION.VERT_BAND,
                ITEM_CAPTION.HORZ_BAND], onClickAction)
        .show();
    }
  }
}

/**
 * A Dialog window with Formula Lines list, preview and editor.
 */
export class FormulaLinesDialog {
  title: string = 'Formula Lines';
  helpUrl: string = '/help/develop/how-to/show-formula-lines.md';

  host: Host;
  preview: Preview;
  editor: Editor;
  viewerTable?: Table;
  dframeTable?: Table;
  creationControl: CreationControl;

  tabs: DG.TabControl;
  dialog: DG.Dialog;

  /** Returns the Table corresponding to the current tab in the tab control */
  get currentTable(): Table {
    return this.tabs.currentPane.name == ITEM_SOURCE.VIEWER
      ? this.viewerTable!
      : this.dframeTable!;
  }

  /** Initializes all parameters and opens the Dialog window */
  constructor(src: DG.DataFrame | DG.Viewer) {
    /** Init Host */
    this.host = new Host(src);

    /** Init CreationControl */
    this.creationControl = new CreationControl(
      () => this.preview.axisCols,
      (item: DG.FormulaLine) => {
        this.currentTable.add(item);
        this.editor.update(0);
        this.preview.update(0);
      });

    /** Init Preview */
    this.preview = new Preview(this.host.viewerItems!, src, this.creationControl.popupMenu);
    this.preview.height = 300;

    /** Init Editor */
    this.editor = new Editor(this.host.viewerItems!, this.preview.dataFrame, (itemIdx: number): boolean => {
      this.currentTable.update(itemIdx);
      return this.preview.update(itemIdx);
    });

    this.tabs = DG.TabControl.create();
    this.tabs.root.style.height = '230px';

    /** Init Viewer Table (in the first tab) */
    if (this.host.viewerItems)
      this.tabs.addPane(ITEM_SOURCE.VIEWER, () => {
        this.viewerTable = new Table(this.host.viewerItems!, (itemIdx: number): boolean => {
          this.editor.update(itemIdx);
          return this.preview.update(itemIdx);
        });
        return this.viewerTable.root;
      });

    /** Init DataFrame Table (in the second tab) */
    if (this.host.dframeItems)
      this.tabs.addPane(ITEM_SOURCE.DATAFRAME, () => {
        this.dframeTable = new Table(this.host.dframeItems!, (itemIdx: number): boolean => {
          this.editor.update(itemIdx);
          return this.preview.update(itemIdx);
        });
        return this.dframeTable.root;
      });

    /** Display "Add new" button */
    this.tabs.header.append(this.creationControl.button);

    /** Change data source when switching tabs */
    this.tabs.onTabChanged.subscribe((_) => {
      this.editor.items = this.currentTable.items;
      this.preview.items = this.currentTable.items;
      this.currentTable.setFirstItemAsCurrent();
    });

    this.dialog = ui.dialog({ title: this.title, helpUrl: this.helpUrl });

    /** Init Dialog layout */
    let layout = ui.div([
      ui.block([
        this.tabs.root,
        this.preview.root
      ], { style: { width: '55%', paddingRight: '20px' } }),
      ui.block([
        this.editor.root
      ], { style: { width: '45%' } })
    ]);

    this.dialog
      .add(layout)
      .onOK(this.onOKAction.bind(this))
      .show({width: 850, height: 650});
  }

  onOKAction() {
    this.host.save();
  }
}
