import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {delay} from '@datagrok-libraries/test/src/test';

/** Formula Line types */
const enum ITEM_TYPE {
  LINE = 'line',
  BAND = 'band',
  AREA_REGION_ANNOTATION = 'area',
  FORMULA_REGION_ANNOTATION = 'formula',
}

/** Formula Line captions */
const enum ITEM_CAPTION {
  LINE = 'Line',
  BAND = 'Band',
  CONST_LINE = 'Constant line',
  VERT_LINE = 'Line - Vertical',
  HORZ_LINE = 'Line - Horizontal',
  HORZ_BAND = 'Band - Horizontal',
  VERT_BAND = 'Band - Vertical',
  FORMULA_REGION = 'Region - Formula Lines',
  RECT_REGION = 'Region - Draw Rectangle',
  POLYGON_REGION = 'Region - Draw Lasso',
}

const enum ITEM_ORIENTATION {
    VERTICAL = 'Vertical',
    HORIZONTAL = 'Horizontal',
    NONE = '',
};

/** Formula Line sources */
const enum ITEM_SOURCE {
  VIEWER = 'Viewer',
  DATAFRAME = 'DataFrame',
}

const enum BTN_CAPTION {
  ADD_NEW = 'Add new',
  CLONE = 'Clone',
  REMOVE = 'Remove',
  FORMULA_LINES_HISTORY = 'Formula Lines History',
  ANNOTATION_REGIONS_HISTORY = 'Annotation Regions History',
  EMPTY = 'Empty',
}

type EditorItem = DG.FormulaLine | DG.AnnotationRegion;

export const DEFAULT_OPTIONS: EditorOptions = {
  allowEditDFLines: true,
};

const HISTORY_KEY = 'formula-lines-dialog';
const HISTORY_KEY_ANNOTATIONS = 'annotation-regions-dialog';
const HISTORY_LENGTH = 12;

/**
 * Returns Formula Line type by its user-friendly [caption].
 * @param {ITEM_CAPTION} caption
 * @return {ITEM_TYPE}
 */
function getItemTypeByCaption(caption: string): string {
  switch (caption) {
    case ITEM_CAPTION.LINE:
    case ITEM_CAPTION.HORZ_LINE:
    case ITEM_CAPTION.VERT_LINE: return ITEM_TYPE.LINE;

    case ITEM_CAPTION.BAND:
    case ITEM_CAPTION.HORZ_BAND:
    case ITEM_CAPTION.VERT_BAND: return ITEM_TYPE.BAND;

    default: throw new Error('Unknown item caption.');
  }
}

const isAnnotationRegionType = (type: string | undefined): boolean =>
  type === ITEM_TYPE.AREA_REGION_ANNOTATION || type === ITEM_TYPE.FORMULA_REGION_ANNOTATION;

const formatAreaFormula = (item: DG.AnnotationRegion): string => {
  if (!item)
    return '';

  if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
    const region = item as DG.AreaAnnotationRegion;
    return `(${region.x}, ${region.y}): ${JSON.stringify(region.area)}`;
  }

  const region = item as DG.FormulaAnnotationRegion;
  return `${region.formula1}; ${region.formula2}`;
}

/**
 * Formula Lines Host Helper.
 * Reads, storages and saves Formula Lines to the host (DataFrame and Viewer).
 */
class Host {
  private dframeFormulaLinesHelper?: DG.FormulaLinesHelper;
  private viewerFormulaLinesHelper?: DG.FormulaLinesHelper;
  private dfAnnotationRegionsHelper?: DG.AnnotationRegionsHelper;
  private viewerAnnotationRegionsHelper?: DG.AnnotationRegionsHelper;

  public dframeFormulaLineItems?: DG.FormulaLine[];
  public viewerFormulaLineItems?: DG.FormulaLine[];
  public dframeAnnotationRegionItems?: DG.AnnotationRegion[];
  public viewerAnnotationRegionItems?: DG.AnnotationRegion[];

  public get dframeItems(): EditorItem[] {
    return [ ...(this.dframeFormulaLineItems ?? []), ...(this.dframeAnnotationRegionItems ?? []) ];
  }

  public get viewerItems(): EditorItem[] {
    return [ ...(this.viewerFormulaLineItems ?? []), ...(this.viewerAnnotationRegionItems ?? []) ];
  }

  constructor(src: DG.DataFrame | DG.Viewer) {
    if (!src)
      throw new Error('Source table/viewer not found.');

    const isViewer = src instanceof DG.Viewer;
    if (isViewer)
      this.initViewerParams(src);

    const df = isViewer ? src.dataFrame : src;
    if (df)
      this.initDfParams(df);
    else if (isViewer)
      throw new Error('Viewer not attached to table.');
  }
  
  private initDfParams(src: DG.DataFrame | DG.Viewer) {
    this.dframeFormulaLinesHelper = src.meta.formulaLines;
    this.dfAnnotationRegionsHelper = src.meta.annotationRegions;
    this.dframeFormulaLineItems = this.dframeFormulaLinesHelper?.items;
    this.dframeAnnotationRegionItems = this.dfAnnotationRegionsHelper?.items;
  }

  private initViewerParams(src: DG.DataFrame | DG.Viewer) {
    this.viewerFormulaLinesHelper = src.meta.formulaLines;
    this.viewerAnnotationRegionsHelper = src.meta.annotationRegions;
    this.viewerFormulaLineItems = this.viewerFormulaLinesHelper.items;
    this.viewerAnnotationRegionItems = this.viewerAnnotationRegionsHelper.items;
  }

  public save() {
    this.dframeFormulaLinesHelper?.clear();
    this.dframeFormulaLinesHelper?.addAll(this.dframeFormulaLineItems!);

    this.dfAnnotationRegionsHelper?.clear();
    this.dfAnnotationRegionsHelper?.addAll(this.dframeAnnotationRegionItems!);

    this.viewerFormulaLinesHelper?.clear();
    this.viewerFormulaLinesHelper?.addAll(this.viewerFormulaLineItems!);

    this.viewerAnnotationRegionsHelper?.clear();
    this.viewerAnnotationRegionsHelper?.addAll(this.viewerAnnotationRegionItems!);
  }
}

/**
 * Table Helper for displaying and navigating Formula Lines list.
 */
class Table {
  public get root(): HTMLElement {return this.grid.root;}
  public get currentItem(): EditorItem | null {
    return this.currentItemIdx < 0
      ? null
      : (this.currentItemIdx >= this.formulaLineItems.length
        ? this.annotationRegionItems[this.currentItemIdx - this.formulaLineItems.length]
        : this.formulaLineItems[this.currentItemIdx]);
  }

  private grid: DG.Grid;

  private get dataFrame(): DG.DataFrame { return this.grid.dataFrame!; }

  public get currentItemIdx(): number { return this.dataFrame.currentRowIdx; }
  public set currentItemIdx(rowIdx: number) { this.dataFrame.currentRowIdx = rowIdx; }

  /** Used to prevent onValuesChanged event when the grid changes itself */
  public notify: boolean = true;

  /** Creates "Delete" button */
  private deleteBtn(itemIdx: number): HTMLElement {
    const btn = ui.button(ui.iconFA('trash-alt'), () => {
      const ifFormulaLine = itemIdx >= 0 && itemIdx < this.formulaLineItems.length;
      (ifFormulaLine ? this.formulaLineItems : this.annotationRegionItems).splice(itemIdx, 1);
      this.dataFrame.rows.removeAt(itemIdx);
      if (this.currentItemIdx > itemIdx)
        this.currentItemIdx--;

      this.onItemChangedAction(this.currentItemIdx);
    });
    btn.style.textAlign = 'center';
    btn.style.height = '20px';
    return btn;
  }

  constructor(
    public formulaLineItems: DG.FormulaLine[],
    public annotationRegionItems: DG.AnnotationRegion[],
    private onItemChangedAction: Function, srcAxes?: AxisNames,
    setCurrentByAxes?: boolean
  ) {
    const dataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'title', []),
      DG.Column.fromList('string', 'formula', []),
      DG.Column.fromList('bool', 'visible', []),
    ]);

    for (let i = 0; i < this.formulaLineItems.length; i++) {
      dataFrame.rows.addNew([this.formulaLineItems[i].title ?? '',
        this.formulaLineItems[i].formula,
        this.formulaLineItems[i].visible,
      ]);
    }

    for (let i = 0; i < this.annotationRegionItems.length; i++) {
      dataFrame.rows.addNew([this.annotationRegionItems[i].header ?? '',
        formatAreaFormula(this.annotationRegionItems[i]),
        !this.annotationRegionItems[i].hidden,
      ]);
    }

    /** Column for "trash" buttons */
    dataFrame.columns.addNewString(BTN_CAPTION.REMOVE);

    this.grid = DG.Grid.create(dataFrame);
    this.grid.setOptions({
      showCurrentRowIndicator: true,
      showSelectedRows: false,
      allowRowResizing: false,
      allowBlockSelection: false,
      allowColSelection: false,
      allowRowReordering: false,
      showRowHeader: false,
      showCellTooltip: false,
      showColumnTooltip: false,
      showCurrentCellOutline: false,
      showContextMenu: false,
      showEditRow: false,
    });

    this.grid.columns.setVisible(['title', 'formula', 'visible', BTN_CAPTION.REMOVE]);
    this.grid.columns.setOrder(['title', 'formula', 'visible', BTN_CAPTION.REMOVE]);

    this.grid.col('title')!.width = 120;
    this.grid.col('formula')!.width = 220;
    this.grid.col('visible')!.width = 40;

    const deleteBtnCol = this.grid.col(BTN_CAPTION.REMOVE)!;
    deleteBtnCol.width = 35;
    deleteBtnCol.cellType = 'html';

    this.grid.onCellPrepare((cell) => {
      if (cell.isColHeader)
        switch (cell.gridColumn.name) {
          case 'title': cell.customText = 'Title'; break;
          case 'formula': cell.customText = 'Formula'; break;
          case 'visible': cell.customText = 'Show'; break;
          case BTN_CAPTION.REMOVE: cell.style.textColor = 0; break;
        }
      else if (cell.isTableCell && cell.gridColumn.name === BTN_CAPTION.REMOVE)
        cell.style.element = this.deleteBtn(cell.gridRow);
    });

    this.dataFrame.onCurrentRowChanged.subscribe((_) => this.onItemChangedAction(this.currentItemIdx));

    this.dataFrame.onValuesChanged.subscribe((ed) => {
      if (!this.notify)
        return;

      if (ed.args?.indexes && ed.args.indexes.length === 1 && ed.args.indexes[0] !== -1)
        this.currentItemIdx = ed.args.indexes[0];

      const item = this.currentItem
      if (item && isAnnotationRegionType(item.type)) {
        const item = this.annotationRegionItems[this.currentItemIdx - this.formulaLineItems.length];
        item.header = this.dataFrame.get('title', this.currentItemIdx);
        item.hidden = !this.dataFrame.get('visible', this.currentItemIdx);
        if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
          const data = this.dataFrame.get('formula', this.currentItemIdx)
            .replace(`(${(item as DG.AreaAnnotationRegion).x}, ${(item as DG.AreaAnnotationRegion).y}): `, '');
          try {
            (item as DG.AreaAnnotationRegion).area = JSON.parse(data || '[]');
          }
          catch {
            (item as DG.AreaAnnotationRegion).area = [];
          }
        } else {
          const data = this.dataFrame.get('formula', this.currentItemIdx).split('; ');
          (item as DG.FormulaAnnotationRegion).formula1 = data[0] ?? '';
          (item as DG.FormulaAnnotationRegion).formula2 = data[1] ?? '';
        }

        this.onItemChangedAction(this.currentItemIdx);
      } else if (item) {
        (item as DG.FormulaLine).title = this.dataFrame.get('title', this.currentItemIdx);
        (item as DG.FormulaLine).formula = this.dataFrame.get('formula', this.currentItemIdx);
        (item as DG.FormulaLine).visible = this.dataFrame.get('visible', this.currentItemIdx);
        this.onItemChangedAction(this.currentItemIdx);
      }

    });

    if (setCurrentByAxes && srcAxes?.x && srcAxes?.y) {
      const checkAxesInFormula = (formula: string) => formula.includes(srcAxes.x!) && formula.includes(srcAxes.y!);
      let itemIdx = formulaLineItems.findIndex((item: DG.FormulaLine) => checkAxesInFormula(item.formula ?? ''));
      if (itemIdx === -1) {
        itemIdx = annotationRegionItems.findIndex((item: DG.AnnotationRegion) => item.type === ITEM_TYPE.AREA_REGION_ANNOTATION
          ? (item as DG.AreaAnnotationRegion).x === srcAxes.x! && (item as DG.AreaAnnotationRegion).y === srcAxes.y!
          : checkAxesInFormula((item as DG.FormulaAnnotationRegion).formula1 ?? '') || checkAxesInFormula((item as DG.FormulaAnnotationRegion).formula2 ?? ''));
        
        if (itemIdx !== -1)
          itemIdx += formulaLineItems.length;
      }

      if (itemIdx !== -1) {
        delay(1).then((_) => {
          this.currentItemIdx = itemIdx;
        });
      }
    }
  }

  public setFirstItemAsCurrent() {
    if (this.formulaLineItems.length > 0 || this.annotationRegionItems.length > 0) {
      this.currentItemIdx = 0;
      this.onItemChangedAction(0);
    } else
      this.onItemChangedAction(-1);
  }

  public update(itemIdx: number, isFormulaLine: boolean = true) {
    const idx = isFormulaLine ? itemIdx : (itemIdx + this.formulaLineItems.length);
    const item: EditorItem = isFormulaLine ? this.formulaLineItems[itemIdx] : this.annotationRegionItems[itemIdx];
    this.notify = false;
    this.dataFrame.set('title', idx, isFormulaLine ? (item as DG.FormulaLine).title : (item as DG.AnnotationRegion).header);
    this.dataFrame.set('formula', idx, isFormulaLine ? (item as DG.FormulaLine).formula : formatAreaFormula(item as DG.AnnotationRegion));
    this.notify = true;
    this.dataFrame.currentRowIdx = idx;
  }

  public add(item: EditorItem, isFormulaLine: boolean = true) {
    const firstIdx = isFormulaLine ? 0 : this.formulaLineItems.length;
    (isFormulaLine ? this.formulaLineItems : this.annotationRegionItems).unshift(item);
    this.dataFrame.rows.insertAt(firstIdx);
    this.notify = false;
    this.dataFrame.set('title', firstIdx, isFormulaLine ? (item as DG.FormulaLine).title : (item as DG.AnnotationRegion).header);
    this.dataFrame.set('formula', firstIdx, isFormulaLine ? (item as DG.FormulaLine).formula : formatAreaFormula(item as DG.AnnotationRegion));
    this.dataFrame.set('visible', firstIdx, isFormulaLine ? (item as DG.FormulaLine).visible : !(item as DG.AnnotationRegion).hidden);
    this.dataFrame.set(BTN_CAPTION.REMOVE, firstIdx, '');
    this.notify = true;
    this.currentItemIdx = firstIdx;
  }
}

interface AxisNames {
  y?: string,
  x?: string,
  yMap?: string,
  xMap?: string,
}

interface AxisColumns {
  y: DG.Column,
  x: DG.Column,
  yMap?: string,
  xMap?: string,
}

export interface EditorOptions {
  allowEditDFLines: boolean,
  [propertyName: string]: any,
}

/**
 * Preview Helper for Formula Lines.
 * Scatter Plot viewer by default.
 */
class Preview {
  public viewer: DG.ScatterPlotViewer | DG.LineChartViewer;
  public dataFrame: DG.DataFrame;
  
  /** Original data frame (used for line chart to validate columns).*/
  private originalDataFrame?: DG.DataFrame;

  /** Source Scatter Plot axes */
  public srcAxes?: AxisNames;

  public set height(h: number) {this.viewer.root.style.height = `${h}px`;}
  public get root(): HTMLElement {return this.viewer.root;}

  /** Returns the current columns pair of the preview Scatter Plot */
  public get axisCols(): AxisColumns {
    let yColName;
    if (this.viewer instanceof DG.LineChartViewer) {
      const yCols: string[] = this.viewer.props.yColumnNames;
      yColName = this.dataFrame.columns.toList().find((col) => yCols.some((n) => col.name.includes(n)))?.name;
    }
    else
      yColName = (this.viewer as DG.ScatterPlotViewer).props.yColumnName;

    const xMap = this.viewer.props.xMap;
    const yMap = this.viewer instanceof DG.ScatterPlotViewer
      ? this.viewer.props.yMap
      : undefined;

    return {
      y: this.dataFrame.getCol(yColName!),
      x: this.dataFrame.getCol(this.viewer.props.xColumnName),
      yMap,
      xMap,
    };
  }

  /** Sets the current axes of the preview Scatter Plot by column names */
  private set axes(names: AxisNames) {
    const options: {  [x: string]: any} = {};
    if (names?.y && this.dataFrame.getCol(names.y))
      if (this.viewer.type === DG.VIEWER.LINE_CHART)
        options['yColumnNames'] = [names.y]
      else {
        options['yColumnName'] = names.y;
        options['yMap'] = names.yMap;
      }

    if (names?.x && this.dataFrame.getCol(names.x)) {
      options['xColumnName'] = names.x;
      options['xMap'] = names.xMap;
    }

    this.viewer.setOptions(options);
  }

  /**
   * Extracts the axes names from the formula. If possible, adjusts the axes
   * of the formula to the axes of the original scatterplot.
   */
  private getItemAxes(axesItem: EditorItem): AxisNames {
    if (isAnnotationRegionType(axesItem.type)) {
      if (axesItem.type === ITEM_TYPE.AREA_REGION_ANNOTATION)
        return {
          y: (axesItem as DG.AreaAnnotationRegion).y,
          x: (axesItem as DG.AreaAnnotationRegion).x,
          yMap: (axesItem as DG.AreaAnnotationRegion).yMap,
          xMap: (axesItem as DG.AreaAnnotationRegion).xMap,
        };
      
        const item = axesItem as DG.FormulaAnnotationRegion;
        const meta1 = item.formula1
          ? DG.FormulaLinesHelper.getMetaByFormula(item.formula1, ITEM_TYPE.LINE)
          : null;

        const meta2 = item.formula2
          ? DG.FormulaLinesHelper.getMetaByFormula(item.formula2, ITEM_TYPE.LINE)
          : null;
        
        return {
          y: meta1?.funcName ?? meta2?.funcName,
          x: meta1?.argName ?? meta2?.argName,
          yMap: item.yMap,
          xMap: item.xMap,
        };
    }

    const item = axesItem as DG.FormulaLine;
    const itemMeta = DG.FormulaLinesHelper.getMeta(item);
    const result: AxisNames = {
      y: item.orientation === ITEM_ORIENTATION.VERTICAL ? itemMeta.argName : itemMeta.funcName,
      x: item.orientation === ITEM_ORIENTATION.VERTICAL ? itemMeta.funcName : itemMeta.argName,
      xMap: item.xMap,
      yMap: item.yMap,
    };
    
    /** If the source axes exist, then we try to set similar axes */
    if (this.srcAxes) {
      result.y ??= this.srcAxes.y;
      result.x ??= this.srcAxes.x;
      result.yMap ??= this.srcAxes.yMap;
      result.xMap ??= this.srcAxes.xMap;

      if (result.x === this.srcAxes.y || result.y === this.srcAxes.x) {
        const tmp = result.x;
        result.x = result.y;
        result.y = tmp;
        
        const tmpMap = result.x;
        result.xMap = result.yMap;
        result.yMap = tmpMap;
      }

      if (result.x === result.y) {
        result.y = this.srcAxes.y;
        result.yMap = this.srcAxes.yMap;
      }
    }

    return result;
  }

  constructor(public formulaLineItems: DG.FormulaLine[],
    public annotationRegionItems: DG.AnnotationRegion[],
    src: DG.DataFrame | DG.Viewer, onContextMenu: Function
  ) {
    if (src instanceof DG.DataFrame)
      this.dataFrame = src;
    else if (src instanceof DG.Viewer) {
      if (src.getOptions()['type'] === DG.VIEWER.LINE_CHART) {
        const viewer = new DG.LineChartViewer(src.dart);
        this.dataFrame = new DG.DataFrame(viewer.activeFrame!);
        this.originalDataFrame = src.dataFrame;
        const yCols: string[] = src.props.yColumnNames;
        const yCol = this.dataFrame.columns.toList()
          .find((col) => col.isNumerical && col.name !== src.props.xColumnName && yCols.some((n) => col.name.includes(n)));
        this.srcAxes = {x: src.props.xColumnName, xMap: src.props.xMap, y: yCol === undefined ? src.props.xColumnName : yCol.name};
      } else if (src.getOptions()['type'] === DG.VIEWER.TRELLIS_PLOT) {
        this.dataFrame = src.dataFrame!;
        const innerLook = src.getOptions()['look']['innerViewerLook'];
        this.srcAxes = {y: innerLook['yColumnName'], x: innerLook['xColumnName'], yMap: innerLook['yMap'], xMap: innerLook['xMap']};
      } else {
        this.dataFrame = src.dataFrame!;
        this.srcAxes = {y: src.props.yColumnName, x: src.props.xColumnName, yMap: src.props.yMap, xMap: src.props.xMap};
      }
    } else
      throw new Error('Host is not DataFrame or Viewer.');

    if (src instanceof DG.Viewer && src.getOptions()['type'] === DG.VIEWER.LINE_CHART)
      this.viewer = DG.Viewer.lineChart(this.dataFrame, {
        yAxisType: src.props.yAxisType,
        xAxisType: src.props.xAxisType,
        splitColumnNames: src.props.splitColumnNames,
        invertXAxis: src.props.invertXAxis,
        showLabels: 'Never',
        showDataframeFormulaLines: false,
        showViewerFormulaLines: true,
        showDataframeAnnotationRegions: false,
        showViewerAnnotationRegions: true,
        showContextMenu: false,
        axesFollowFilter: false,
        axisFont: '11px Arial',
        legendVisibility: DG.VisibilityMode.Never,
        xAxisHeight: 25,
      });
    else {
      this.viewer = DG.Viewer.scatterPlot(this.dataFrame, {
        ...(src as DG.ScatterPlotViewer).props,
        yAxisType: src instanceof DG.Viewer && src.getOptions()['type'] === DG.VIEWER.SCATTER_PLOT ? src.props.yAxisType : 'linear',
        xAxisType: src instanceof DG.Viewer && src.getOptions()['type'] === DG.VIEWER.SCATTER_PLOT ? src.props.xAxisType : 'linear',
        invertXAxis: src instanceof DG.Viewer && src.getOptions()['type'] === DG.VIEWER.SCATTER_PLOT ? src.props.invertXAxis : false,
        invertYAxis: src instanceof DG.Viewer && src.getOptions()['type'] === DG.VIEWER.SCATTER_PLOT ? src.props.invertYAxis : false,
        showDataframeFormulaLines: false,
        showViewerFormulaLines: true,
        showDataframeAnnotationRegions: false,
        showViewerAnnotationRegions: true,
        showSizeSelector: false,
        showColorSelector: false,
        showContextMenu: false,
        axesFollowFilter: false,
        showMinMaxTickmarks: false,
        showMouseOverPoint: false,
        showCurrentPoint: false,
        zoomAndFilter: 'no action',
        axisFont: '11px Arial',
        legendVisibility: DG.VisibilityMode.Never,
        xAxisHeight: 25,
      });
    }

    if (this.srcAxes)
      this.axes = this.srcAxes;

    /**
     * Creates special context menu for preview Scatter Plot.
     * Before opening the menu, it calculates the world coordinates of the click point.
     */
    this.viewer.root.addEventListener('contextmenu', (event: MouseEvent) => {
      event.preventDefault();
      if (this.viewer instanceof DG.ScatterPlotViewer || this.viewer instanceof DG.LineChartViewer) {
        const worldPoint = this.viewer.screenToWorld(event.offsetX, event.offsetY);
        onContextMenu(worldPoint.y, worldPoint.x);
      } 
    });
  }

  /**
   * Shows a line with [itemIdx] index on the Scatter Plot.
   * Returns true if the rendering was successful, false otherwise.
   */
  public update(itemIdx: number, isFormulaLine: boolean = true): boolean {
    /** If there are no lines, try to set the axes as in the original Scatter Plot. */
    if (itemIdx < 0 && this.srcAxes)
      this.axes = this.srcAxes;
    
    const clearMeta = (): void => {
      this.viewer.meta.annotationRegions.clear();
      this.viewer.meta.formulaLines.clear();  
    }

    if (isFormulaLine) {
      try {
        /** Duplicate the original item to display it even if it's hidden */
        const item = this.formulaLineItems[itemIdx];
        const previewItem = structuredClone(item);
        previewItem.visible = true;
  
        /** Trying to show the item */
        clearMeta();
        this.viewer.meta.formulaLines.add(previewItem);
        this.axes = this.getItemAxes(previewItem);
        return true;
      } catch {
        clearMeta();
        return false;
      }
    } else {
      try {
        /** Duplicate the original item to display it even if it's hidden */
        const item = this.annotationRegionItems[itemIdx];
        const previewItem = structuredClone(item);
        previewItem.hidden = false;
  
        /** Trying to show the item */
        clearMeta();
        this.viewer.meta.annotationRegions.add(previewItem);
        this.axes = this.getItemAxes(previewItem);
        return true;
      } catch {
        clearMeta();
        return false;
      }
    }
  }
}

/**
 * Editor Helper for Formula Lines (form with corresponding inputs).
 */
class Editor {
  public get root(): HTMLElement { return this.form; }

  private form: HTMLElement;
  private columnInput: DG.InputBase<DG.Column | null> | undefined;

  /** The title must be accessible from other inputs, because it may depend on the formula */
  private ibTitle?: DG.InputBase;

  private setTitleIfEmpty(oldFormula: string, newFormula: string) {
    if ([oldFormula, this.getTitleFromFormula(oldFormula)].includes((this.ibTitle!.input as HTMLInputElement).placeholder) ||
      newFormula === this.ibTitle!.value)
      this.ibTitle!.value = newFormula;
  }

  constructor(
    public formulaLineItems: DG.FormulaLine[],
    public annotationRegionItems: DG.AnnotationRegion[],
    private dataFrame: DG.DataFrame,
    public onItemChangedAction: (itemIdx: number, isFormulaLine: boolean) => boolean,
    private onFormulaValidation: (isValid: boolean) => void,
  ) {
    this.form = ui.form([]);
  }

  /** Creates and fills editor for given Formula Line */
  public update(itemIdx: number, isFormulaLine: boolean = true) {
    const newForm = itemIdx >= 0
      ? (isFormulaLine ? this.createFormulaLineForm(itemIdx) : this.annotationRegionForm(itemIdx))
      : ui.div(['No formula line or annotation region selected, add one to edit.'], { classes: 'ui-form', style: {
        marginLeft: '-14px',
        overflowX: 'auto',
        textAlign: 'center',
        marginTop: '6px',
      }});

    this.form.replaceWith(newForm);
    this.form = newForm;
  }

  private annotationRegionForm(itemIdx: number): HTMLElement {
    const item = itemIdx >= 0 ? this.annotationRegionItems[itemIdx] : { type: ITEM_TYPE.AREA_REGION_ANNOTATION };
    const mainPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});
    const formatPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});
    const descriptionPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});

    /** Preparing the "Main" panel */
    if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
      mainPane.append(this.inputAreaColumn(itemIdx, 'x'));
      mainPane.append(this.inputAreaColumn(itemIdx, 'y'));
      mainPane?.append(this.areaPointsInput(itemIdx));
    } else {
      mainPane.append(this.inputAnnotationFormula(itemIdx, 'formula1'));
      mainPane.append(this.inputAnnotationFormula(itemIdx, 'formula2'));
    }

    /** Preparing the "Format" panel */
    formatPane.append(this.areaInputColor(itemIdx, 'Region Color', 'fillColor', DG.Color.toHtml(DG.Color.gray)));
    formatPane.append(this.inputOpacity(itemIdx, false));
    formatPane.append(this.areaInputColor(itemIdx, 'Outline Color', 'outlineColor', DG.Color.toHtml(DG.Color.gray)));
    formatPane.append(this.inputLineWidth(itemIdx));
    
    /** Preparing the "Description" panel */
    descriptionPane.append(this.inputAnnotationHeader(itemIdx));
    descriptionPane.append(this.areaInputColor(itemIdx, 'Header Color', 'headerColor'));
    descriptionPane.append(this.inputDescription(itemIdx, false));

    /** Creating the accordion */
    const combinedPanels = ui.accordion();
    combinedPanels.addPane(item.type === ITEM_TYPE.AREA_REGION_ANNOTATION ? 'Area' : 'Formula', () => mainPane, true);
    combinedPanels.addPane('Format', () => formatPane, true);
    combinedPanels.addPane('Description', () => descriptionPane, true);

    return ui.div([combinedPanels.root]);
  }

  private createFormulaLineForm(itemIdx: number): HTMLElement {
    const item = itemIdx >= 0 ? this.formulaLineItems[itemIdx] : {type: ITEM_TYPE.BAND};
    let caption = ITEM_CAPTION.BAND;
    let [itemY, itemX, expression] = ['', '', ''];

    if (itemIdx >= 0 && item.type !== ITEM_TYPE.BAND) {
      const itemMeta = DG.FormulaLinesHelper.getMeta(item);
      [itemY, itemX, expression] = [itemMeta.funcName!, itemMeta.argName!, itemMeta.expression!];
      caption = itemX ? ITEM_CAPTION.LINE : ITEM_CAPTION.CONST_LINE;
    }

    const mainPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});
    const formatPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});
    const tooltipPane = ui.div([], {classes: 'ui-form', style: {marginLeft: '-20px', overflowX: 'auto'}});

    if (itemIdx >= 0) {
      this.columnInput = undefined;
      /** Preparing the "Main" panel */
      mainPane.append(caption === ITEM_CAPTION.CONST_LINE ?
        this.inputConstant(itemIdx, itemY, expression) :
        this.inputFormula(itemIdx));
      if (caption === ITEM_CAPTION.BAND)
        mainPane.append(this.inputColumn2(itemIdx));

      /** Preparing the "Format" panel */
      formatPane.append(this.inputColor(itemIdx));
      formatPane.append(this.inputOpacity(itemIdx));
      if (caption !== ITEM_CAPTION.BAND)
        formatPane.append(this.inputStyle(itemIdx));
      formatPane.append(this.inputRange(itemIdx));
      formatPane.append(this.inputArrange(itemIdx));

      /** Preparing the "Tooltip" panel */
      tooltipPane.append(this.inputTitle(itemIdx));
      tooltipPane.append(this.inputShowLabels(itemIdx));
      tooltipPane.append(this.inputShowDescriptionInTooltip(itemIdx));
      tooltipPane.append(this.inputDescription(itemIdx));
    }

    /** Creating the accordion */
    const combinedPanels = ui.accordion();
    combinedPanels.addPane(itemIdx >= 0 ? caption : ITEM_CAPTION.LINE, () => mainPane, true);
    combinedPanels.addPane('Format', () => formatPane, true);
    combinedPanels.addPane('Tooltip', () => tooltipPane, true);

    return ui.div([combinedPanels.root]);
  }

  /** Creates textarea for item formula */
  private inputFormula(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const ibFormula = ui.input.textArea('', {value: item.formula ?? '',
      onValueChanged: (value) => {
        const oldFormula = item.formula!;
        item.formula = value;
        const resultOk = this.onItemChangedAction(itemIdx, true);
        this.setFormulaValidationResult(resultOk, value, ibFormula, item.type === ITEM_TYPE.BAND);
        this.setTitleIfEmpty(oldFormula, item.formula);
      }});

    const elFormula = ibFormula.input as HTMLInputElement;
    elFormula.placeholder = 'Formula';
    //elFormula.setAttribute('style', 'width: 360px; height: 60px; margin-right: -6px;');
    elFormula.setAttribute('style', 'height: 60px;');

    ui.tools.initFormulaAccelerators(ibFormula, this.dataFrame);

    return ibFormula.root;
  }

  /** Creates color picker for item color */
  private inputColor(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const ibColor = ui.input.color('Color', {value: item.color ?? '#000000',
      onValueChanged: (value) => {
        item.color = value;
        this.onItemChangedAction(itemIdx, true);
      }});

    const elColor = ibColor.input as HTMLInputElement;
    elColor.placeholder = '#000000';
    // elColor.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibColor.root]);
  }

  private areaInputColor(itemIdx: number, header: string = 'Color', key: keyof DG.AnnotationRegion, defaultColor: string = '#000000'): HTMLElement {
    const item = this.annotationRegionItems[itemIdx] as DG.AnnotationRegion;
    const ibColor = ui.input.color(header, {value: item[key] ? DG.Color.toHtml(item[key] as number) : defaultColor,
      onValueChanged: (value) => {
        (item as any)[key] = DG.Color.fromHtml(value);
        this.onItemChangedAction(itemIdx, false);
      }});
    const elColor = ibColor.input as HTMLInputElement;
    elColor.placeholder = defaultColor;
    return ui.divH([ibColor.root]);
  }

  private inputLineWidth(itemIdx: number): HTMLElement {
    const item = this.annotationRegionItems[itemIdx] as DG.AnnotationRegion;
    
    const elOpacity = ui.element('input');
    elOpacity.type = 'range';
    elOpacity.min = 0;
    elOpacity.max = 20;
    elOpacity.value = item.outlineWidth ?? 1;
    elOpacity.addEventListener('input', () => {
      item.outlineWidth = parseInt(elOpacity.value);
      this.onItemChangedAction(itemIdx, false);
    });
    // elOpacity.setAttribute('style', 'width: 204px; margin-top: 6px; margin-left: 0px;');
    elOpacity.setAttribute('style', 'margin-top: 6px; width: 100%;');

    const label = ui.label('Outline Width', 'ui-label ui-input-label');

    return ui.divH([ui.div([label, elOpacity], 'ui-input-root')]);
  }

  /** Creates range slider for item opacity */
  private inputOpacity(itemIdx: number, isFormulaLine: boolean = true): HTMLElement {
    const item = isFormulaLine ? this.formulaLineItems[itemIdx] : this.annotationRegionItems[itemIdx];
    const elOpacity = ui.element('input');
    elOpacity.type = 'range';
    elOpacity.min = 0;
    elOpacity.max = 100;
    elOpacity.value = item.opacity ?? 30;
    elOpacity.addEventListener('input', () => {
      item.opacity = parseInt(elOpacity.value);
      this.onItemChangedAction(itemIdx, isFormulaLine);
    });
    // elOpacity.setAttribute('style', 'width: 204px; margin-top: 6px; margin-left: 0px;');
    elOpacity.setAttribute('style', 'margin-top: 6px; width: 100%;');

    const label = ui.label('Opacity', 'ui-label ui-input-label');

    return ui.divH([ui.div([label, elOpacity], 'ui-input-root')]);
  }


  /** Creates combobox for item line style and text input for item width */
  private inputStyle(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const ibStyle = ui.input.choice('Style', {value: item.style ?? 'solid',
      items: ['solid', 'dotted', 'dashed', 'longdash', 'dotdash'], onValueChanged: (value) => {
        item.style = value;
        this.onItemChangedAction(itemIdx , true);
      }});

    const elStyle = ibStyle.input as HTMLInputElement;
    //elStyle.style.width = '135px';

    const ibWidth = ui.input.int('', {value: item.width ?? 1,
      onValueChanged: (value) => {
        item.width = value;
        this.onItemChangedAction(itemIdx, true);
      }});
    ibWidth.addPostfix('px');

    const elWidth = ibWidth.input as HTMLInputElement;
    elWidth.placeholder = '1';
    elWidth.setAttribute('style', 'width: 61px; padding-right: 24px;');

    // const unit = ui.divText('px', {style: {marginTop: '10px', marginLeft: '-24px', zIndex: '1000'} });

    return ui.divH([ibStyle.root, ibWidth.root]);
  }

  /** Creates text inputs for min-max values of item */
  private inputRange(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const ibMin = ui.input.string('Range', {value: `${item.min ?? ''}`,
      onValueChanged: (value) => {
        item.min = value.length === 0 ? undefined : Number(value);
        this.onItemChangedAction(itemIdx, true);
      }});

    const elMin = ibMin.input as HTMLInputElement;
    elMin.placeholder = 'min';
    elMin.setAttribute('style', 'width: 98px;');

    const ibMax = ui.input.string('', {value: `${item.max ?? ''}`,
      onValueChanged: (value) => {
        item.max = value.length === 0 ? undefined : Number(value);
        this.onItemChangedAction(itemIdx, true);
      }});

    const elMax = ibMax.input as HTMLInputElement;
    elMax.placeholder = 'max';
    elMax.setAttribute('style', 'width: 98px;');

    return ui.divH([ibMin.root, ibMax.root]);
  }

  /** Creates combobox for item position (z-index) */
  private inputArrange(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const ibArrange = ui.input.choice('Arrange', {
      value: item.zIndex && item.zIndex > 0 ? 'above markers' : 'below markers', items: ['above markers', 'below markers'],
      onValueChanged: (value) => {
        item.zIndex = value === 'above markers' ? 100 : -100;
        this.onItemChangedAction(itemIdx, true);
      }});

    const elArrange = ibArrange.input as HTMLInputElement;
    elArrange.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibArrange.root]);
  }

  private getTitleFromFormula(formula: string): string {
    let title = formula;
    if (title) {
      const regexp = /\${(.+?)}/gi;
      const matches: string[][] = [...title.matchAll(regexp)].map((m) => [m[0], m[1]]);
      for (const match of matches)
        title = title!.replace(match[0], match[1]);
    }
    return title;
  }

  private inputAnnotationFormula(itemIdx: number, title: keyof DG.FormulaAnnotationRegion): HTMLElement {
    const item = this.annotationRegionItems[itemIdx] as DG.FormulaAnnotationRegion;

    const ibHeader = ui.input.string(title === 'formula1' ? 'Formula 1' : 'Formula 2', {value: item[title] ? item[title] as string : '',
      onValueChanged: (value) => {
        (item as any)[title] = value;
        const resultOk = this.onItemChangedAction(itemIdx, false);
        this.setFormulaValidationResult(resultOk, value, ibHeader);
      }});

    const elHeader = ibHeader.input as HTMLInputElement;
    elHeader.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibHeader.root]);
  }

  private setFormulaValidationResult(resultOk: boolean, value: string, ibHeader: DG.InputBase<string>, isBand: boolean = false): void {
    const elHeader = ibHeader.input as HTMLInputElement;
    const validateValue = (): string => {
      if (isBand)
        return resultOk ? '' : 'Invalid formula syntax';

      // Must contain exactly one '='
      const parts = value.split('=');
      if (parts.length !== 2)
        return 'Line formula should be in format: ${x or y column} = expression';

      const [lhsRaw, rhsRaw] = parts.map(p => p.trim());

      // LHS must be exactly one ${...}
      const lhsMatch = /^\$\{([^}]+)\}$/.exec(lhsRaw);
      if (!lhsMatch)
        return 'Left side must be a single column in format ${column}';

      const lhsColumn = lhsMatch[1].trim();

      // Extract all ${...} on RHS
      const rhsColumnRegex = /\$\{([^}]+)\}/g;
      const rhsColumns: string[] = [];

      let match;
      while ((match = rhsColumnRegex.exec(rhsRaw)) !== null) {
        rhsColumns.push(match[1].trim());
      }

      // If RHS references the same column → invalid
      if (rhsColumns.includes(lhsColumn))
        return 'Line formula should contain different columns on both sides of the "=" sign';

      // Expression syntax validation comes last
      return resultOk ? '' : 'Invalid formula syntax';
    };
    
    const validationTooltip = validateValue();
    resultOk = resultOk && validationTooltip === '';
    ibHeader.setTooltip(validationTooltip);
    elHeader.classList.toggle('d4-forced-invalid', !resultOk);
    this.onFormulaValidation(resultOk);
  }

  /** Creates text input for item title */
  private inputTitle(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;
    if (!item.title || item.title === item.formula)
      item.title = this.getTitleFromFormula(item.formula!);
    const thisElem = this;

    function formTitleValue(value: string) {
      const potentialValue = thisElem.getTitleFromFormula(item.formula!);
      if (value === '' || value === item.formula || value === potentialValue) {
        elTitle.placeholder = potentialValue;
        elTitle.value = '';
        return potentialValue;
      }
      elTitle.placeholder = '';
      return value;
    }

    this.ibTitle = ui.input.string('Title', {value: item.title ?? '',
      onValueChanged: (value) => {
        item.title = formTitleValue(value);
        this.onItemChangedAction(itemIdx, true);
      }});

    const elTitle = this.ibTitle.input as HTMLInputElement;
    elTitle.setAttribute('style', 'width: 204px; max-width: none;');
    formTitleValue(elTitle.value);

    return ui.divH([this.ibTitle.root]);
  }

  private inputAnnotationHeader(itemIdx: number): HTMLElement {
    const item = this.annotationRegionItems[itemIdx];

    const ibHeader = ui.input.string('Title', {value: item.header ?? '',
      onValueChanged: (value) => {
        item.header = value;
        this.onItemChangedAction(itemIdx, false);
      }});

    const elHeader = ibHeader.input as HTMLInputElement;
    elHeader.setAttribute('style', 'width: 204px; max-width: none;');

    return ui.divH([ibHeader.root]);
  }

  /** Creates show on plot bool */
  private inputShowLabels(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const iShowLabels = ui.input.bool('Show on plot', {value: item.showOnPlot ?? true,
      onValueChanged: (value) => {
        item.showOnPlot = value;
        this.onItemChangedAction(itemIdx, true);
      }});


    return iShowLabels.root;
  }

  /** Creates show on plot bool */
  inputShowDescriptionInTooltip(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    const iShowLabels = ui.input.bool('Show on tooltip', {value: item.showOnTooltip ?? true,
      onValueChanged: (value) => {
        item.showOnTooltip = value;
        this.onItemChangedAction(itemIdx, true);
      }});


    return iShowLabels.root;
  }

  private setPointsValidationResult(resultOk: boolean, value: string, ibPoints: DG.InputBase<string>): void {
    const elPoints = ibPoints.input as HTMLInputElement;
    const tooltipWarning = 'Points should be a list of [x, y] number pairs, e.g., [1, 2], [3, 4]';
    const validateValue = (): string => {
      if (resultOk)
        return '';
      
      const parsed = JSON.parse(`[${value}]`);
      return Array.isArray(parsed) && parsed.length < 3 ? 'Area must have at least 3 points' : tooltipWarning;
    };
    
    try {
      ibPoints.setTooltip(validateValue());
    } catch {
      ibPoints.setTooltip(tooltipWarning);
    }

    elPoints.classList.toggle('d4-forced-invalid', !resultOk);
    this.onFormulaValidation(resultOk);
  }

  private areaPointsInput(itemIdx: number): HTMLElement {
    const item = this.annotationRegionItems[itemIdx] as DG.AreaAnnotationRegion;
    const value = item.area ? JSON.stringify(item.area).replaceAll(',', ', ') :  '';
    const textArea = ui.input.textArea('Points', {
      value: value.substring(1, value.length - 1),
      onValueChanged: (value) => {
        let resultOk = true;
        try {
          var parsed = JSON.parse(`[${value}]`);
          if (!Array.isArray(parsed)) {
            resultOk = false;
          } else {
            const filtered = parsed.filter((p) => Array.isArray(p) && p.length === 2
              && typeof p[0] === 'number' && typeof p[1] === 'number');
            if (filtered.length !== parsed.length || filtered.length < 3) {
              resultOk = false;
            } else {
              item.area = filtered;
              this.onItemChangedAction(itemIdx, false);
            }
          }
        } catch {
          resultOk = false;
        }

        this.setPointsValidationResult(resultOk, value, textArea);
      }});

    const elDescription = textArea.input as HTMLInputElement;
    elDescription.setAttribute('style',
      'height: 75px; font-family: inherit; font-size: inherit;');
    // 'width: 194px; height: 40px; padding-left: 6px; margin-right: -8px; font-family: inherit; font-size: inherit;');

    return textArea.root;
  }
  

  /** Creates textarea for item description */
  private inputDescription(itemIdx: number, isFormulaLine: boolean = true): HTMLElement {
    const item = isFormulaLine ? this.formulaLineItems[itemIdx] : this.annotationRegionItems[itemIdx];

    const ibDescription = ui.input.textArea('Description', {value: item.description ?? '',
      onValueChanged: (value) => {
        item.description = value;
        this.onItemChangedAction(itemIdx, isFormulaLine);
      }});

    const elDescription = ibDescription.input as HTMLInputElement;
    elDescription.setAttribute('style',
      'height: 40px; font-family: inherit; font-size: inherit;');
    // 'width: 194px; height: 40px; padding-left: 6px; margin-right: -8px; font-family: inherit; font-size: inherit;');

    return ibDescription.root;
  }

   /** Creates column input for band second column */
  private inputAreaColumn(itemIdx: number, type: 'x' | 'y'): HTMLElement {
    const item = this.annotationRegionItems[itemIdx] as DG.AreaAnnotationRegion;

    const mapKey = type + 'Map' as keyof DG.AreaAnnotationRegion;
    const itemCol = item[type] ?? '';
    const colName = item[mapKey] && itemCol && itemCol.endsWith(item[mapKey] as string)
      ? itemCol.substring(0, itemCol.length - (item[mapKey] as string).length - 1)
      : item[type];

    const ibColumn2 = ui.input.column(type.toUpperCase() + ' column', {
      nullable: false,
      table: this.dataFrame,
      value: this.dataFrame.col(colName ?? '') ?? undefined,
      onValueChanged: (value) => {
        item[type] = item[mapKey] && value?.name ? `${value.name} ${item[mapKey]}` : value?.name;
        this.onItemChangedAction(itemIdx, false);
      }});
      
    this.columnInput = ibColumn2;
    const elColumn2 = ibColumn2.input as HTMLInputElement;
    //elColumn2.setAttribute('style', 'width: 204px; max-width: none;');
    
    return ui.divH([ibColumn2.root]);
  }

  public inputColumn2Changing: boolean = false;

  /** Creates column input for band second column */
  private inputColumn2(itemIdx: number): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    //@ts-ignore
    const ibColumn2 = ui.input.column('Adjacent column', {
      nullable: false,
      table: this.dataFrame,
      value: item.column2 ? this.dataFrame.col(item.column2) ?? undefined : undefined,
      onValueChanged: (value) => {
        this.inputColumn2Changing = true;
        item.column2 = value.name;
        this.onItemChangedAction(itemIdx, true);
        this.inputColumn2Changing = false;
      }});
      
    this.columnInput = ibColumn2;
    const elColumn2 = ibColumn2.input as HTMLInputElement;
    //elColumn2.setAttribute('style', 'width: 204px; max-width: none;');
    
    return ui.divH([ibColumn2.root]);
  }

  /** Creates column input and text input for constant item */
  private inputConstant(itemIdx: number, colName: string, value: string): HTMLElement {
    const item = this.formulaLineItems[itemIdx] as DG.FormulaLine;

    //@ts-ignore
    const ibColumn = ui.input.column('Column', {
      nullable: false,
      table: this.dataFrame,
      value: colName ? this.dataFrame.col(colName) ?? undefined : undefined,
      onValueChanged: (value) => {
        const oldFormula = item.formula!;
        item.formula = '${' + value + '} = ' + ibValue.value;
        this.onItemChangedAction(itemIdx, true);
        this.setTitleIfEmpty(oldFormula, item.formula);
      }});

    this.columnInput = ibColumn;
    const elColumn = ibColumn.input as HTMLInputElement;
    //elColumn.setAttribute('style', 'width: 204px; max-width: none; margin-right: -10px;');

    const ibValue = ui.input.string('Value', {value: value, onValueChanged: (value) => {
      const oldFormula = item.formula!;
      item.formula = '${' + ibColumn.value + '} = ' + value;
      this.onItemChangedAction(itemIdx, true);
      this.setTitleIfEmpty(oldFormula, item.formula);
    }});
    ibValue.nullable = false;

    const elValue = ibValue.input as HTMLInputElement;
    //elValue.setAttribute('style', 'width: 204px; max-width: none; margin-right: -10px;');

    return ui.div([ibColumn.root, ibValue.root], 'ui-form');
  }
}

/**
 * Helper that implements the logic of creating a Formula Line item of a given type.
 */
class CreationControl {
  public popupMenu: Function;        // Opens a popup menu with predefined new Formula Line item types

  /** Items for History menu group */
  public formulaLinesHistoryItems: DG.FormulaLine[];           // Stores session history
  public formulaLinesJustCreatedItems: DG.FormulaLine[] = [];  // Stores history of the currently open dialog

  public annotationRegionsHistoryItems: DG.AnnotationRegion[];
  public annotationRegionsJustCreatedItems: DG.AnnotationRegion[] = [];

  public loadFormulaLinesHistory(): DG.FormulaLine[] {
    return localStorage[HISTORY_KEY] ? JSON.parse(localStorage[HISTORY_KEY]) : [];
  }

  public loadAnnotationRegionsHistory(): DG.AnnotationRegion[] {
    return localStorage[HISTORY_KEY_ANNOTATIONS] ? JSON.parse(localStorage[HISTORY_KEY_ANNOTATIONS]) : []
  }

  public saveHistory() {
    const compareItems = (a: EditorItem, b: EditorItem) => JSON.stringify(a) === JSON.stringify(b);
    /** Remove duplicates from just created items (object comparison via JSON.stringify) */
    this.formulaLinesJustCreatedItems = this.formulaLinesJustCreatedItems.filter((val, ind, arr) =>
      arr.findIndex((t) => compareItems(t, val)) === ind);

    /** Remove identical older items from history */
    this.formulaLinesHistoryItems = this.formulaLinesHistoryItems.filter((arr) =>
      !this.formulaLinesJustCreatedItems.find((val) => compareItems(val, arr)));

    const newHistoryItems = this.formulaLinesJustCreatedItems.concat(this.formulaLinesHistoryItems);
    newHistoryItems.splice(HISTORY_LENGTH);

    localStorage[HISTORY_KEY] = JSON.stringify(newHistoryItems);

    /** Repeat for annotation regions */
    this.annotationRegionsJustCreatedItems = this.annotationRegionsJustCreatedItems.filter((val, ind, arr) =>
      arr.findIndex((t) => compareItems(t, val)) === ind);

    this.annotationRegionsHistoryItems = this.annotationRegionsHistoryItems.filter((arr) =>
      !this.annotationRegionsJustCreatedItems.find((val) => compareItems(val, arr)));

    const newAnnotationHistoryItems = this.annotationRegionsJustCreatedItems.concat(this.annotationRegionsHistoryItems);
    newAnnotationHistoryItems.splice(HISTORY_LENGTH);

    localStorage[HISTORY_KEY_ANNOTATIONS] = JSON.stringify(newAnnotationHistoryItems);
  }

  /** Creates a button and binds an item creation menu to it */
  public get button(): HTMLElement {
    const btn = ui.bigButton(BTN_CAPTION.ADD_NEW, this.popupMenu);
    return ui.div([btn], {style: {width: '100%', textAlign: 'right'}});
  }

  constructor(
    getCols: () => AxisColumns,                              // Used to create constant lines passing through the mouse click point on the Scatter Plot
    getCurrentItem: () => EditorItem | null,                 // Used to create clone
    private onItemCreatedAction: (item: EditorItem) => void,                          // Updates the Table, Preview and Editor states after item creation
    createArea: (lassoMode?: boolean) => Promise<DG.AnnotationRegion | null>  // Used to create area annotation regions
  ) {
    this.formulaLinesHistoryItems = this.loadFormulaLinesHistory();
    this.annotationRegionsHistoryItems = this.loadAnnotationRegionsHistory();

    this.popupMenu = (valY?: number, valX?: number) => {
      const onClickAction = (itemCaption: string) => {
        if (itemCaption === ITEM_CAPTION.POLYGON_REGION || itemCaption === ITEM_CAPTION.RECT_REGION) {
          createArea?.(itemCaption === ITEM_CAPTION.POLYGON_REGION).then((region: DG.AnnotationRegion | null) => {
            if (!region)
              return;

            this.annotationRegionsJustCreatedItems.unshift(region);
            /** Update the Table, Preview and Editor states */
            onItemCreatedAction(region);
          });
          return;
        }

        const cols: AxisColumns = getCols();
        const colY = cols.y;
        const colX = cols.x;
        if (itemCaption === ITEM_CAPTION.FORMULA_REGION) {
          const item: DG.FormulaAnnotationRegion = {
            type: ITEM_TYPE.FORMULA_REGION_ANNOTATION,
            formula1: '${' + colY.name + '} = ${' + colX.name + '}' + ' + ' + colY.stats.q2.toFixed(1),
            formula2: '${' + colY.name + '} = ${' + colX.name + '}' + ' - ' + colY.stats.q2.toFixed(1),
          };

          this.annotationRegionsJustCreatedItems.unshift(item);
          /** Update the Table, Preview and Editor states */
          onItemCreatedAction(item);
          return;
        }

        let item: DG.FormulaLine = {};
        
        /** Fill the item with the necessary data */
        switch (itemCaption) {
          case ITEM_CAPTION.LINE:
            item.formula = '${' + colY.name + '} = ${' + colX.name + '}';
            break;

          case ITEM_CAPTION.VERT_LINE:
            const vertVal = (valX ?? colX.stats.q2).toFixed(1);
            item.orientation = ITEM_ORIENTATION.VERTICAL;
            item.formula = '${' + colX.name + '} = ' + vertVal;
            break;

          case ITEM_CAPTION.HORZ_LINE:
            const horzVal = (valY ?? colY.stats.q2).toFixed(1);
            item.orientation = ITEM_ORIENTATION.HORIZONTAL;
            item.formula = '${' + colY.name + '} = ' + horzVal;
            break;

          case ITEM_CAPTION.VERT_BAND:
            const left = colX.stats.q1.toFixed(1);
            const right = colX.stats.q3.toFixed(1);
            item.formula = '${' + colX.name + '} in(' + left + ', ' + right + ')';
            item.orientation = ITEM_ORIENTATION.VERTICAL;
            item.column2 = colY.name;
            break;

          case ITEM_CAPTION.HORZ_BAND:
            const bottom = colY.stats.q1.toFixed(1);
            const top = colY.stats.q3.toFixed(1);
            item.formula = '${' + colY.name + '} in(' + bottom + ', ' + top + ')';
            item.orientation = ITEM_ORIENTATION.HORIZONTAL;
            item.column2 = colX.name;
            break;

          case BTN_CAPTION.CLONE:
            item = getCurrentItem() as DG.FormulaLine;
            break;
        }

        item.type ??= getItemTypeByCaption(itemCaption);
        
        item = DG.FormulaLinesHelper.setDefaults(item);

        this.formulaLinesJustCreatedItems.unshift(item);

        /** Update the Table, Preview and Editor states */
        onItemCreatedAction(item);
      };

      /** Construct popup menu */
      const menu = DG.Menu.popup().items([
        ITEM_CAPTION.LINE,
        ITEM_CAPTION.VERT_LINE,
        ITEM_CAPTION.HORZ_LINE,
        ITEM_CAPTION.VERT_BAND,
        ITEM_CAPTION.HORZ_BAND,
      ], onClickAction);

      const cols = getCols();
      if (cols.x?.isNumerical && !cols.xMap && cols.y?.isNumerical && !cols.yMap) {
        const regionItems: Record<string, string> = {
          [ITEM_CAPTION.FORMULA_REGION]: 'Adds a new area defined by two formula lines.',
          [ITEM_CAPTION.RECT_REGION]: 'Draws a rectangle area.',
          [ITEM_CAPTION.POLYGON_REGION]: 'Draws a polygon area using lasso tool.',
        }
  
        for (const itemCaption in regionItems)
          menu.item(itemCaption, () => onClickAction(itemCaption), null, {
            description: regionItems[itemCaption]
          });
      }

      /** Add separator only if other menu items exist */
      if (getCurrentItem() || this.formulaLinesHistoryItems?.length || this.annotationRegionsHistoryItems.length)
        menu.separator();

      /**
       * Add "Clone" menu if the current table line exists.
       * TODO: The best option is to make the menu item enabled/disabled. But there is no such API yet.
       */
      if (getCurrentItem())
        menu.items([BTN_CAPTION.CLONE], onClickAction);

      /**
       * Add "History" menu group.
       * TODO: The best option is to make the menu item enabled/disabled. But there is no such API yet.
       */
      if (this.formulaLinesHistoryItems?.length)
        this.fillHistoryGroup(menu.group(BTN_CAPTION.FORMULA_LINES_HISTORY),
          this.formulaLinesHistoryItems, this.formulaLinesJustCreatedItems, (item: DG.FormulaLine) => item.formula!);

      if (this.annotationRegionsHistoryItems?.length)
        this.fillHistoryGroup(menu.group(BTN_CAPTION.ANNOTATION_REGIONS_HISTORY),
          this.annotationRegionsHistoryItems, this.annotationRegionsJustCreatedItems, formatAreaFormula);

      menu.show();
    };
  }

  private fillHistoryGroup(group: DG.Menu, items: EditorItem[], justCreatedItems: EditorItem[], getTitle: (item: EditorItem) => string): void {
    for (const item in items) {
      const title = getTitle(items[item]);
      if (title)
        group.item(title.length > 50 ? title.substring(0, 50) + '...' : title, () => {
          const newItem = structuredClone(items[item]);
          justCreatedItems.unshift(newItem);
          this.onItemCreatedAction(newItem);
        });
    }

    group.endGroup();
  }
}

/**
 * A Dialog window with Formula Lines list, preview and editor.
 */
export class FormulaLinesDialog {
  public dialog: DG.Dialog = ui.dialog({
    title: 'Formula Lines',
    helpUrl: '/help/develop/how-to/viewers/show-formula-lines.md',
  });

  private host: Host;
  private preview: Preview;
  private editor: Editor;
  private viewerTable?: Table;
  private dframeTable?: Table;
  private creationControl: CreationControl;
  private tabs: DG.TabControl;

  /** Returns the Table corresponding to the current tab in the tab control */
  public get currentTable(): Table {
    return this.tabs.currentPane.name === ITEM_SOURCE.VIEWER ? this.viewerTable! : this.dframeTable!;
  }

  /** Initializes all parameters and opens the Dialog window */
  constructor(
    src: DG.DataFrame | DG.Viewer,
    private options: EditorOptions = DEFAULT_OPTIONS,
    private showValueOnOpen?: { index?: number, isDataFrame?: boolean, isAnnotationArea?: boolean })
  {
    /** Init Helpers */
    this.host = this.initHost(src);
    this.creationControl = this.initCreationControl();
    this.preview = this.initPreview(src);
    this.editor = this.initEditor();
    this.tabs = this.initTabs();
    this.dialog.sub(this.dialog.onClose.subscribe(() => this.dialog.detach()));

    /** Init Dialog layout */
    const layout = ui.div([
      ui.block([this.tabs.root, this.preview.root], {style: {width: '55%', paddingRight: '20px'}}),
      ui.block([this.editor.root], {style: {width: '45%'}}),
    ]);

    const width = Math.min(1000, Math.floor(document.body.clientWidth / 1.3));
    const height = Math.min(800, Math.floor(document.body.clientHeight / 1.5));
    this.dialog
      .add(layout)
      .onOK(this.onOKAction.bind(this), {closeOnEnter: false})
      .show({
        resizable: true,
        width,
        height,
        x: Math.floor((window.innerWidth - width) / 2),
        y: Math.floor((window.innerHeight - height) / 2),
      });

    this.initDefaultOnOpenState();

    this.dialog.sub(this.preview.viewer.onPropertyValueChanged.subscribe((typeArgs) => {
      const currentItem = this.currentTable?.currentItem;
      if (!this.editor || this.editor.inputColumn2Changing || !this.preview?.viewer || currentItem?.type !== ITEM_TYPE.BAND)
        return;

      const band = currentItem as DG.FormulaLine;
      const { property } = typeArgs.args as unknown as { property: DG.Property };
      if (!['xColumnName', 'yColumnName', 'yColumnNames'].includes(property.name))
        return;
      
      const isHorz = band.orientation === ITEM_ORIENTATION.HORIZONTAL;
      if (isHorz && property.name === 'xColumnName') {
        band.column2 = this.preview.axisCols.x.name;
        this.editor.update(this.currentTable.currentItemIdx, true);
      } else if (!isHorz && (property.name === 'yColumnName' || property.name === 'yColumnNames')) {
        band.column2 = this.preview.axisCols.y.name;
        this.editor.update(this.currentTable.currentItemIdx, true);
      }
    }));
  }
    
  private initDefaultOnOpenState(): void {      
      if (!this.showValueOnOpen)
        return;

      this.tabs.currentPane = this.tabs.getPane(this.showValueOnOpen.isDataFrame ? ITEM_SOURCE.DATAFRAME : ITEM_SOURCE.VIEWER);
      if (this.showValueOnOpen.index) {
        const isFormulaLine = this.preview.formulaLineItems.length > 0 && this.showValueOnOpen.index < this.preview.formulaLineItems.length;
        this.currentTable.update(this.showValueOnOpen.index, isFormulaLine);
      } else {
        this.currentTable.setFirstItemAsCurrent();
      }
  }

  private initHost(src: DG.DataFrame | DG.Viewer): Host {
    return new Host(src);
  }

  private initPreview(src: DG.DataFrame | DG.Viewer): Preview {
    const preview = new Preview(this.host.viewerFormulaLineItems! ?? this.host.dframeFormulaLineItems!,
      this.host.viewerAnnotationRegionItems! ?? this.host.dframeAnnotationRegionItems!,
      src, this.creationControl.popupMenu);
    preview.height = 310;
    return preview;
  }

  private initCreationControl(): CreationControl {
    return new CreationControl(
      () => this.preview.axisCols,
      () => this.currentTable.currentItem,
      (item: EditorItem) =>
        this.onItemCreatedAction(item, !isAnnotationRegionType(item?.type ?? '')),
      (lassoMode?: boolean) => new Promise<DG.AnnotationRegion | null>((resolve) => {
        if (this.preview.viewer instanceof DG.ScatterPlotViewer || this.preview.viewer instanceof DG.LineChartViewer) {
          this.preview.viewer.disableAnnotationRegionDrawing();
          this.preview.viewer.setOptions({ annotationRegions: '[]', formulaLines: '[]' });
          this.editor.update(-1, false);
          this.preview.viewer.enableAnnotationRegionDrawing(lassoMode, (region: { [key: string]: unknown }) => {
            region['isDataFrameRegion'] = this.tabs.currentPane.name === ITEM_SOURCE.DATAFRAME;
            const props = this.preview.viewer.props as DG.IScatterPlotSettings;
            const annotationRegions = JSON.parse(props.annotationRegions || '[]');
            annotationRegions.push(region);
            props.annotationRegions = JSON.stringify(annotationRegions);
            resolve(region as DG.AnnotationRegion);
          });
        } else
          resolve(null);
      })
    );
  }

  private initEditor(): Editor {
    return new Editor(this.host.viewerFormulaLineItems! ?? this.host.dframeFormulaLineItems!,
      this.host.viewerAnnotationRegionItems! ?? this.host.dframeAnnotationRegionItems!,
      this.preview.dataFrame,
      (itemIdx: number, isFormulaLine: boolean = true): boolean => {
        this.currentTable.update(itemIdx, isFormulaLine);
        return this.preview.update(itemIdx, isFormulaLine);
      },
      (isValid: boolean): void => {
        isValid ? this.dialog.getButton('OK').classList.remove('disabled') : this.dialog.getButton('OK').classList.add('disabled');
      });
  }

  private initTabs(): DG.TabControl {
    const tabs = DG.TabControl.create();
    tabs.root.style.height = '230px';

    /** Init Viewer Table (in the first tab) */
    if (this.host.viewerFormulaLineItems || this.host.viewerAnnotationRegionItems) {
      tabs.addPane(ITEM_SOURCE.VIEWER, () => {
        this.viewerTable = this.initTable(this.host.viewerFormulaLineItems ?? [], this.host.viewerAnnotationRegionItems ?? []);
        return this.viewerTable.root;
      });
    }

    /** Init DataFrame Table (in the second tab) */
    if (this.options.allowEditDFLines && (this.host.dframeFormulaLineItems || this.host.dframeAnnotationRegionItems)) {
      tabs.addPane(ITEM_SOURCE.DATAFRAME, () => {
        this.dframeTable = this.initTable(this.host.dframeFormulaLineItems ?? [], this.host.dframeAnnotationRegionItems ?? []);
        return this.dframeTable.root;
      });
    }
    // Overrides the standard component logic that hides the header containing only one tab
    tabs.header.style.removeProperty('display');

    /** Display "Add new" button */
    tabs.header.append(this.creationControl.button);

    /** Change data source when switching tabs */
    tabs.onTabChanged.subscribe((_) => {
      this.editor.formulaLineItems = this.currentTable.formulaLineItems;
      this.preview.formulaLineItems = this.currentTable.formulaLineItems;
      this.editor.annotationRegionItems = this.currentTable.annotationRegionItems;
      this.preview.annotationRegionItems = this.currentTable.annotationRegionItems;
      this.currentTable.setFirstItemAsCurrent();
    });

    return tabs;
  }

  private initTable(formulaLineItems: DG.FormulaLine[], annotationRegionItems: DG.AnnotationRegion[]): Table {
    return new Table(formulaLineItems, annotationRegionItems,
      (itemIdx: number): boolean => {
        const isFormulaLine = this.preview.formulaLineItems.length > 0 && itemIdx < this.preview.formulaLineItems.length;
        itemIdx-= isFormulaLine ? 0 : this.preview.formulaLineItems.length;
        this.editor.update(itemIdx, isFormulaLine);
        return this.preview.update(itemIdx, isFormulaLine);
      }, this.preview.srcAxes, !this.showValueOnOpen);
  }

  private onOKAction() {
    this.host.save();
    this.creationControl.saveHistory();
  }

  private onItemCreatedAction(item: DG.FormulaLine, isFormulaLine: boolean = true): void {
    this.currentTable.add(item, isFormulaLine);
    this.editor.update(0, isFormulaLine);
    this.preview.update(0, isFormulaLine);
  }
}
