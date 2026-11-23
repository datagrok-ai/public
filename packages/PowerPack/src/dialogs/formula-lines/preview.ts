import * as DG from 'datagrok-api/dg';

import {formatAreaFormula, validateBandFormula} from './formulas';
import {ItemList} from './item-list';
import {createPreviewViewer} from './preview-viewer';
import {AxisColumns, AxisNames, EditorItem, isAnnotationRegionType, ITEM_ORIENTATION, ITEM_TYPE} from './types';
import {axesToOptions, previewDataFrame, readAxisColumns, sourceAxes} from './viewer-axes';

export type DrawingViewer = DG.ScatterPlotViewer | DG.LineChartViewer | DG.HistogramViewer | DG.BarChartViewer |
  DG.BoxPlot | DG.DensityPlotViewer;

/** True for the viewers an annotation region can be drawn on. */
export const supportsRegionDrawing = (viewer: DG.Viewer): viewer is DrawingViewer =>
  viewer instanceof DG.ScatterPlotViewer || viewer instanceof DG.LineChartViewer ||
  viewer instanceof DG.HistogramViewer || viewer instanceof DG.BarChartViewer ||
  viewer instanceof DG.BoxPlot || viewer instanceof DG.DensityPlotViewer;

/** How much the other formula lines are dimmed next to the line being edited. */
const GHOST_ALPHA = 0.4;

/** Returns `color` with its alpha scaled by `alpha` (0..1), or undefined if it can't be parsed. */
function dimColor(color: string | undefined | null, alpha: number): string | undefined {
  if (!color)
    return undefined;
  try {
    const c = DG.Color.fromHtml(color.trim());
    if (c == null || Number.isNaN(c))
      return undefined;
    const a = DG.Color.a(c);
    return DG.Color.toHtml(DG.Color.setAlpha(c, Math.min(a, Math.max(8, Math.round(a * alpha)))));
  } catch {
    return undefined;
  }
}

/**
 * Preview Helper for Formula Lines.
 * Scatter Plot viewer by default.
 */
export class Preview {
  public readonly viewer: DG.Viewer;
  public readonly dataFrame: DG.DataFrame;

  /** Lines from the inactive tab; rendered as ghosts only */
  public otherFormulaLines: DG.FormulaLine[] = [];

  /** Source viewer axes */
  public readonly srcAxes?: AxisNames;

  // Text of the item on the preview, reported in the viewer status
  private currentItem = '';

  public get root(): HTMLElement {return this.viewer.root;}

  /** Returns the current columns pair of the preview viewer. */
  public get axisCols(): AxisColumns {return readAxisColumns(this.viewer, this.dataFrame);}

  constructor(
    public items: ItemList,
    src: DG.DataFrame | DG.Viewer,
    onContextMenu: (valY?: number, valX?: number) => void,
  ) {
    if (src instanceof DG.Viewer) {
      this.dataFrame = previewDataFrame(src);
      this.srcAxes = sourceAxes(src, this.dataFrame);
    } else
      this.dataFrame = src;

    this.viewer = createPreviewViewer(src, this.dataFrame);

    if (this.srcAxes)
      this.setAxes(this.srcAxes);

    this.viewer.addStatusProvider('formula-lines-preview', () => ({values: {'current item': this.currentItem}}));

    /**
     * Creates special context menu for preview Scatter Plot.
     * Before opening the menu, it calculates the world coordinates of the click point.
     */
    this.viewer.root.addEventListener('contextmenu', (event: MouseEvent) => {
      event.preventDefault();
      if (this.viewer instanceof DG.ScatterPlotViewer || this.viewer instanceof DG.LineChartViewer) {
        const worldPoint = this.viewer.screenToWorld(event.offsetX, event.offsetY);
        onContextMenu(worldPoint.y, worldPoint.x);
      } else if (supportsRegionDrawing(this.viewer))
        onContextMenu();
    });
  }

  /**
   * Shows the item of the [rowIdx] row on the preview viewer.
   * Returns true if the rendering was successful, false otherwise.
   */
  public update(rowIdx: number): boolean {
    /** If there are no lines, try to set the axes as in the original viewer. */
    if (rowIdx < 0 && this.srcAxes)
      this.setAxes(this.srcAxes);

    try {
      const item = this.items.at(rowIdx)!;
      return this.items.isFormulaLine(rowIdx) ? this.showFormulaLine(item) : this.showAnnotationRegion(item);
    } catch {
      this.clear();
      return false;
    }
  }

  private showFormulaLine(item: DG.FormulaLine): boolean {
    if (item.type === ITEM_TYPE.BAND && validateBandFormula(item.formula ?? '') !== '') {
      this.clear();
      return false;
    }

    /** Duplicate the original item to display it even if it's hidden */
    const previewItem = structuredClone(item);
    previewItem.visible = true;

    this.clear();
    // Render the other formula lines as ghosts so the user keeps context for the line being edited.
    // Bands and regions are not ghosted; ghosts are only added when the current item is a line.
    if (item.type !== ITEM_TYPE.BAND) {
      for (const line of [...this.items.formulaLines, ...this.otherFormulaLines]) {
        if (line !== item)
          this.addLineGhost(line);
      }
    }
    this.viewer.meta.formulaLines.add(previewItem);
    this.setAxes(this.getItemAxes(previewItem));
    this.currentItem = previewItem.formula ?? '';
    return true;
  }

  private showAnnotationRegion(item: DG.AnnotationRegion): boolean {
    /** Duplicate the original item to display it even if it's hidden */
    const previewItem = structuredClone(item);
    previewItem.hidden = false;

    this.clear();
    this.viewer.meta.annotationRegions.add(previewItem);
    this.setAxes(this.getItemAxes(previewItem));
    this.currentItem = formatAreaFormula(previewItem);
    return true;
  }

  private addLineGhost(line: DG.FormulaLine | undefined): void {
    if (!line || !line.visible || line.type === ITEM_TYPE.BAND)
      return;
    try {
      const g = DG.FormulaLinesHelper.setDefaults(structuredClone(line));
      const dimmed = dimColor(g.color, GHOST_ALPHA);
      if (!dimmed)
        return;
      g.color = dimmed;
      this.viewer.meta.formulaLines.add(g);
    } catch {}
  }

  private clear(): void {
    this.viewer.meta.annotationRegions.clear();
    this.viewer.meta.formulaLines.clear();
    this.currentItem = '';
  }

  private setAxes(names: AxisNames): void {
    this.viewer.setOptions(axesToOptions(this.viewer, this.dataFrame, names));
  }

  // The axes named in the formula, adjusted to the axes of the original viewer where possible
  private getItemAxes(axesItem: EditorItem): AxisNames {
    const result = isAnnotationRegionType(axesItem.type) ?
      this.getRegionAxes(axesItem) :
      this.getFormulaLineAxes(axesItem);
    if (axesItem.type === ITEM_TYPE.AREA_REGION_ANNOTATION)
      return result;

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

        const tmpMap = result.xMap;
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

  private getRegionAxes(region: DG.AnnotationRegion): AxisNames {
    if (region.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
      const area = region as DG.AreaAnnotationRegion;
      return {y: area.y, x: area.x, yMap: area.yMap, xMap: area.xMap};
    }

    const item = region as DG.FormulaAnnotationRegion;
    const meta1 = item.formula1 ?
      DG.FormulaLinesHelper.getMetaByFormula(item.formula1, ITEM_TYPE.LINE) :
      null;

    const meta2 = item.formula2 ?
      DG.FormulaLinesHelper.getMetaByFormula(item.formula2, ITEM_TYPE.LINE) :
      null;

    // For `axisCol = const` formulas argName is undefined, so the srcAxes swap decides which axis funcName belongs to.
    return {
      y: meta1?.funcName ?? meta2?.funcName,
      x: meta1?.argName ?? meta2?.argName,
      yMap: item.yMap,
      xMap: item.xMap,
    };
  }

  private getFormulaLineAxes(item: DG.FormulaLine): AxisNames {
    const itemMeta = DG.FormulaLinesHelper.getMeta(item);
    return {
      y: item.orientation === ITEM_ORIENTATION.VERTICAL ? itemMeta.argName : itemMeta.funcName,
      x: item.orientation === ITEM_ORIENTATION.VERTICAL ? itemMeta.funcName : itemMeta.argName,
      xMap: item.xMap,
      yMap: item.yMap,
    };
  }
}
