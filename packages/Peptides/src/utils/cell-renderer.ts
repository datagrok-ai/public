import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import * as type from './types';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {calculateMonomerPositionStatistics} from './algorithms';
import * as rxjs from 'rxjs';
import {showTooltipAt, TooltipOptions} from './tooltips';
import {MonomerPositionStats, PositionStats} from './statistics';
import {CLUSTER_TYPE} from '../viewers/logo-summary';
import {SARViewer} from '../viewers/sar-viewer';

export function renderCellSelection(canvasContext: CanvasRenderingContext2D, bound: DG.Rect): void {
  canvasContext.strokeStyle = DG.Color.toHtml(DG.Color.selectedRows);
  canvasContext.lineWidth = 3;
  canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
}

/** A function that sets amino acid residue cell renderer to the specified column */
export function setMonomerRenderer(col: DG.Column, alphabet: string): void {
  col.semType = C.SEM_TYPES.MONOMER;
  col.setTag(DG.TAGS.CELL_RENDERER, C.SEM_TYPES.MONOMER);
  col.setTag(C.TAGS.ALPHABET, alphabet);
}

export function renderMutationCliffCell(canvasContext: CanvasRenderingContext2D, currentMonomer: string,
  currentPosition: string, viewer: SARViewer, bound: DG.Rect): void {
  const positionStats = viewer.monomerPositionStats[currentPosition];
  const pVal = positionStats![currentMonomer]!.pValue;
  const currentMeanDifference = positionStats![currentMonomer]!.meanDifference;

  // Transform p-value to increase intensity for smaller values and decrease for larger values
  const maxPValComplement = 1 - positionStats!.general.maxPValue;
  const minPValComplement = 1 - positionStats!.general.minPValue;
  const pValCentering = Math.min(maxPValComplement, minPValComplement);
  const centeredMaxPValComplement = maxPValComplement - pValCentering;
  const centeredMinPValComplement = minPValComplement - pValCentering;
  const centeredPValLimit = Math.max(centeredMaxPValComplement, centeredMinPValComplement);
  const pValComplement = pVal === null ? 0 : 1 - pVal - pValCentering;

  const x = currentMeanDifference >= 0 ? pValComplement : -pValComplement;
  const coef = DG.Color.toHtml(pVal === null ? DG.Color.lightLightGray :
    DG.Color.scaleColor(x, -centeredPValLimit, centeredPValLimit));

  const halfWidth = bound.width / 2;
  const maxMeanDifference = Math.max(Math.abs(viewer.monomerPositionStats.general.minMeanDifference),
    viewer.monomerPositionStats.general.maxMeanDifference);
  const rCoef = Math.abs(currentMeanDifference) / maxMeanDifference;
  const maxRadius = 0.9 * halfWidth / 2; // Fill at most 90% of the half of the cell width
  const radius = Math.floor(maxRadius * rCoef);

  const midX = Math.ceil(bound.x + 1 + halfWidth);
  const midY = Math.ceil(bound.y + 1 + bound.height / 2);
  canvasContext.beginPath();
  canvasContext.fillStyle = coef;
  canvasContext.arc(midX - halfWidth / 2, midY, radius < 3 || pVal === null ? 3 : radius, 0, Math.PI * 2, true);
  canvasContext.closePath();
  canvasContext.fill();

  canvasContext.textBaseline = 'middle';
  canvasContext.textAlign = 'end';
  canvasContext.fillStyle = '#606060';
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.shadowBlur = 5;
  canvasContext.shadowColor = DG.Color.toHtml(DG.Color.white);
  const uniqueValues = new Set<number>();
  const substitutions = viewer.mutationCliffs?.get(currentMonomer)?.get(currentPosition)?.entries() ?? null;
  if (substitutions !== null) {
    for (const [key, value] of substitutions) {
      uniqueValues.add(key);
      for (const val of value)
        uniqueValues.add(val);
    }
  }
  if (uniqueValues.size !== 0)
    canvasContext.fillText(uniqueValues.size.toString(), midX + halfWidth - 5, midY, halfWidth - 5);

  const monomerSelection = viewer.mutationCliffsSelection[currentPosition];
  if (monomerSelection && monomerSelection.includes(currentMonomer))
    renderCellSelection(canvasContext, bound);
}

export function renderInvariantMapCell(canvasContext: CanvasRenderingContext2D, currentMonomer: string,
  currentPosition: string, invariantMapSelection: type.Selection, cellValue: number, bound: DG.Rect,
  color: number): void {
  //FIXME: This is a hack, because `color` value sometimes comes incomplete. E.g. we found that here `color` value is
  // 255 and its contrast color would be black, which is not visible on blue (color code) background. The full number
  // is actually 4278190335.
  color = DG.Color.fromHtml(DG.Color.toHtml(color));
  canvasContext.fillStyle = DG.Color.toHtml(color);
  canvasContext.fillRect(bound.x, bound.y, bound.width, bound.height);
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(color));
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  const monomerSelection = invariantMapSelection[currentPosition];
  if (monomerSelection && monomerSelection.includes(currentMonomer))
    renderCellSelection(canvasContext, bound);
}

export function renderLogoSummaryCell(canvasContext: CanvasRenderingContext2D, cellValue: string,
  clusterSelection: type.Selection, bound: DG.Rect): void {
  canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
  canvasContext.textAlign = 'center';
  canvasContext.textBaseline = 'middle';
  canvasContext.fillStyle = '#000';
  canvasContext.fillText(cellValue.toString(), bound.x + (bound.width / 2), bound.y + (bound.height / 2), bound.width);

  if (clusterSelection[CLUSTER_TYPE.CUSTOM].includes(cellValue) ||
    clusterSelection[CLUSTER_TYPE.ORIGINAL].includes(cellValue))
    renderCellSelection(canvasContext, bound);
}

export function drawLogoInBounds(ctx: CanvasRenderingContext2D, bounds: DG.Rect, stats: PositionStats, position: string,
  sortedOrder: string[], rowCount: number, cp: SeqPalette, monomerSelectionStats: {
    [monomer: string]: number
  } = {},
  drawOptions: type.DrawOptions = {}): { [monomer: string]: DG.Rect } {
  const pr = window.devicePixelRatio;
  drawOptions.symbolStyle ??= '16px Roboto, Roboto Local, sans-serif';
  drawOptions.upperLetterHeight ??= 12.2;
  drawOptions.upperLetterAscent ??= 0.25;
  drawOptions.marginVertical ??= 2;
  drawOptions.marginHorizontal ??= 2;
  drawOptions.selectionWidth ??= 1;
  drawOptions.textHeight ??= 13;
  drawOptions.headerStyle ??= `bold ${drawOptions.textHeight * pr}px Roboto, Roboto Local, sans-serif`;

  const totalSpace = (sortedOrder.length - 1) * drawOptions.upperLetterAscent; // Total space between letters
  const barHeight = (bounds.height - 2 * drawOptions.marginVertical - totalSpace - 1.25 * drawOptions.textHeight) * pr;
  const leftShift = drawOptions.marginHorizontal * 2;
  const barWidth = (bounds.width - (leftShift + drawOptions.marginHorizontal)) * pr;
  const xStart = (bounds.x + leftShift) * pr;
  const selectionWidth = Math.min(drawOptions.selectionWidth * pr, 1);
  const xSelection = (bounds.x + 1) * pr;
  let currentY = (bounds.y + drawOptions.marginVertical) * pr;

  const monomerBounds: { [monomer: string]: DG.Rect } = {};
  for (const monomer of sortedOrder) {
    const monomerHeight = barHeight * (stats[monomer]!.count / rowCount);
    const selectionHeight = barHeight * ((monomerSelectionStats[monomer] ?? 0) / rowCount);
    monomerBounds[monomer] = new DG.Rect(xStart / pr, currentY / pr, barWidth / pr, monomerHeight / pr);

    ctx.resetTransform();
    if (monomer !== '-' && monomer !== '') {
      const monomerTxt = monomerToShort(monomer, 5);
      const mTm: TextMetrics = ctx.measureText(monomerTxt);

      // Filling selection
      ctx.lineWidth = selectionWidth;
      ctx.line(xSelection, currentY, xSelection, currentY + selectionHeight, DG.Color.rowSelection);

      ctx.fillStyle = cp.get(monomer) ?? cp.get('other');
      ctx.textAlign = 'left';
      ctx.textBaseline = 'top';
      ctx.font = drawOptions.symbolStyle;
      // Hacks to scale uppercase characters to target rectangle
      const widthTransform = barWidth / mTm.width;
      const heightTransform = monomerHeight / drawOptions.upperLetterHeight;
      ctx.setTransform(widthTransform, 0, 0, heightTransform, xStart, currentY);
      ctx.fillText(monomerTxt, 0, 0);
    }
    currentY += monomerHeight + drawOptions.upperLetterAscent * pr;
  }

  // Drawing column header
  ctx.resetTransform();
  ctx.fillStyle = DG.Color.toHtml(DG.Color.black);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.font = drawOptions.headerStyle;
  ctx.fillText(position, (bounds.x + bounds.width / 2) * pr, (bounds.y + bounds.height - drawOptions.textHeight) * pr);

  return monomerBounds;
}

export type CellRendererOptions = {
  isSelectionTable?: boolean, headerSelectedMonomers?: type.SelectionStats, webLogoBounds: WebLogoBounds,
  cachedWebLogoTooltip: type.CachedWebLogoTooltip, colorPalette: SeqPalette, unhighlightCallback?: () => void,
  highlightCallback?: (mp: type.SelectionItem, dataFrame: DG.DataFrame, stats: MonomerPositionStats) => void,
  selectionCallback?: (monomerPosition: type.SelectionItem, options: type.SelectionOptions) => void,
};
export type WebLogoBounds = { [position: string]: { [monomer: string]: DG.Rect } };

export function setWebLogoRenderer(grid: DG.Grid, monomerPositionStats: MonomerPositionStats,
  positionColumns: DG.Column<string>[], activityCol: DG.Column<number>, options: CellRendererOptions,
  tooltipOptions: TooltipOptions = {
    x: 0, y: 0, mpStats: {} as MonomerPositionStats,
    monomerPosition: {} as type.SelectionItem,
  }): void {
  options.isSelectionTable ??= false;
  if (Object.keys(tooltipOptions.mpStats).length == 0)
    tooltipOptions.mpStats = monomerPositionStats;
  if (options.isSelectionTable && (!options.webLogoBounds || !options.cachedWebLogoTooltip)) {
    throw new Error('Peptides: Cannot set WebLogo renderer for selection table without `headerSelectedMonomers`, ' +
      '`webLogoBounds` and `cachedWebLogoTooltip` options.');
  }

  const df = grid.dataFrame;
  grid.setOptions({'colHeaderHeight': 130});
  const headerRenderer = (gcArgs: DG.GridCellRenderArgs): void => {
    const ctx = gcArgs.g;
    const bounds = gcArgs.bounds;
    const col = gcArgs.cell.tableColumn;

    ctx.save();
    try {
      ctx.beginPath();
      ctx.rect(bounds.x, bounds.y, bounds.width, bounds.height);
      ctx.clip();

      //TODO: optimize
      if (gcArgs.cell.isColHeader && col?.semType === C.SEM_TYPES.MONOMER) {
        const isDfFiltered = df.filter.anyFalse;
        const stats = (isDfFiltered ?
          calculateMonomerPositionStatistics(activityCol, df.filter, positionColumns, {
            isFiltered: true,
            columns: [col.name],
          }) :
          monomerPositionStats)[col.name];
        if (!stats)
          return;
        //TODO: precalc on stats creation
        const sortedStatsOrder = Object.keys(stats).sort((a, b) => {
          if (a === '' || a === '-')
            return +1;
          else if (b === '' || b === '-')
            return -1;
          return 0;
        }).filter((v) => v !== 'general');

        options.webLogoBounds![col.name] = drawLogoInBounds(ctx, bounds, stats, col.name,
          sortedStatsOrder, df.filter.trueCount, options.colorPalette,
          options.headerSelectedMonomers ? options.headerSelectedMonomers[col.name] : {});
        gcArgs.preventDefault();
      }
    } catch (e) {
      console.warn(`PeptidesHeaderLogoError: couldn't render WebLogo for column \`${col!.name}\`. ` +
        `See original error below.`);
      console.warn(e);
    } finally {
      ctx.restore();
    }
  };
  grid.onCellRender.subscribe((gcArgs) => headerRenderer(gcArgs));

  const eventAction = (ev: MouseEvent): void => {
    const cell = grid.hitTest(ev.offsetX, ev.offsetY);
    if (cell?.isColHeader && cell.tableColumn?.semType === C.SEM_TYPES.MONOMER) {
      const monomerPosition = findWebLogoMonomerPosition(cell, ev, (options.webLogoBounds));
      if (monomerPosition === null) {
        if (!options.isSelectionTable && options.unhighlightCallback != null)
          options.unhighlightCallback();
        return;
      }
      tooltipOptions.monomerPosition = monomerPosition;
      requestBarchartAction(ev, monomerPosition, df, activityCol, options, tooltipOptions);
      if (!options.isSelectionTable && options.highlightCallback != null)
        options.highlightCallback(monomerPosition, df, monomerPositionStats);
    }
  };

  // The following events makes the barchart interactive
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'mousemove').subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  rxjs.fromEvent<MouseEvent>(grid.overlay, 'click').subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
}

function requestBarchartAction(ev: MouseEvent, monomerPosition: type.SelectionItem, df: DG.DataFrame,
  activityCol: DG.Column<number>, options: CellRendererOptions, tooltipOptions: TooltipOptions): void {
  if (ev.type === 'click' && !options.isSelectionTable && options.selectionCallback != null)
    options.selectionCallback(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
  else {
    const bar = `${monomerPosition.positionOrClusterType} = ${monomerPosition.monomerOrCluster}`;
    if (options.cachedWebLogoTooltip!.bar === bar)
      ui.tooltip.show(options.cachedWebLogoTooltip!.tooltip!, ev.clientX, ev.clientY);
    else {
      options.cachedWebLogoTooltip!.bar = bar;
      tooltipOptions.x = ev.clientX;
      tooltipOptions.y = ev.clientY;
      tooltipOptions.monomerPosition = monomerPosition;
      options.cachedWebLogoTooltip!.tooltip = showTooltipAt(df, activityCol, {}, tooltipOptions);
    }
  }
}

function findWebLogoMonomerPosition(cell: DG.GridCell, ev: MouseEvent, webLogoBounds: WebLogoBounds,
): type.SelectionItem | null {
  const barCoords = webLogoBounds[cell.tableColumn!.name];
  for (const [monomer, coords] of Object.entries(barCoords)) {
    const isIntersectingX = ev.offsetX >= coords.x && ev.offsetX <= coords.x + coords.width;
    const isIntersectingY = ev.offsetY >= coords.y && ev.offsetY <= coords.y + coords.height;
    if (isIntersectingX && isIntersectingY)
      return {monomerOrCluster: monomer, positionOrClusterType: cell.tableColumn!.name};
  }

  return null;
}
