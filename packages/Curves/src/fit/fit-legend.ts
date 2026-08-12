import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {IFitChartData, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {ColorType, getSeriesColor} from './render-utils';
import {isLegendShown} from './fit-layout';

/** Whether the curves get named at all: the plot has to have room for it, and it can be turned off. */
export function isLegendVisible(data: IFitChartData, screenBounds: DG.Rect): boolean {
  return (data.chartOptions?.showLegend ?? true) && isLegendShown(screenBounds);
}

/** A legend row: one curve, or the name of a column that holds several of them. */
export interface LegendEntry {
  text: string;
  color: string;
  column?: string;
  pointColor?: string;
  marker?: DG.MARKER_TYPE;
}

interface LegendRow {
  entry: LegendEntry;
  /** What is drawn, shortened to the width the legend allows itself. */
  shown: string;
  box: DG.Rect;
}

const LEGEND_FONT = '11px Roboto, "Roboto Local"';
/** The legend annotates the plot, so it may not take more than a corner of it. */
const MAX_LEGEND_WIDTH_FRACTION = 0.55;
const MAX_LEGEND_HEIGHT_FRACTION = 0.5;
const LEGEND_PX_PADDING = 3;
const CORNER_PX_INSET = 4;
/** Ink this close to a corner still crowds it, which is what tells a free corner from a tight one. */
const CORNER_PX_CLEARANCE = 14;
/** Below this a corner counts as free enough, so a stray point cannot move the legend. */
const CORNER_FREE_ENOUGH = 2;
/** Once a name is shown whole it takes this much crowding to start shortening it, so that resizing
 * a cell by a few pixels does not flip the legend between the two. */
const NAME_STICKY_CROWDING = 12;
const LEGEND_BACKGROUND = 'rgba(255, 255, 255, 0.75)';
/** The hover may fill most of the window; a tooltip taller than that is off the screen anyway. */
const TOOLTIP_HEIGHT_FRACTION = 0.75;
/** A row of a listed legend is a label; the whole of a name is what hovering that one row is for. */
const MAX_TOOLTIP_ROW_PX = 400;

/** How many rows the hover lists before it counts the rest instead. */
export function maxTooltipRows(): number {
  return Math.max(10, Math.floor(window.innerHeight * TOOLTIP_HEIGHT_FRACTION /
    FitConstants.LEGEND_RECORD_PX_HEIGHT));
}

let measureCanvas: CanvasRenderingContext2D | undefined;
const lastLayout = new WeakMap<IFitChartData, {box: DG.Rect, rows: LegendRow[], moreBox?: DG.Rect,
  entries: LegendEntry[], corner: number, whole: boolean}>();

/** A context of its own, so a row can be measured while hit testing, which has no canvas to draw on. */
function measure(text: string): number {
  measureCanvas ??= document.createElement('canvas').getContext('2d')!;
  measureCanvas.font = LEGEND_FONT;
  return measureCanvas.measureText(text).width;
}

function ellipsize(text: string, maxWidth: number): string {
  if (measure(text) <= maxWidth)
    return text;
  let trimmed = text;
  while (trimmed.length > 1 && measure(`${trimmed}…`) > maxWidth)
    trimmed = trimmed.slice(0, -1);
  return `${trimmed}…`;
}

/** The marker a series draws its points with, which is what makes a legend row recognizable. */
function seriesMarker(series: IFitSeries): DG.MARKER_TYPE {
  const first = series.points?.length ? series.points[0].marker : null;
  const shared = series.points?.every((p) => p.marker === first && p.marker != null) ? first : null;
  return (shared ?? series.markerType ?? DG.MARKER_TYPE.CIRCLE) as DG.MARKER_TYPE;
}

function legendName(series: IFitSeries, useAuxLegendNames?: boolean): string | undefined {
  const name = useAuxLegendNames ? series.auxLegendName ?? series.name : series.name;
  return name === '' || name === null ? undefined : name;
}

/** The rows the legend would draw, in order. Colors are taken by the series' own index, so a row
 * matches the curve it names. */
export function legendEntries(data: IFitChartData): LegendEntry[] {
  const series = data.series ?? [];
  const columns = [...new Set(series.map((s) => s.columnName))].filter((name) => name !== null && name !== undefined);
  const entries: LegendEntry[] = [];
  for (const column of columns) {
    const inColumn = series.map((s, idx) => ({series: s, idx}))
      .filter(({series}) => series.columnName === column && legendName(series, data.chartOptions?.useAuxLegendNames));
    if (inColumn.length === 0)
      continue;
    // a column holding a single curve says nothing that the curve's own row does not
    if (data.chartOptions?.showColumnLabel && inColumn.length > 1)
      entries.push({text: column!, color: 'black'});
    for (const {series, idx} of inColumn) {
      entries.push({text: legendName(series, data.chartOptions?.useAuxLegendNames)!, column: column,
        color: getSeriesColor(series, idx, ColorType.FIT_LINE),
        pointColor: getSeriesColor(series, idx, ColorType.POINT), marker: seriesMarker(series)});
    }
  }

  // one name in two columns identifies nothing on its own; the same name twice in one column would
  // only get the same prefix twice, so it is left as the data has it
  const columnsByName = new Map<string, Set<string>>();
  for (const entry of entries.filter((e) => e.marker !== undefined && e.column))
    (columnsByName.get(entry.text) ?? columnsByName.set(entry.text, new Set()).get(entry.text)!).add(entry.column!);
  for (const entry of entries) {
    if (entry.marker !== undefined && entry.column && columnsByName.get(entry.text)!.size > 1)
      entry.text = `${entry.column}: ${entry.text}`;
  }
  return entries;
}

/** How far a point is from a box, 0 when inside it. */
function distanceTo(box: DG.Rect, p: DG.Point): number {
  const dx = Math.max(box.x - p.x, 0, p.x - box.maxX);
  const dy = Math.max(box.y - p.y, 0, p.y - box.maxY);
  return Math.sqrt(dx * dx + dy * dy);
}

/** How much ink a corner has to live with: what falls inside it, or close enough to crowd it. */
export function cornerCrowding(corner: DG.Rect, drawnAt?: DG.Point[]): number {
  if (!drawnAt?.length)
    return 0;
  const near = corner.inflate(CORNER_PX_CLEARANCE, CORNER_PX_CLEARANCE);
  return drawnAt.filter((p) => near.contains(p.x, p.y)).length;
}

/** The corner the legend takes. The top right one unless another is clearly emptier, so that the
 * legend does not hop about as rows are hovered or the grid scrolls. */
export function chooseLegendCorner(dataBox: DG.Rect, width: number, height: number,
  drawnAt?: DG.Point[], previous?: number): {rect: DG.Rect, index: number} {
  const corners = [
    new DG.Rect(dataBox.maxX - width - CORNER_PX_INSET, dataBox.y + CORNER_PX_INSET, width, height),
    new DG.Rect(dataBox.x + CORNER_PX_INSET, dataBox.y + CORNER_PX_INSET, width, height),
    new DG.Rect(dataBox.maxX - width - CORNER_PX_INSET, dataBox.maxY - height - CORNER_PX_INSET, width, height),
    new DG.Rect(dataBox.x + CORNER_PX_INSET, dataBox.maxY - height - CORNER_PX_INSET, width, height),
  ];
  // where it already is, or the top right one to begin with
  const held = previous ?? 0;
  if (!drawnAt?.length)
    return {rect: corners[held], index: held};
  const room = corners.map((corner) => ({taken: cornerCrowding(corner, drawnAt),
    clearance: Math.min(...drawnAt.map((p) => distanceTo(corner, p)))}));
  // among corners the curves crowd equally, the one they keep furthest away
  let best = 0;
  for (let i = 1; i < corners.length; i++) {
    if (room[i].taken < room[best].taken ||
      (room[i].taken === room[best].taken && room[i].clearance > room[best].clearance))
      best = i;
  }
  return room[held].taken <= CORNER_FREE_ENOUGH || room[best].taken * 2 >= room[held].taken ?
    {rect: corners[held], index: held} : {rect: corners[best], index: best};
}

/** Where each row goes, in the corner the curves leave freest. */
function legendLayout(data: IFitChartData, dataBox: DG.Rect, ratio: number,
  drawnAt?: DG.Point[]): {rows: LegendRow[], hidden: number, moreBox?: DG.Rect, entries: LegendEntry[],
  corner: number, whole: boolean} {
  const previous = lastLayout.get(data);
  const entries = legendEntries(data);
  if (entries.length === 0)
    return {rows: [], hidden: 0, entries: entries, corner: 0, whole: true};

  const glyph = glyphWidth(ratio);
  if (dataBox.width * MAX_LEGEND_WIDTH_FRACTION - glyph <= 0)
    return {rows: [], hidden: 0, entries: entries, corner: 0, whole: true};

  const fits = Math.floor(dataBox.height * MAX_LEGEND_HEIGHT_FRACTION / FitConstants.LEGEND_RECORD_PX_HEIGHT);
  // the last row goes to saying how much was left out, rather than to one more curve out of many
  const drawn = entries.length > fits ? entries.slice(0, Math.max(1, fits - 1)) : entries;
  const hidden = entries.length - drawn.length;
  const height = FitConstants.LEGEND_RECORD_PX_HEIGHT * (drawn.length + (hidden > 0 ? 1 : 0));

  // a name is worth showing whole when the corner it would go in has the room to spare; the cap is
  // for the charts where it does not
  const natural = Math.min(glyph + Math.max(...drawn.map((entry) => measure(entry.text)),
    hidden > 0 ? measure(`+${hidden} more`) : 0), dataBox.width - 2 * CORNER_PX_INSET);
  const naturalCorner = chooseLegendCorner(dataBox, natural, height, drawnAt, previous?.corner);
  // what it took to shorten the name last time, so a few pixels of resizing do not flip it back
  const endures = previous?.whole ? NAME_STICKY_CROWDING : CORNER_FREE_ENOUGH;
  const whole = cornerCrowding(naturalCorner.rect, drawnAt) <= endures;
  const allowed = whole ? natural : Math.min(natural, dataBox.width * MAX_LEGEND_WIDTH_FRACTION);

  const shown = drawn.map((entry) => ellipsize(entry.text, allowed - glyph));
  const width = glyph + Math.max(...shown.map(measure));
  const pick = whole ? naturalCorner : chooseLegendCorner(dataBox, width, height, drawnAt, previous?.corner);
  const block = pick.rect;
  const rows: LegendRow[] = drawn.map((entry, i) => ({entry: entry, shown: shown[i],
    box: new DG.Rect(block.x, block.y + FitConstants.LEGEND_RECORD_PX_HEIGHT * i,
      width, FitConstants.LEGEND_RECORD_PX_HEIGHT)}));
  const moreBox = hidden === 0 ? undefined : new DG.Rect(block.x,
    block.y + FitConstants.LEGEND_RECORD_PX_HEIGHT * rows.length, width, FitConstants.LEGEND_RECORD_PX_HEIGHT);
  return {rows: rows, hidden: hidden, moreBox: moreBox, entries: entries, corner: pick.index, whole: whole};
}

function glyphWidth(ratio: number): number {
  return Math.floor(FitConstants.POINT_PX_SIZE * ratio * 1.5) + FitConstants.LEGEND_RECORD_LINE_PX_WIDTH +
    FitConstants.LEGEND_RECORD_LINE_RIGHT_PX_MARGIN;
}

/** The line and the marker that say which curve a row names. */
function drawGlyph(g: CanvasRenderingContext2D, entry: LegendEntry, x: number, middle: number, ratio: number): void {
  if (entry.marker === undefined)
    return;
  const markerSize = Math.floor(FitConstants.POINT_PX_SIZE * ratio * 1.5);
  g.beginPath();
  g.strokeStyle = entry.color;
  g.lineWidth = 2 * ratio;
  g.moveTo(x + markerSize, middle);
  g.lineTo(x + markerSize + FitConstants.LEGEND_RECORD_LINE_PX_WIDTH, middle);
  g.stroke();
  DG.Paint.marker(g, entry.marker, x + markerSize / 2, middle, entry.pointColor!, markerSize);
}

/** The legend as a tooltip: the same glyphs and names, for what the plot had no room to say. */
export function legendTooltipElement(entries: LegendEntry[]): HTMLElement {
  // hovering one row is there to give back the whole of a name the plot shortened, so it is said in
  // full however long; a list of many is a legend, and a legend's rows are labels
  const whole = entries.length === 1;
  const shown = entries.slice(0, maxTooltipRows());
  const rows = shown.map((entry) => {
    const canvas = ui.canvas(glyphWidth(1), FitConstants.LEGEND_RECORD_PX_HEIGHT);
    drawGlyph(canvas.getContext('2d')!, entry, 0, FitConstants.LEGEND_RECORD_PX_HEIGHT / 2, 1);
    const row = ui.divH([canvas, ui.divText(whole ? entry.text : ellipsize(entry.text, MAX_TOOLTIP_ROW_PX),
      {style: {color: entry.color}})]);
    row.style.alignItems = 'center';
    return row;
  });
  if (entries.length > shown.length)
    rows.push(ui.divText(`+${entries.length - shown.length} more`));
  return ui.divV(rows);
}

/** Draws the legend in the corner of the plot, over a backdrop so the curves underneath
 * do not make it unreadable. */
export function renderLegend(g: CanvasRenderingContext2D, data: IFitChartData, dataBox: DG.Rect, ratio: number,
  drawnAt?: DG.Point[]): void {
  const layout = legendLayout(data, dataBox, ratio, drawnAt);
  const {rows, hidden} = layout;
  if (rows.length === 0)
    return;
  // hovering happens long after the chart was drawn, and only the render knows which corner it chose
  lastLayout.set(data, {box: dataBox, rows: rows, moreBox: layout.moreBox, entries: layout.entries,
    corner: layout.corner, whole: layout.whole});

  g.save();
  g.font = LEGEND_FONT;
  // the axes labels rendered just before leave textAlign centred
  g.textAlign = 'left';
  const block = rows[0].box;
  const height = FitConstants.LEGEND_RECORD_PX_HEIGHT * (rows.length + (hidden > 0 ? 1 : 0));
  g.fillStyle = LEGEND_BACKGROUND;
  g.fillRect(block.x - LEGEND_PX_PADDING, block.y - LEGEND_PX_PADDING,
    block.width + 2 * LEGEND_PX_PADDING, height + 2 * LEGEND_PX_PADDING);

  for (const row of rows) {
    const baseline = row.box.maxY - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN;
    drawGlyph(g, row.entry, row.box.x, baseline - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN, ratio);
    g.fillStyle = row.entry.color;
    g.fillText(row.shown, row.box.x + (row.entry.marker === undefined ? 0 : glyphWidth(ratio)), baseline);
  }

  if (hidden > 0) {
    g.fillStyle = 'black';
    g.fillText(`+${hidden} more`, block.x, rows[rows.length - 1].box.maxY +
      FitConstants.LEGEND_RECORD_PX_HEIGHT - FitConstants.LEGEND_RECORD_LINE_BOTTOM_PX_MARGIN);
  }
  g.restore();
}

/** What the legend could not say at that point: the whole of a shortened name, or the whole legend
 * when the pointer is on the row counting what did not fit. */
export function legendTooltip(data: IFitChartData, dataBox: DG.Rect, x: number,
  y: number): LegendEntry[] | null {
  const last = lastLayout.get(data);
  // the same data is drawn in the grid and in the enlarged chart, so the rows only apply to the plot
  // they were laid out for
  if (!last || Math.round(last.box.x) !== Math.round(dataBox.x) ||
    Math.round(last.box.width) !== Math.round(dataBox.width))
    return null;
  if (last.moreBox?.contains(x, y))
    return last.entries;
  const row = last.rows.find((r) => r.box.contains(x, y));
  return row && row.shown !== row.entry.text ? [row.entry] : null;
}
