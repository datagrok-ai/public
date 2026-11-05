/* What the GPU renderers draw: the core scatter plot's rows, positions, colors, sizes and shapes, read once per frame
 * and rebuilt only when their inputs change. A renderer receives the arrays through its sink and draws a frame from
 * the returned description; it never reads the viewer itself, so both backends draw the same thing. */
import * as DG from 'datagrok-api/dg';
import {isDotRenderer, MARKER_TYPES, plotCanvases, readLook, ScatterLook} from './scatterplot-look';
const CIRCLE_BORDER = MARKER_TYPES.indexOf('circle border');
/** The size the core draws a row with an empty Size value at, when *Show Markers With Empty Size* is on. */
const EMPTY_SIZE_MARKER_SIZE = 6;
/** Positions are float32 offsets from an origin near the viewport; once the viewport is this many of its widths
 * away from the origin, they are rebuilt around it, so that zooming deep in keeps sub-pixel precision. */
const REBASE_DISTANCE = 100;
const MIN_LOG = 1e-30;
const NO_ROWS = new Uint32Array(0);

/** The rows of each pass: filtered-out rows, filtered rows, and selected rows once more, on top. */
export type RowList = 'filtered' | 'filteredOut' | 'selected';

/** Receives the arrays that changed; a renderer copies them to the GPU and keeps no reference. */
export interface ScatterSink {
  setPositions(axis: 'x' | 'y', values: Float32Array): void;
  setRowValues(kind: 'colors' | 'sizes' | 'types', values: Uint32Array | Float32Array): void;
  /** The selection as bitset words, least significant bit first. */
  setSelection(bits: Int32Array): void;
  setRows(list: RowList, rows: Uint32Array): void;
}

/** One frame of the plot. Row values come from the arrays last set on the sink, or are one value for every row. */
export interface Frame {
  /** The core's pixel-buffer dot renderer: one opaque device pixel per row, no shapes, no selected rows on top. */
  dots: boolean;
  dpr: number;
  viewBox: DG.Rect;
  /** A row at (x, y) is at (xform[0] + xform[1] * x, xform[2] + xform[3] * y) on the viewer canvas, in CSS pixels. */
  xform: [number, number, number, number];
  /** The GPU canvas covers the plot area, snapped to the device pixels of the viewer canvas. */
  origin: [number, number];
  size: [number, number];
  rowColors: boolean;
  /** ARGB, as the core gives colors. */
  color: number;
  rowSizes: boolean;
  markerSize: number;
  rowTypes: boolean;
  markerType: number;
  /** Selected rows take [selectedColor], as the core paints them. */
  selection: boolean;
  selectedColor: number;
  /** The color of every filtered-out row, or null when they keep their row colors. */
  filteredOutColor: number | null;
  opacity: number;
  /** The width of the darker border the core strokes around markers; 0 for none. */
  border: number;
  /** `Marker.draw` caches the markers of a canvas this large and blits them at whole CSS pixels. */
  markerCache: boolean;
  largestMarker: number;
  rows: Record<RowList, number>;
}

const objectIds = new WeakMap<object, number>();
let lastObjectId = 0;

/** Tells apart columns and bitsets that share a name and a version: versions start over on every new object. */
function objectId(x: {dart: object} | null | undefined): number {
  if (x == null)
    return 0;
  let id = objectIds.get(x.dart);
  if (id === undefined)
    objectIds.set(x.dart, id = ++lastObjectId);
  return id;
}

/** The raw value of an empty cell; `col.getRawData()` of a datetime column is float64, of a double column float32. */
function noneValue(col: DG.Column, raw: ArrayLike<number>): number {
  if (col.type == DG.COLUMN_TYPE.INT)
    return DG.INT_NULL;
  return raw instanceof Float32Array ? Math.fround(DG.FLOAT_NULL) : DG.FLOAT_NULL;
}

export function largestMarker(sc: DG.ScatterPlotViewer,
  look: Pick<ScatterLook, 'sizeColumnName' | 'markerMaxSize'>): number {
  return Math.max(look.sizeColumnName ? look.markerMaxSize : sc.getMarkerSize(0), EMPTY_SIZE_MARKER_SIZE);
}


/** One axis: row positions as float32 offsets from [origin], in log space on a log axis. */
class Axis {
  private key = '';
  private origin = 0;
  private log = false;
  /** 1 where the row has no position: empty, infinite, or not positive on a log axis. */
  invalid = new Uint8Array(0);
  version = 0;

  constructor(private readonly name: 'x' | 'y') {}

  /** Rebuilds the positions when the column changed, or when the viewport [min, max] drifted far from the origin. */
  update(sink: ScatterSink, col: DG.Column, log: boolean, min: number, max: number): void {
    const categorical = col.type == DG.COLUMN_TYPE.STRING;
    const key = `${objectId(col)}|${col.version}|${col.length}|${log}`;
    this.log = log;
    const [lo, hi] = this.range(min, max);
    const center = (lo + hi) / 2;
    if (key == this.key && (categorical || !(Math.abs(center - this.origin) > REBASE_DISTANCE * (hi - lo))))
      return;

    this.key = key;
    this.origin = categorical || !Number.isFinite(center) ? 0 : center;
    this.version++;
    const n = col.length;
    const raw = col.getRawData();
    const values = new Float32Array(n);
    const invalid = this.invalid.length == n ? this.invalid : new Uint8Array(n);
    if (categorical) {
      for (let i = 0; i < n; i++)
        values[i] = raw[i];
      invalid.fill(0);
    }
    else {
      const none = noneValue(col, raw);
      for (let i = 0; i < n; i++) {
        const v = raw[i] as number;
        const bad = v === none || !Number.isFinite(v) || (log && v <= 0);
        invalid[i] = bad ? 1 : 0;
        values[i] = bad ? 0 : (log ? Math.log(Math.max(v, MIN_LOG)) : v) - this.origin;
      }
    }
    this.invalid = invalid;
    sink.setPositions(this.name, values);
  }

  /** [start + scale * position] places a row on the screen as the core's `Transform.worldToScreen` does. */
  transform(min: number, max: number, start: number, length: number, inverse: boolean): [number, number] {
    const [lo, hi] = this.range(min, max);
    const scale = length / (hi - lo);
    const offset = (this.origin - lo) * scale;
    return inverse ? [start + length - offset, -scale] : [start + offset, scale];
  }

  private range(min: number, max: number): [number, number] {
    return this.log ? [Math.log(Math.max(min, MIN_LOG)), Math.log(max)] : [min, max];
  }
}


export class ScatterData {
  private readonly x = new Axis('x');
  private readonly y = new Axis('y');
  private colorsKey = '';
  /** The selection when the row colors were taken: the core painted these rows in the selection color. */
  private colorsSelection: Int32Array | null = null;
  private sizesKey = '';
  private typesKey = '';
  private selectionKey = '';
  /** Size column state: 1 where the size is empty, 2 where it is not positive on a log scale. */
  private sizeFlags = new Uint8Array(0);
  private sizeFlagsKey = '';
  private sizeFlagsVersion = 0;
  private emptySizes = false;
  private filteredKey = '';
  private filtered: Uint32Array = NO_ROWS;
  private filteredVersion = 0;
  private filteredOutCount = 0;
  private selectedKey = '';
  private selectedCount = 0;
  private dataFrame: DG.DataFrame | null = null;
  private metadataSub: {unsubscribe(): void} | null = null;
  private metadataVersion = 0;

  constructor(private readonly sink: ScatterSink) {}

  dispose(): void {
    this.metadataSub?.unsubscribe();
  }

  update(sc: DG.ScatterPlotViewer): Frame {
    const look = readLook(sc);
    const df = sc.dataFrame;
    const dpr = window.devicePixelRatio;
    const viewBox = sc.viewBox;
    const viewport = sc.viewport;
    this.watchMetadata(df);

    const xCol = df.col(look.xColumnName)!;
    const yCol = df.col(look.yColumnName)!;
    const xLog = xCol.isNumerical && look.xAxisType == DG.AxisType.logarithmic;
    const yLog = yCol.isNumerical && look.yAxisType == DG.AxisType.logarithmic;
    this.x.update(this.sink, xCol, xLog, viewport.x, viewport.x + viewport.width);
    this.y.update(this.sink, yCol, yLog, viewport.y, viewport.y + viewport.height);
    const xform = [
      ...this.x.transform(viewport.x, viewport.x + viewport.width, viewBox.x, viewBox.width, look.invertXAxis),
      // the core's Y axis grows upwards unless inverted
      ...this.y.transform(viewport.y, viewport.y + viewport.height, viewBox.y, viewBox.height, !look.invertYAxis),
    ] as [number, number, number, number];

    const dots = isDotRenderer(look);
    const sizeCol = look.sizeColumnName ? df.col(look.sizeColumnName) : null;
    this.updateSizeFlags(sizeCol, sizeCol != null && sizeCol.isNumerical &&
      look.markerSizeScaling == DG.AxisType.logarithmic);
    // the dot renderer ignores sizes and types, the hit test does not
    const [rowSizes, markerSize] = this.updateSizes(sc, look, sizeCol);
    const [rowTypes, markerType] = this.updateTypes(sc, look, df);
    const selectionColors = look.showSelectedRows && look.rowSource != 'Selected' &&
      look.rowSource != 'FilteredSelected';
    const selection = selectionColors && this.updateSelection(df.selection);
    const rowColors = this.updateColors(sc, look, df, selectionColors);
    this.updateRows(sc, look, df, dots, selectionColors);

    const left = Math.floor(viewBox.x * dpr);
    const top = Math.floor(viewBox.y * dpr);
    const viewer = plotCanvases(sc).canvas;
    return {
      dots, dpr, viewBox, xform, origin: [left, top],
      size: [Math.max(1, Math.ceil((viewBox.x + viewBox.width) * dpr) - left),
        Math.max(1, Math.ceil((viewBox.y + viewBox.height) * dpr) - top)],
      rowColors, color: look.filteredRowsColor, rowSizes, markerSize, rowTypes, markerType,
      selection, selectedColor: look.selectedRowsColor,
      filteredOutColor: look.rowSource == 'All' ? null : look.filteredOutRowsColor,
      opacity: Math.trunc(Math.min(Math.max(look.markerOpacity, 0), 100) * 255 / 100) / 255,
      border: look.markerDrawBorder ? look.markerBorderWidth : 0,
      markerCache: viewer.width > 250 && viewer.height > 250,
      largestMarker: largestMarker(sc, look),
      rows: {filtered: this.filtered.length, filteredOut: this.filteredOutCount, selected: this.selectedCount},
    };
  }

  /** Column color coding lives in tags, which change no column version. */
  private watchMetadata(df: DG.DataFrame): void {
    if (this.dataFrame?.dart === df.dart)
      return;
    this.metadataSub?.unsubscribe();
    this.dataFrame = df;
    this.metadataSub = df.onMetadataChanged.subscribe(() => this.metadataVersion++);
    this.metadataVersion++;
  }

  private updateSizeFlags(sizeCol: DG.Column | null, log: boolean): void {
    const key = sizeCol ? `${objectId(sizeCol)}|${sizeCol.version}|${sizeCol.length}|${log}` : '';
    if (key == this.sizeFlagsKey)
      return;
    this.sizeFlagsKey = key;
    this.sizeFlagsVersion++;
    if (!sizeCol) {
      this.sizeFlags = new Uint8Array(0);
      this.emptySizes = false;
      return;
    }
    const n = sizeCol.length;
    const flags = new Uint8Array(n);
    const raw = sizeCol.getRawData();
    const fast = sizeCol.type == DG.COLUMN_TYPE.INT || sizeCol.type == DG.COLUMN_TYPE.FLOAT;
    const none = noneValue(sizeCol, raw as ArrayLike<number>);
    let emptySizes = false;
    for (let i = 0; i < n; i++) {
      const empty = fast ? raw[i] === none : sizeCol.isNone(i);
      emptySizes ||= empty;
      flags[i] = empty ? 1 : log && (fast ? raw[i] as number : sizeCol.getNumber(i)) <= 0 ? 2 : 0;
    }
    this.sizeFlags = flags;
    this.emptySizes = emptySizes;
  }

  private updateSizes(sc: DG.ScatterPlotViewer, look: ScatterLook, sizeCol: DG.Column | null): [boolean, number] {
    if (!sizeCol) {
      this.sizesKey = '';
      return [false, sc.getMarkerSize(0)];
    }
    const key = `${this.sizeFlagsKey}|${look.markerMinSize}|${look.markerMaxSize}`;
    if (key != this.sizesKey) {
      this.sizesKey = key;
      this.sink.setRowValues('sizes', sc.getMarkerSizes());
    }
    return [true, 0];
  }

  private updateTypes(sc: DG.ScatterPlotViewer, look: ScatterLook, df: DG.DataFrame): [boolean, number] {
    const markersCol = look.markersColumnName ? df.col(look.markersColumnName) : null;
    const emptySizes = look.showMarkersWithEmptySize && this.emptySizes;
    const type = Math.max(MARKER_TYPES.indexOf(look.markerType), 0);
    if (!markersCol && !emptySizes) {
      this.typesKey = '';
      return [false, type];
    }
    const key = [objectId(markersCol), markersCol?.version, look.markersMap, this.metadataVersion, type, emptySizes,
      this.sizeFlagsVersion].join('|');
    if (key != this.typesKey) {
      this.typesKey = key;
      const types = markersCol ? sc.getMarkerTypes() : new Uint32Array(df.rowCount).fill(type);
      if (emptySizes) {
        for (let i = 0; i < types.length; i++) {
          if (this.sizeFlags[i] == 1)
            types[i] = CIRCLE_BORDER;
        }
      }
      this.sink.setRowValues('types', types);
    }
    return [true, 0];
  }

  private updateSelection(selection: DG.BitSet): boolean {
    if (selection.trueCount == 0)
      return false;
    const key = `${objectId(selection)}|${selection.version}|${selection.length}`;
    if (key != this.selectionKey) {
      this.selectionKey = key;
      this.sink.setSelection(selection.getBuffer().subarray(0, (selection.length + 31) >>> 5));
    }
    return true;
  }

  /** The core paints selected rows in the selection color; the renderers do it from the selection instead, so a new
   * selection takes the core's colors again only when it deselects a row the core painted selected. */
  private updateColors(sc: DG.ScatterPlotViewer, look: ScatterLook, df: DG.DataFrame,
    selectionColors: boolean): boolean {
    if (!look.colorColumnName) {
      this.colorsKey = '';
      return false;
    }
    const selection = df.selection;
    const colorCol = df.col(look.colorColumnName);
    const key = [objectId(colorCol), colorCol?.version, look.colorColumnName, look.colorMap, look.showSelectedRows,
      look.rowSource, look.selectedRowsColor, look.filteredRowsColor, look.missingValueColor, look.colorAxisType,
      look.colorMin, look.colorMax, look.invertColorScheme, look.linearColorScheme, look.categoricalColorScheme,
      this.metadataVersion, objectId(selection), df.rowCount].join('|');
    const bits = selection.getBuffer();
    const taken = this.colorsSelection;
    let deselected = false;
    for (let w = 0; taken && w < taken.length && !deselected; w++)
      deselected = (taken[w] & ~bits[w]) != 0;
    if (key != this.colorsKey || deselected) {
      this.colorsKey = key;
      this.sink.setRowValues('colors', sc.getMarkerColors());
      this.colorsSelection = selectionColors ? bits.slice(0, (selection.length + 31) >>> 5) : null;
    }
    return true;
  }

  /** The rows of each pass: the filtered-out ones (when shown), the filtered ones, and the selected ones once more,
   * on top, as the canvas renderer draws them (the dot renderer does not). */
  private updateRows(sc: DG.ScatterPlotViewer, look: ScatterLook, df: DG.DataFrame, dots: boolean,
    selectionColors: boolean): void {
    const filter = sc.filter;
    const skipEmptySize = !dots && !look.showMarkersWithEmptySize;
    const key = [objectId(filter), filter.version, filter.length, this.x.version, this.y.version,
      this.sizeFlagsVersion, look.showFilteredOutPoints, skipEmptySize].join('|');
    if (key != this.filteredKey) {
      this.filteredKey = key;
      const n = filter.length;
      const bits = filter.getBuffer();
      const xInvalid = this.x.invalid;
      const yInvalid = this.y.invalid;
      const sizeFlags = this.sizeFlags;
      const hasSize = sizeFlags.length == n;
      const filtered = new Uint32Array(n);
      const filteredOut = new Uint32Array(look.showFilteredOutPoints ? n : 0);
      let f = 0;
      let o = 0;
      for (let i = 0; i < n; i++) {
        if (xInvalid[i] || yInvalid[i] || (hasSize && (sizeFlags[i] == 2 || (skipEmptySize && sizeFlags[i] == 1))))
          continue;
        if ((bits[i >>> 5] >>> (i & 31)) & 1)
          filtered[f++] = i;
        else if (filteredOut.length)
          filteredOut[o++] = i;
      }
      this.filtered = filtered.subarray(0, f);
      this.filteredVersion++;
      this.filteredOutCount = o;
      this.sink.setRows('filtered', this.filtered);
      this.sink.setRows('filteredOut', filteredOut.subarray(0, o));
    }

    const selection = df.selection;
    const onTop = !dots && selectionColors && selection.trueCount > 0;
    const selectedKey = onTop ? `${this.filteredVersion}|${objectId(selection)}|${selection.version}` : '';
    if (selectedKey != this.selectedKey) {
      this.selectedKey = selectedKey;
      const bits = selection.getBuffer();
      const rows = new Uint32Array(onTop ? this.filtered.length : 0);
      let count = 0;
      for (const i of rows.length ? this.filtered : NO_ROWS) {
        if ((bits[i >>> 5] >>> (i & 31)) & 1)
          rows[count++] = i;
      }
      this.selectedCount = count;
      this.sink.setRows('selected', rows.subarray(0, count));
    }
  }
}
