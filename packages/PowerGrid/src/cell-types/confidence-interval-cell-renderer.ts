import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

const CELL_TYPE = 'ConfidenceInterval';
const PADDING = 4;
const subscribedGrids = new WeakSet<DG.Grid>();
const DEFAULT_COLOR = '#4A90D9';
const WARNING_COLOR = '#E87C3E';


enum CIType {
  ThreeColumn = 'Value + Min Max',
  TwoColumn = 'Value + Margin',
}

enum ScaleType {
  Global = 'Global',
  PerRow = 'Per row',
  Custom = 'Custom',
  Symmetric = 'Symmetric around zero',
}

enum CenterMark {
  Dot = 'Dot',
  Diamond = 'Diamond',
  Line = 'Line',
}

interface CISettings {
  type: CIType;
  estimateColumn: string;
  lowerColumn: string;      // lower bound (3-col) or margin of error (2-col)
  upperColumn: string;      // upper bound (3-col only)
  scaleType: ScaleType;
  customMin: number;
  customMax: number;
  referenceLine: number | null;
  logScale: boolean;
  centerMark: CenterMark;
  showSerifs: boolean;
  color: string;
  colorColumnName: string;
}

interface RowValues {
  estimate: number | null;
  lower: number | null;
  upper: number | null;
}

/** Samples up to `count` row indices spread uniformly across [0, total). */
function sampleRows(total: number, count: number): Int32Array {
  if (count >= total) {
    const rows = new Int32Array(total);
    for (let i = 0; i < total; i++) rows[i] = i;
    return rows;
  }
  const step = total / count;
  const rows = new Int32Array(count);
  for (let i = 0; i < count; i++)
    rows[i] = Math.floor(i * step + Math.random() * step);
  return rows;
}

/**
 * Detects three numerical columns where c1 < c2 < c3 holds for all sampled rows.
 * Sorts columns by min value, then verifies consecutive pairs — O(n × sampleSize).
 */
function detectCIColumns(df: DG.DataFrame): {lower: string, estimate: string, upper: string} | null {
  const sorted = [...df.columns.numerical]
    .filter((c) => c.type !== DG.TYPE.DATE_TIME && isFinite(c.min))
    .sort((a, b) => a.min - b.min);
  const n = sorted.length;
  if (n < 3 || df.rowCount === 0) return null;

  const rows = sampleRows(df.rowCount, Math.min(100, df.rowCount));

  // For each consecutive pair, check if strict ordering holds across all sampled rows
  const holds = new Uint8Array(n - 1).fill(1);
  for (let si = 0; si < rows.length; si++) {
    const row = rows[si];
    for (let p = 0; p < n - 1; p++) {
      if (!holds[p]) continue;
      if (sorted[p].isNone(row) || sorted[p + 1].isNone(row)) continue;
      if (sorted[p].getNumber(row) >= sorted[p + 1].getNumber(row))
        holds[p] = 0;
    }
  }

  // Find two consecutive holding pairs → a valid triple
  for (let p = 0; p < n - 2; p++) {
    if (holds[p] && holds[p + 1])
      return {lower: sorted[p].name, estimate: sorted[p + 1].name, upper: sorted[p + 2].name};
  }
  return null;
}

function getSettings(gc: DG.GridColumn): CISettings {
  gc.settings ??= {};
  const isNew = gc.settings[CELL_TYPE] === undefined;
  const s: CISettings = gc.settings[CELL_TYPE] ??= {};
  const df = gc.grid.dataFrame;

  if (isNew) {
    const detected = detectCIColumns(df);
    if (detected) {
      s.estimateColumn = detected.estimate;
      s.lowerColumn = detected.lower;
      s.upperColumn = detected.upper;
    }
  }

  const numCols = [...df.columns.numerical].filter((c) => c.type !== DG.TYPE.DATE_TIME);
  s.type ??= CIType.ThreeColumn;
  s.estimateColumn ??= numCols[0]?.name ?? '';
  s.lowerColumn ??= numCols[1]?.name ?? '';
  s.upperColumn ??= numCols[2]?.name ?? '';
  s.scaleType ??= ScaleType.Global;
  s.customMin ??= 0;
  s.customMax ??= 100;
  if (s.referenceLine === undefined) s.referenceLine = null;
  s.logScale ??= false;
  s.centerMark ??= CenterMark.Dot;
  s.showSerifs ??= true;
  s.color ??= DEFAULT_COLOR;
  s.colorColumnName ??= '';
  return s;
}

function getCol(df: DG.DataFrame, name: string): DG.Column | null {
  if (!name) return null;
  try { return df.columns.byName(name); }
  catch { return null; }
}

function getRowValues(df: DG.DataFrame, settings: CISettings, row: number): RowValues {
  const estCol = getCol(df, settings.estimateColumn);
  const estimate = estCol && !estCol.isNone(row) ? estCol.getNumber(row) : null;

  if (settings.type === CIType.ThreeColumn) {
    const loCol = getCol(df, settings.lowerColumn);
    const hiCol = getCol(df, settings.upperColumn);
    return {
      estimate,
      lower: loCol && !loCol.isNone(row) ? loCol.getNumber(row) : null,
      upper: hiCol && !hiCol.isNone(row) ? hiCol.getNumber(row) : null,
    };
  }
  else {
    const marginCol = getCol(df, settings.lowerColumn);
    const margin = marginCol && !marginCol.isNone(row) ? marginCol.getNumber(row) : null;
    return {
      estimate,
      lower: estimate !== null && margin !== null ? estimate - margin : null,
      upper: estimate !== null && margin !== null ? estimate + margin : null,
    };
  }
}

function computeGlobalRange(gc: DG.GridColumn, settings: CISettings): {min: number, max: number} {
  const cached = gc.temp['ci_range'];
  if (cached) return cached;

  const df = gc.grid.dataFrame;
  let min = Infinity, max = -Infinity;

  if (settings.type === CIType.ThreeColumn) {
    for (const name of [settings.estimateColumn, settings.lowerColumn, settings.upperColumn]) {
      const col = getCol(df, name);
      if (!col) continue;
      if (col.min < min) min = col.min;
      if (col.max > max) max = col.max;
    }
  }
  else {
    const estCol = getCol(df, settings.estimateColumn);
    const marginCol = getCol(df, settings.lowerColumn);
    if (estCol && marginCol) {
      for (let i = 0; i < df.rowCount; i++) {
        if (estCol.isNone(i) || marginCol.isNone(i)) continue;
        const est = estCol.getNumber(i);
        const margin = marginCol.getNumber(i);
        if (est - margin < min) min = est - margin;
        if (est + margin > max) max = est + margin;
      }
    }
  }

  if (!isFinite(min) || !isFinite(max)) { min = 0; max = 1; }
  const result = {min, max};
  gc.temp['ci_range'] = result;
  return result;
}

function getScale(gc: DG.GridColumn, settings: CISettings, row: number): {min: number, max: number} {
  let min: number, max: number;

  if (settings.scaleType === ScaleType.Custom) {
    min = settings.customMin;
    max = settings.customMax;
  }
  else if (settings.scaleType === ScaleType.PerRow) {
    const vals = getRowValues(gc.grid.dataFrame, settings, row);
    const nums = [vals.estimate, vals.lower, vals.upper].filter((v) => v !== null && isFinite(v)) as number[];
    if (nums.length === 0) return {min: 0, max: 1};
    min = Math.min(...nums);
    max = Math.max(...nums);
    const range = max - min || Math.abs(max) * 0.2 || 1;
    min -= range * 0.15;
    max += range * 0.15;
  }
  else {
    const range = computeGlobalRange(gc, settings);
    min = range.min;
    max = range.max;
    if (settings.scaleType === ScaleType.Symmetric) {
      const absMax = Math.max(Math.abs(min), Math.abs(max));
      min = -absMax;
      max = absMax;
    }
  }

  if (settings.logScale) {
    min = min > 0 ? Math.log10(min) : 0;
    max = max > 0 ? Math.log10(max) : 1;
  }
  return {min, max};
}

function mapX(value: number, scale: {min: number, max: number}, logScale: boolean,
  left: number, width: number): number {
  const v = logScale && value > 0 ? Math.log10(value) : value;
  if (scale.max === scale.min) return left + width / 2;
  return left + ((v - scale.min) / (scale.max - scale.min)) * width;
}

function getCellColor(settings: CISettings, gc: DG.GridColumn, row: number): string {
  if (settings.colorColumnName) {
    const col = getCol(gc.grid.dataFrame, settings.colorColumnName);
    if (col) {
      const c = col.meta.colors.getColor(row);
      if (c !== 0)
        return DG.Color.toHtml(c);
    }
  }
  return settings.color;
}

function drawArrowHead(g: CanvasRenderingContext2D, x: number, y: number,
  size: number, direction: number, color: string) {
  g.fillStyle = color;
  g.beginPath();
  g.moveTo(x, y - size);
  g.lineTo(x + direction * size * 0.8, y);
  g.lineTo(x, y + size);
  g.closePath();
  g.fill();
}


function getHeaderScale(gc: DG.GridColumn, settings: CISettings): {min: number, max: number} | null {
  if (settings.scaleType === ScaleType.Custom)
    return {min: settings.customMin, max: settings.customMax};

  const range = computeGlobalRange(gc, settings);
  let {min, max} = range;
  if (!isFinite(min) || !isFinite(max)) return null;

  if (settings.scaleType === ScaleType.Symmetric) {
    const absMax = Math.max(Math.abs(min), Math.abs(max));
    min = -absMax;
    max = absMax;
  }

  if (settings.logScale) {
    min = min > 0 ? Math.log10(min) : 0;
    max = max > 0 ? Math.log10(max) : 1;
  }
  return {min, max};
}

function renderHeader(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  gc: DG.GridColumn): void {
  const settings = getSettings(gc);
  const scale = getHeaderScale(gc, settings);

  g.save();
  g.rect(x, y, w, h);
  g.clip();

  g.fillStyle = '#333333';
  g.font = '12px Roboto, "Roboto Local"';
  g.textBaseline = 'top';
  g.fillText(gc.name, x + PADDING, y + 4);

  if (scale) {
    const axisH = 24;
    DG.Paint.horzAxis(g, scale.min, scale.max,
      x + PADDING, y + h - axisH, w - 2 * PADDING, axisH,
      settings.logScale, false);
  }

  g.restore();
}

export class ConfidenceIntervalCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Confidence Interval'; }
  get cellType() { return CELL_TYPE; }

  getDefaultSize(gridColumn: DG.GridColumn): {width?: number | null, height?: number | null} {
    return {width: 150, height: null};
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const grid = gridCell.grid;
    // if (!subscribedGrids.has(grid)) {
    //   subscribedGrids.add(grid);
    //   grid.onCellRender.subscribe((args) => {
    //     if (!args.cell.isColHeader || args.cell.gridColumn.cellType !== CELL_TYPE)
    //       return;
    //     renderHeader(args.g, args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height,
    //       args.cell.gridColumn);
    //     args.preventDefault();
    //   });
    // }

    const df = gridCell.grid.dataFrame;
    if (w < 20 || h < 10 || !df) return;

    const settings = getSettings(gridCell.gridColumn);
    const row = gridCell.cell.row.idx;
    const vals = getRowValues(df, settings, row);

    if (vals.estimate === null && vals.lower === null && vals.upper === null)
      return;

    const left = x + PADDING;
    const bw = w - 2 * PADDING;
    const cy = y + h / 2;
    const scale = getScale(gridCell.gridColumn, settings, row);
    const mx = (v: number) => mapX(v, scale, settings.logScale, left, bw);
    const baseColor = getCellColor(settings, gridCell.gridColumn, row);
    const inverted = vals.lower !== null && vals.upper !== null && vals.lower > vals.upper;
    const drawColor = inverted ? WARNING_COLOR : baseColor;

    // Reference line
    if (settings.referenceLine !== null) {
      const rx = mx(settings.referenceLine);
      if (rx >= left && rx <= left + bw) {
        g.strokeStyle = '#B0B0B0';
        g.lineWidth = 1;
        g.setLineDash([2, 2]);
        g.beginPath();
        g.moveTo(rx, y + 2);
        g.lineTo(rx, y + h - 2);
        g.stroke();
        g.setLineDash([]);
      }
    }

    // CI range
    if (vals.lower !== null && vals.upper !== null) {
      const lo = Math.min(vals.lower, vals.upper);
      const hi = Math.max(vals.lower, vals.upper);
      let loX = mx(lo);
      let hiX = mx(hi);
      const loClipped = !isFinite(lo) || loX < left;
      const hiClipped = !isFinite(hi) || hiX > left + bw;
      loX = Math.max(loX, left);
      hiX = Math.min(hiX, left + bw);

      // Whisker line
      g.strokeStyle = drawColor;
      g.lineWidth = 1.5;
      if (inverted) g.setLineDash([3, 2]);
      g.beginPath();
      g.moveTo(loX, cy);
      g.lineTo(hiX, cy);
      g.stroke();
      if (inverted) g.setLineDash([]);

      // Serifs (end caps)
      if (settings.showSerifs) {
        const serifH = Math.min(h * 0.3, 8);
        g.beginPath();
        if (!loClipped) {
          g.moveTo(loX, cy - serifH);
          g.lineTo(loX, cy + serifH);
        }
        if (!hiClipped) {
          g.moveTo(hiX, cy - serifH);
          g.lineTo(hiX, cy + serifH);
        }
        g.stroke();
      }
      g.lineWidth = 1;

      // Clipping/infinite arrows
      const arrowSize = Math.min(h * 0.15, 5);
      if (loClipped)
        drawArrowHead(g, left + arrowSize * 0.8, cy, arrowSize, -1, drawColor);
      if (hiClipped)
        drawArrowHead(g, left + bw - arrowSize * 0.8, cy, arrowSize, 1, drawColor);
    }

    // Center mark (point estimate)
    if (vals.estimate !== null) {
      const estX = mx(vals.estimate);
      if (estX >= left && estX <= left + bw) {
        g.fillStyle = drawColor;
        g.strokeStyle = drawColor;
        const markSize = Math.min(h * 0.2, 5);

        switch (settings.centerMark) {
        case CenterMark.Dot:
          g.beginPath();
          g.arc(estX, cy, markSize, 0, Math.PI * 2);
          g.fill();
          break;
        case CenterMark.Diamond:
          g.beginPath();
          g.moveTo(estX, cy - markSize * 1.3);
          g.lineTo(estX + markSize, cy);
          g.lineTo(estX, cy + markSize * 1.3);
          g.lineTo(estX - markSize, cy);
          g.closePath();
          g.fill();
          break;
        case CenterMark.Line:
          g.lineWidth = 2;
          g.beginPath();
          g.moveTo(estX, cy - Math.min(h * 0.3, 8));
          g.lineTo(estX, cy + Math.min(h * 0.3, 8));
          g.stroke();
          g.lineWidth = 1;
          break;
        }
      }
    }
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const settings = getSettings(gridCell.gridColumn);
    const vals = getRowValues(gridCell.grid.dataFrame, settings, gridCell.cell.row.idx);

    if (vals.estimate === null && vals.lower === null && vals.upper === null) {
      ui.tooltip.hide();
      return;
    }

    const keysDiv = ui.divV([], {style: {marginRight: '10px', fontWeight: 'bold', textAlign: 'right'}});
    const valuesDiv = ui.divV([], {style: {fontWeight: 'normal'}});
    const add = (label: string, value: number | null) => {
      if (value === null) return;
      keysDiv.appendChild(ui.divText(label));
      valuesDiv.appendChild(ui.divText(String(parseFloat(value.toPrecision(4)))));
    };

    add('Estimate', vals.estimate);
    add('Lower', vals.lower);
    add('Upper', vals.upper);
    if (vals.lower !== null && vals.upper !== null)
      add('CI Width', Math.abs(vals.upper - vals.lower));
    if (vals.lower !== null && vals.upper !== null && vals.lower > vals.upper)
      keysDiv.appendChild(ui.divText('⚠ Inverted bounds'));

    ui.tooltip.show(ui.divH([keysDiv, valuesDiv], {style: {display: 'flex'}}), e.x + 16, e.y + 16);
  }

  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    ui.tooltip.hide();
  }

  renderSettings(gc: DG.GridColumn): Element {
    const settings = getSettings(gc);
    const df = gc.grid.dataFrame;
    const numColNames = [...df.columns.numerical]
      .filter((c) => c.type !== DG.TYPE.DATE_TIME)
      .map((c) => c.name);
    const invalidate = () => {
      delete gc.temp['ci_range'];
      gc.grid.invalidate();
    };

    const typeInput = ui.input.choice('Type', {
      value: settings.type,
      items: [CIType.ThreeColumn, CIType.TwoColumn],
      onValueChanged: (v) => { settings.type = v!; updateVisibility(); invalidate(); },
      nullable: false,
    });

    const estimateInput = ui.input.choice('Estimate', {
      value: settings.estimateColumn,
      items: numColNames,
      onValueChanged: (v) => { settings.estimateColumn = v ?? ''; invalidate(); },
    });

    const lowerInput = ui.input.choice('Lower', {
      value: settings.lowerColumn,
      items: numColNames,
      onValueChanged: (v) => { settings.lowerColumn = v ?? ''; invalidate(); },
    });

    const upperInput = ui.input.choice('Upper', {
      value: settings.upperColumn,
      items: numColNames,
      onValueChanged: (v) => { settings.upperColumn = v ?? ''; invalidate(); },
    });

    const scaleInput = ui.input.choice('Scale', {
      value: settings.scaleType,
      items: [ScaleType.Global, ScaleType.PerRow, ScaleType.Custom, ScaleType.Symmetric],
      onValueChanged: (v) => { settings.scaleType = v!; updateVisibility(); invalidate(); },
      nullable: false,
    });

    const customMinInput = ui.input.float('Min', {
      value: settings.customMin,
      onValueChanged: (v) => { if (v !== null) { settings.customMin = v; invalidate(); } },
    });

    const customMaxInput = ui.input.float('Max', {
      value: settings.customMax,
      onValueChanged: (v) => { if (v !== null) { settings.customMax = v; invalidate(); } },
    });

    const refLineInput = ui.input.float('Reference', {
      value: settings.referenceLine ?? undefined as any,
      nullable: true,
      onValueChanged: (v) => { settings.referenceLine = v; invalidate(); },
    });

    const logScaleInput = ui.input.bool('Log Scale', {
      value: settings.logScale,
      onValueChanged: (v) => { settings.logScale = v; invalidate(); },
    });

    const centerMarkInput = ui.input.choice('Center Mark', {
      value: settings.centerMark,
      items: [CenterMark.Dot, CenterMark.Diamond, CenterMark.Line],
      onValueChanged: (v) => { settings.centerMark = v!; invalidate(); },
      nullable: false,
    });

    const showSerifsInput = ui.input.bool('Serifs', {
      value: settings.showSerifs,
      onValueChanged: (v) => { settings.showSerifs = v; invalidate(); },
    });

    const colorInput = ui.input.color('Color', {
      value: settings.color,
      onValueChanged: (v) => { settings.color = v; invalidate(); },
    });

    const allColNames = ['', ...df.columns.names()];
    const colorColInput = ui.input.choice('Color Column', {
      value: settings.colorColumnName,
      items: allColNames,
      onValueChanged: (v) => { settings.colorColumnName = v ?? ''; invalidate(); },
    });

    function updateVisibility() {
      const isTwoCol = settings.type === CIType.TwoColumn;
      upperInput.root.style.display = isTwoCol ? 'none' : '';
      lowerInput.caption = isTwoCol ? 'Margin' : 'Lower';

      const isCustom = settings.scaleType === ScaleType.Custom;
      customMinInput.root.style.display = isCustom ? '' : 'none';
      customMaxInput.root.style.display = isCustom ? '' : 'none';
    }
    updateVisibility();

    return ui.inputs([
      typeInput, estimateInput, lowerInput, upperInput,
      scaleInput, customMinInput, customMaxInput, refLineInput, logScaleInput,
      centerMarkInput, showSerifsInput,
      colorInput, colorColInput,
    ]);
  }
}
