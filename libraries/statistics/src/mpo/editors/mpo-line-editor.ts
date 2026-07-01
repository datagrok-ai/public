/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type Konva from 'konva';
import {DesirabilityLine, NumericalDesirability, toScale, fromScale, domainMinX, domainMaxX, MpoScale} from '../mpo';
import {Subject} from 'rxjs';

let _konva: typeof Konva | undefined;
let _konvaPromise: Promise<typeof Konva> | undefined;

async function getKonva(): Promise<typeof Konva> {
  if (_konva)
    return _konva;
  _konvaPromise ??= import('konva').then((m) => m.default);
  _konva = await _konvaPromise;
  return _konva;
}

type Point = [number, number];

// Constants for the editor layout
const EDITOR_PADDING = {top: 10, right: 10, bottom: 20, left: 30};
const POINT_RADIUS = 3;

const COLORS = {
  line: DG.Color.toHtml(DG.Color.filteredRows),
  handle: DG.Color.toHtml(DG.Color.selectedRows),
  barFill: DG.Color.toHtml(DG.Color.histogramBar),
  barStroke: DG.Color.toHtml(DG.Color.lightGray),
  warning: DG.Color.toHtml(DG.Color.red),
};

class CoordMapper {
  private plotWidth: number;
  private plotHeight: number;
  private scaleX: number;
  private scaleY: number;
  private tMin: number;
  private tMax: number;

  constructor(
    private minX: number,
    private maxX: number,
    private width: number,
    private height: number,
    private inverted = false,
    private log = false,
  ) {
    this.plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
    this.plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;
    this.tMin = toScale(minX, log);
    this.tMax = toScale(maxX, log);
    this.scaleX = (this.tMax - this.tMin === 0) ? 1 : this.plotWidth / (this.tMax - this.tMin);
    this.scaleY = this.plotHeight;
  }

  toCanvasCoords(p: Point): {x: number, y: number} {
    const yDisp = this.inverted ? 1 - p[1] : p[1];
    const px = this.log ? Math.max(this.minX, p[0]) : p[0];
    const canvasX = EDITOR_PADDING.left + (toScale(px, this.log) - this.tMin) * this.scaleX;
    const canvasY = EDITOR_PADDING.top + this.plotHeight - (yDisp * this.scaleY);

    return {x: canvasX, y: canvasY};
  }

  toDataCoords(canvasX: number, canvasY: number): {x: number, y: number} {
    const tx = this.tMin + (canvasX - EDITOR_PADDING.left) / this.scaleX;
    let dataX = fromScale(tx, this.log);
    let dataY = (EDITOR_PADDING.top + this.plotHeight - canvasY) / this.scaleY;

    dataX = Math.max(this.minX, Math.min(this.maxX, dataX));
    dataY = Math.max(0, Math.min(1, dataY));
    if (this.inverted)
      dataY = 1 - dataY;

    return {x: dataX, y: dataY};
  }
}

export class MpoDesirabilityLineEditor {
  root = ui.div([], 'statistics-mpo-line-editor');
  onChanged = new Subject<DesirabilityLine>();
  supportsModeDialog: boolean = true;

  private _prop: NumericalDesirability;
  private barsLayer?: Konva.Layer;
  private pendingBarValues?: number[];

  private stage?: Konva.Stage;
  private layer?: Konva.Layer;

  private konvaLine?: Konva.Line;
  private pointsGroup?: Konva.Group;
  private barValues?: number[];

  private redrawFn!: (notify?: boolean) => void;

  private specialHandle?: Konva.Circle;
  onParamsChanged?: (prop: NumericalDesirability) => void;

  // Flag to prevent touchpad right-click from adding a new point
  private ignoreNextClick = false;
  private dragScaleY = 0;
  private _width: number;
  private _height: number;

  constructor(prop: NumericalDesirability, width: number, height: number) {
    this._prop = prop;
    this._width = width;
    this._height = height;
    this.ensureDefaultLine();
    this.root.style.width = `${width}px`;
    this.root.style.height = `${height}px`;
    this.root.style.position = 'relative';

    requestAnimationFrame(() => this.initKonva(width, height));
  }

  private ensureDefaultLine(): void {
    if (this._prop.line.length > 0)
      return;

    const min = this._prop.min ?? 0;
    const max = this._prop.max ?? 1;
    this._prop.line = [[min, 0.5], [max, 0.5]];
  }

  private updateDragScales(): void {
    const plotHeight = this._height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;
    this.dragScaleY = 1 / plotHeight;
  }

  private isInPlotArea(pos: {x: number, y: number}, width: number, height: number) {
    return (
      pos.x >= EDITOR_PADDING.left &&
      pos.x <= width - EDITOR_PADDING.right &&
      pos.y >= EDITOR_PADDING.top &&
      pos.y <= height - EDITOR_PADDING.bottom
    );
  }

  private async initKonva(width: number, height: number) {
    if (!this.root.parentElement) {
      console.warn('Konva container not attached to DOM yet.');
      return;
    }

    const Konva = await getKonva();
    this.barsLayer = new Konva.Layer();

    this.stage = new Konva.Stage({
      container: this.root,
      width,
      height,
    });

    this.stage.add(this.barsLayer);
    this.layer = new Konva.Layer();
    this.stage.add(this.layer);

    const minX = this.getMinX();
    const maxX = this.getMaxX();

    // --- Draw Axes ---
    this.drawAxes(minX, maxX, width, height);

    // --- Draw Line and Points ---
    this.konvaLine = new Konva.Line({
      points: [],
      stroke: COLORS.line,
      strokeWidth: 2,
      lineCap: 'round',
      lineJoin: 'round',
    });
    this.layer.add(this.konvaLine);

    this.pointsGroup = new Konva.Group();
    this.layer.add(this.pointsGroup);

    // --- Redraw Function ---
    this.redrawFn = (notify: boolean = true) => {
      const minX = this.getMinX();
      const maxX = this.getMaxX();

      if (this._prop.mode === 'freeform' && this._prop.freeformLine)
        this._prop.line = this._prop.freeformLine;

      if (this._prop.mode !== 'freeform') {
        if (!this._prop.freeformLine)
          this._prop.freeformLine = [...this._prop.line];
        this._prop.line = this.computeLine();
      }

      this.pointsGroup!.destroyChildren();
      const konvaPoints: number[] = [];

      const mapper = new CoordMapper(minX, maxX, this._width, this._height, !!this._prop.inverted, this.isLog);
      const sortedIndices = [...this._prop.line.keys()]
        .sort((a, b) => this._prop.line[a][0] - this._prop.line[b][0]);
      const sortedLine = sortedIndices.map((idx) => this._prop.line[idx]);

      sortedLine.forEach((p, index) => {
        const coords = mapper.toCanvasCoords([p[0], p[1]]);
        konvaPoints.push(coords.x, coords.y);

        if (this._prop.mode !== 'freeform')
          return;

        const pointCircle = new Konva.Circle({
          x: coords.x,
          y: coords.y,
          radius: POINT_RADIUS,
          fill: 'white',
          stroke: COLORS.line,
          strokeWidth: 1,
          draggable: true,
          hitStrokeWidth: 5,
        });

        // Store index directly on the node for easy access
        pointCircle.setAttr('_dataIndex', sortedIndices[index]);

        pointCircle.on('dragmove', (evt) => {
          const circle = evt.target as Konva.Circle;
          const pos = circle.position();
          const dataIndex = circle.getAttr('_dataIndex');
          const sortedPos = sortedIndices.indexOf(dataIndex);

          // Constrain dragging horizontally between neighbors (or bounds)
          const prevX = sortedPos > 0 ? this._prop.line[sortedIndices[sortedPos - 1]][0] : minX;
          const nextX = sortedPos < sortedIndices.length - 1 ?
            this._prop.line[sortedIndices[sortedPos + 1]][0] :
            maxX;

          // Add a small buffer to avoid points overlapping exactly
          const buffer = (maxX - minX === 0) ? 0 : 0.001 * (maxX - minX); // Avoid NaN if minX === maxX
          const minCanvasX = mapper.toCanvasCoords([prevX + (sortedPos > 0 ? buffer : 0), 0]).x;
          const maxCanvasX = mapper.toCanvasCoords([nextX - (sortedPos < sortedIndices.length - 1 ? buffer : 0), 0]).x;

          pos.x = Math.max(minCanvasX, Math.min(maxCanvasX, pos.x));

          const plotTop = EDITOR_PADDING.top;
          const plotBottom = this._height - EDITOR_PADDING.bottom;
          pos.y = Math.max(plotTop, Math.min(plotBottom, pos.y));

          circle.position(pos);

          const dataCoords = mapper.toDataCoords(pos.x, pos.y);
          this._prop.line[dataIndex][0] = dataCoords.x;
          this._prop.line[dataIndex][1] = dataCoords.y;

          const currentKonvaPoints = this._prop.line.map((pData, idx) => {
            // Use dragged circle position directly for the point being dragged
            if (idx === dataIndex)
              return [pos.x, pos.y];
            else {
              const c = mapper.toCanvasCoords([pData[0], pData[1]]);
              return [c.x, c.y];
            }
          }).flat();

          this.konvaLine!.points(currentKonvaPoints);
          this.layer!.batchDraw();
        });

        pointCircle.on('dragend', () => {
          this._prop.line.sort((a, b) => a[0] - b[0]);
          this.redrawFn?.();
        });

        pointCircle.on('contextmenu', (evt) => {
          evt.evt.preventDefault();
          if (this._prop.line.length <= 2) {
            grok.shell.warning('Cannot remove points, minimum of 2 required.');
            return;
          }

          this.ignoreNextClick = true;
          const circle = evt.target as Konva.Circle;
          const dataIndex = circle.getAttr('_dataIndex') as number;

          if (dataIndex >= 0) {
            if (this._prop.min == null)
              this._prop.min = this.getMinX();
            if (this._prop.max == null)
              this._prop.max = this.getMaxX();
            this._prop.line.splice(dataIndex, 1);
            this.redrawFn?.();
          }
        });

        pointCircle.on('mouseenter', (evt) => {
          this.stage!.container().style.cursor = 'pointer';
          const circle = evt.target as Konva.Circle;
          const pos = circle.position();
          const dataCoords = mapper.toDataCoords(pos.x, pos.y);
          const tooltipText = `X: ${dataCoords.x.toFixed(2)}, Y: ${dataCoords.y.toFixed(2)}<br><br>Drag to move, double-click to edit, right-click to delete`;
          ui.tooltip.show(tooltipText, evt.evt.clientX, evt.evt.clientY);
        });

        pointCircle.on('mouseleave', () => {
          this.stage!.container().style.cursor = 'default';
          ui.tooltip.hide();
        });

        pointCircle.on('dblclick dbltap', (evt) => {
          this.ignoreNextClick = true;
          ui.tooltip.hide();
          const dataIndex = (evt.target as Konva.Circle).getAttr('_dataIndex') as number;
          this.showPointEditor(dataIndex, evt.evt.clientX, evt.evt.clientY);
        });

        this.pointsGroup!.add(pointCircle);
      });

      this.konvaLine!.points(konvaPoints);
      this.layer!.batchDraw();

      // --- Add special handle (Gaussian peak / Sigmoid inflection) ---
      this.addSpecialHandle(this._width, this._height);

      if (notify)
        this.onChanged.next(this._prop.line);
    };

    // --- Left-click to Add Point ---
    this.stage.on('click tap', (evt) => {
      if (this.ignoreNextClick) {
        this.ignoreNextClick = false;
        return;
      }

      if (evt.target instanceof Konva.Circle || evt.evt.button !== 0)
        return;

      const pos = this.stage?.getPointerPosition();
      if (!pos)
        return;

      const minX = this.getMinX();
      const maxX = this.getMaxX();

      if (!this.isInPlotArea(pos, this.stage!.width(), this.stage!.height()))
        return;

      const mapper = new CoordMapper(minX, maxX, this.stage!.width(), this.stage!.height(), !!this._prop.inverted, this.isLog);
      const dataCoords = mapper.toDataCoords(pos.x, pos.y);

      this._prop.line.push([dataCoords.x, dataCoords.y]);
      this._prop.line.sort((a, b) => a[0] - b[0]);
      this.redrawFn?.();
    });

    // Cursor change
    this.stage.on('mousemove', (evt) => {
      if (!this.stage)
        return;
      const pos = this.stage.getPointerPosition();
      if (!pos)
        return;

      if (this._prop.mode !== 'freeform') {
        if (this.isInPlotArea(pos, this._width, this._height))
          this.stage.container().style.cursor = 'grab';
        else
          this.stage.container().style.cursor = 'default';
      }
    });

    this.stage.on('mouseout', () => ui.tooltip.hide());

    // Enable curve drag (smooth)
    this.enableCurveDrag();
    this.updateDragScales();

    ui.onSizeChanged(this.root).subscribe((_) => {
      const {clientWidth: w, clientHeight: h} = this.root;
      if (w > 0 && h > 0)
        this.resize(w, h);
    });

    // Initial draw
    if (this.pendingBarValues) {
      this.drawBars(this.pendingBarValues);
      this.pendingBarValues = undefined;
    }
    this.redrawFn?.(false);
  }

  get line(): DesirabilityLine {
    return this._prop.line;
  }

  private computeLine(): DesirabilityLine {
    if (this._prop.mode === 'freeform')
      return this._prop.line;

    const log = this.isLog;
    const tMin = toScale(this.getMinX(), log);
    const tMax = toScale(this.getMaxX(), log);
    const n = 60;
    const line: DesirabilityLine = [];

    for (let i = 0; i <= n; ++i) {
      const u = tMin + (tMax - tMin) * (i / n);
      let y = 0;

      if (this._prop.mode === 'gaussian') {
        const meanVal = this._prop.mean ?? this.getDefaultMean();
        const sigma = this._prop.sigma ?? this.getDefaultSigma();
        const mean = toScale(log ? Math.max(this.getMinX(), meanVal) : meanVal, log);
        const z = (u - mean) / sigma;
        y = Math.exp(-0.5 * z * z);
      }

      if (this._prop.mode === 'sigmoid') {
        const x0Val = this._prop.x0 ?? this.getDefaultX0();
        const k = this._prop.k ?? this.getDefaultK();
        const x0 = toScale(log ? Math.max(this.getMinX(), x0Val) : x0Val, log);
        y = 1 / (1 + Math.exp(-k * (u - x0)));
      }

      line.push([fromScale(u, log), y]);
    }

    return line;
  }

  private drawAxes(minX: number, maxX: number, width: number, height: number) {
    const axisColor = getComputedStyle(document.documentElement).getPropertyValue('--grey-2').trim() || '#DBDCDF';
    this.layer!.add(
      new _konva!.Line({
        points: [EDITOR_PADDING.left, height - EDITOR_PADDING.bottom, width - EDITOR_PADDING.right, height - EDITOR_PADDING.bottom],
        stroke: axisColor,
        strokeWidth: 1,
      }),
      new _konva!.Line({
        points: [EDITOR_PADDING.left, EDITOR_PADDING.top, EDITOR_PADDING.left, height - EDITOR_PADDING.bottom],
        stroke: axisColor,
        strokeWidth: 1,
      }),
      new _konva!.Text({x: EDITOR_PADDING.left, y: height - EDITOR_PADDING.bottom + 3, text: this.formatTick(minX), fontSize: 9, fill: 'grey'}),
      new _konva!.Text({x: width - EDITOR_PADDING.right - 15, y: height - EDITOR_PADDING.bottom + 3, text: this.formatTick(maxX), fontSize: 9, fill: 'grey'}),
    );

    const hidden = this.isLog && this.barValues ? this.barValues.filter((v) => Number.isFinite(v) && v <= 0).length : 0;
    if (hidden > 0) {
      this.layer!.add(new _konva!.Text({x: EDITOR_PADDING.left, y: height - EDITOR_PADDING.bottom + 12,
        width: width - EDITOR_PADDING.left - EDITOR_PADDING.right, align: 'center',
        text: `${hidden} non-positive value${hidden > 1 ? 's' : ''} not shown`, fontSize: 8, fill: COLORS.warning}));
    }
  }

  private get isLog(): boolean {
    return this._prop.scale === MpoScale.Log;
  }

  private formatTick(v: number): string {
    if (v === 0)
      return '0';
    const a = Math.abs(v);
    if (a < 1e-4 || a >= 1e5)
      return v.toExponential(0);
    return parseFloat(v.toPrecision(3)).toString();
  }

  getMinX(): number {
    return domainMinX(this._prop);
  }

  getMaxX(): number {
    return domainMaxX(this._prop);
  }

  getDefaultMean(): number {
    return fromScale((toScale(this.getMinX(), this.isLog) + toScale(this.getMaxX(), this.isLog)) / 2, this.isLog);
  }

  getDefaultSigma(): number {
    return Math.max(0.01, (toScale(this.getMaxX(), this.isLog) - toScale(this.getMinX(), this.isLog)) / 6);
  }

  getDefaultX0(): number {
    return fromScale((toScale(this.getMinX(), this.isLog) + toScale(this.getMaxX(), this.isLog)) / 2, this.isLog);
  }

  getDefaultK(): number {
    return 10;
  }

  redrawAll(notify: boolean = true): void {
    if (!this.stage || !this.layer || !this.redrawFn)
      return;

    const width = this.stage.width();
    const height = this.stage.height();
    const minX = this.getMinX();
    const maxX = this.getMaxX();

    this.konvaLine?.remove();
    this.pointsGroup?.remove();
    this.layer.destroyChildren();
    this.specialHandle = undefined;

    this.drawAxes(minX, maxX, width, height);

    this.layer.add(this.konvaLine!, this.pointsGroup!);

    this.redrawFn(notify);

    if (this.barValues)
      this.drawBars();
  }

  resize(width: number, height: number): void {
    if (width === this._width && height === this._height)
      return;

    this._width = width;
    this._height = height;

    if (this.stage) {
      this.stage.width(width);
      this.stage.height(height);
      this.updateDragScales();
      this.redrawAll(false);
    }
  }

  drawBars(values?: number[]) {
    if (values)
      this.barValues = values;

    if (!this.barsLayer) {
      this.pendingBarValues = values;
      return;
    }
    this.barsLayer.destroyChildren();

    if (!this.barValues || this.barValues.length === 0)
      return;

    const stage = this.barsLayer.getStage();
    if (!stage) {
      this.pendingBarValues = values;
      return;
    }

    const width = stage.width();
    const height = stage.height();

    const minX = this.getMinX();
    const maxX = this.getMaxX();

    if (maxX === minX)
      return;

    const log = this.isLog;
    const tMin = toScale(minX, log);
    const tMax = toScale(maxX, log);
    const numBins = 20;
    const plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
    const plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;

    const binWidth = Math.max(1e-9, (tMax - tMin) / numBins);
    const bins = new Array(numBins).fill(0);

    (this.barValues ?? []).forEach((v) => {
      if (!Number.isFinite(v) || (log && v <= 0))
        return;
      const idx = Math.min(numBins - 1, Math.max(0, Math.floor((toScale(v, log) - tMin) / binWidth)));
      ++bins[idx];
    });

    const maxCount = Math.max(...bins) || 1;

    bins.forEach((count, i) => {
      const barX = EDITOR_PADDING.left + (i / numBins) * plotWidth;
      const barW = plotWidth / numBins - 1;
      const barH = (count / maxCount) * plotHeight;

      const rect = new _konva!.Rect({
        x: barX,
        y: EDITOR_PADDING.top + plotHeight - barH,
        width: barW,
        height: barH,
        fill: COLORS.barFill,
        opacity: 0.25,
        stroke: COLORS.barStroke,
        strokeWidth: 0.5,
      });

      this.barsLayer!.add(rect);
    });

    this.barsLayer!.batchDraw();
  }

  private enableCurveDrag() {
    if (!this.stage)
      return;

    let dragging = false;
    let startPointer: {x: number, y: number} | null = null;
    let startMean = 0;
    let startSigma = 0;
    let startX0 = 0;
    let startK = 0;

    this.stage.on('mousedown touchstart', (evt) => {
      if (this._prop.mode === 'freeform')
        return;
      if (!evt.evt)
        return;

      const pos = this.stage!.getPointerPosition();
      if (!pos)
        return;

      if (!this.isInPlotArea(pos, this._width, this._height))
        return;

      dragging = true;
      startPointer = pos;

      startMean = this._prop.mean ?? this.getDefaultMean();
      startSigma = this._prop.sigma ?? this.getDefaultSigma();
      startX0 = this._prop.x0 ?? this.getDefaultX0();
      startK = this._prop.k ?? this.getDefaultK();

      this.stage!.container().style.cursor = 'grabbing';
    });

    this.stage.on('mousemove touchmove', (evt) => {
      if (!dragging || !startPointer)
        return;

      const pos = this.stage!.getPointerPosition();
      if (!pos)
        return;

      const dx = pos.x - startPointer.x;
      const dy = pos.y - startPointer.y;

      const log = this.isLog;
      const tSpan = toScale(this.getMaxX(), log) - toScale(this.getMinX(), log);
      const plotWidth = this._width - EDITOR_PADDING.left - EDITOR_PADDING.right;
      const uPerPx = plotWidth === 0 ? 0 : tSpan / plotWidth;

      if (this._prop.mode === 'gaussian') {
        this._prop.mean = fromScale(toScale(log ? Math.max(this.getMinX(), startMean) : startMean, log) + dx * uPerPx, log);
        this._prop.sigma = Math.max(0.01, startSigma + (-dy) * this.dragScaleY * tSpan);
      }

      if (this._prop.mode === 'sigmoid') {
        this._prop.x0 = fromScale(toScale(log ? Math.max(this.getMinX(), startX0) : startX0, log) + dx * uPerPx, log);
        this._prop.k = Math.max(0.1, startK + (-dy) * this.dragScaleY * 50);
      }

      this._prop.line = this.computeLine();
      this.redrawFn?.();
      this.onParamsChanged?.(this._prop);
    });

    this.stage.on('mouseup touchend', () => {
      dragging = false;
      startPointer = null;
      this.stage!.container().style.cursor = 'default';
    });
  }

  private addSpecialHandle(width: number, height: number) {
    const clamp = (v: number, min: number, max: number) => Math.max(min, Math.min(max, v));

    if (this._prop.mode === 'freeform')
      return;

    const minX = this.getMinX();
    const maxX = this.getMaxX();
    const mapper = new CoordMapper(minX, maxX, width, height, !!this._prop.inverted, this.isLog);

    let x = (minX + maxX) / 2;
    let y = 0.5;

    if (this._prop.mode === 'gaussian') {
      x = this._prop.mean ?? this.getDefaultMean();
      y = 1;
    } else if (this._prop.mode === 'sigmoid') {
      x = this._prop.x0 ?? this.getDefaultX0();
      y = 0.5;
    }

    const coords = mapper.toCanvasCoords([x, y]);

    if (!this.specialHandle) {
      this.specialHandle = new _konva!.Circle({
        x: coords.x,
        y: coords.y,
        radius: 7,
        fill: COLORS.handle,
        draggable: true,
        hitStrokeWidth: 15,
      });

      this.specialHandle.on('dragmove', (evt) => {
        const m = new CoordMapper(this.getMinX(), this.getMaxX(), this._width, this._height, !!this._prop.inverted, this.isLog);
        const data = m.toDataCoords(evt.target.x(), evt.target.y());

        if (this._prop.mode === 'gaussian') {
          this._prop.mean = data.x;
          this._prop.sigma = Math.max(0.01, Math.abs(data.y - 1));
        } else if (this._prop.mode === 'sigmoid') {
          this._prop.x0 = data.x;
          this._prop.k = clamp(Math.abs(data.y - 0.5) * 30, 0.1, 30);
        }

        this._prop.line = this.computeLine();
        this.konvaLine!.points(this._prop.line.flatMap((p) => {
          const c = m.toCanvasCoords([p[0], p[1]]);
          return [c.x, c.y];
        }));
        this.layer!.batchDraw();
        this.onParamsChanged?.(this._prop);
      });

      this.specialHandle.on('dragend', () => this.redrawFn?.());

      this.layer!.add(this.specialHandle);
    } else {
      if (!this.specialHandle.getLayer())
        this.layer!.add(this.specialHandle);
      this.specialHandle.position(coords);
    }
    this.layer!.batchDraw();
  }

  private showPointEditor(dataIndex: number, clientX: number, clientY: number): void {
    const point = this._prop.line[dataIndex];
    const xInput = ui.input.float('X', {value: point[0], min: this.getMinX(), max: this.getMaxX(), format: '#0.00', step: 0.01});
    const yInput = ui.input.float('Y', {value: this._prop.inverted ? 1 - point[1] : point[1], min: 0, max: 1, format: '#0.00', step: 0.01});

    const close = () => {
      content.removeEventListener('keydown', onKey);
      popup.remove();
    };

    const apply = () => {
      const x = xInput.value;
      const y = yInput.value;
      if (x == null || y == null || isNaN(x) || isNaN(y))
        return;
      this._prop.line[dataIndex] = [x, this._prop.inverted ? 1 - y : y];
      this.redrawAll();
    };

    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Enter') {
        apply();
        close();
      }
      if (e.key === 'Escape')
        close();
    };

    const content = ui.inputs([xInput, yInput]);
    content.style.overflow = 'hidden';
    content.addEventListener('keydown', onKey);

    const rootRect = this.root.getBoundingClientRect();
    const popup = ui.showPopup(content, this.root, {dx: clientX - rootRect.left, dy: clientY - rootRect.bottom, smart: false});
  }

  setColumn(col: DG.Column | null): void {
    if (!col)
      return;

    const values = col.toList();
    this.drawBars(values);
  }
}
