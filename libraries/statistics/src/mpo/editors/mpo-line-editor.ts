/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import Konva from 'konva';
import {DesirabilityLine, NumericalDesirability} from '../mpo';
import {Subject} from 'rxjs';

type Point = [number, number];

// Constants for the editor layout
const EDITOR_PADDING = {top: 10, right: 10, bottom: 20, left: 30};
const POINT_RADIUS = 3;

const COLORS = {
  line: '#2077b4',
  point: '#d72f30',
  barFill: 'rgba(160,196,255,0.5)',
};

class CoordMapper {
  private plotWidth: number;
  private plotHeight: number;
  private scaleX: number;
  private scaleY: number;

  constructor(
    private minX: number,
    private maxX: number,
    private width: number,
    private height: number,
  ) {
    this.plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
    this.plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;
    this.scaleX = (this.maxX - this.minX === 0) ? 1 : this.plotWidth / (this.maxX - this.minX);
    this.scaleY = this.plotHeight;
  }

  toCanvasCoords(p: Point): {x: number, y: number} {
    const canvasX = EDITOR_PADDING.left + (p[0] - this.minX) * this.scaleX;
    const canvasY = EDITOR_PADDING.top + this.plotHeight - (p[1] * this.scaleY);

    return {x: canvasX, y: canvasY};
  }

  toDataCoords(canvasX: number, canvasY: number): {x: number, y: number} {
    let dataX = this.minX + (canvasX - EDITOR_PADDING.left) / this.scaleX;
    let dataY = (EDITOR_PADDING.top + this.plotHeight - canvasY) / this.scaleY;

    dataX = Math.max(this.minX, Math.min(this.maxX, dataX));
    dataY = Math.max(0, Math.min(1, dataY));

    return {x: dataX, y: dataY};
  }
}

export class MpoDesirabilityLineEditor {
  root = ui.div();
  onChanged = new Subject<DesirabilityLine>();
  supportsModeDialog: boolean = true;

  private _prop: NumericalDesirability;
  private barsLayer: Konva.Layer;
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

  constructor(prop: NumericalDesirability, width: number, height: number) {
    this._prop = prop;
    this.ensureDefaultLine();
    this.barsLayer = new Konva.Layer();
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

  private isInPlotArea(pos: {x: number, y: number}, width: number, height: number) {
    return (
      pos.x >= EDITOR_PADDING.left &&
      pos.x <= width - EDITOR_PADDING.right &&
      pos.y >= EDITOR_PADDING.top &&
      pos.y <= height - EDITOR_PADDING.bottom
    );
  }

  private initKonva(width: number, height: number) {
    if (!this.root.parentElement) {
      console.warn('Konva container not attached to DOM yet.');
      return;
    }

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

      const mapper = new CoordMapper(minX, maxX, width, height);
      const sortedIndices = [...this._prop.line.keys()]
        .sort((a, b) => this._prop.line[a][0] - this._prop.line[b][0]);
      const sortedLine = sortedIndices.map((idx) => this._prop.line[idx]);

      sortedLine.forEach((p, index) => {
        const coords = mapper.toCanvasCoords([p[0], p[1]]);
        konvaPoints.push(coords.x, coords.y);

        const pointCircle = new Konva.Circle({
          x: coords.x,
          y: coords.y,
          radius: POINT_RADIUS,
          fill: COLORS.point,
          stroke: 'black',
          strokeWidth: 1,
          draggable: this._prop.mode === 'freeform',
          hitStrokeWidth: 5,
        });

        // Store index directly on the node for easy access
        pointCircle.setAttr('_dataIndex', sortedIndices[index]);

        pointCircle.on('dragmove', (evt) => {
          if (this._prop.mode !== 'freeform')
            return;

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
          const plotBottom = height - EDITOR_PADDING.bottom;
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
          if (this._prop.mode !== 'freeform')
            return;

          this._prop.line.sort((a, b) => a[0] - b[0]);
          this.redrawFn?.();
        });

        pointCircle.on('contextmenu', (evt) => {
          if (this._prop.mode !== 'freeform')
            return;

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
          if (this._prop.mode !== 'freeform')
            return;

          this.stage!.container().style.cursor = 'pointer';
          const circle = evt.target as Konva.Circle;
          const pos = circle.position();
          const dataCoords = mapper.toDataCoords(pos.x, pos.y);
          const tooltipText = `X: ${dataCoords.x.toFixed(2)}, Y: ${dataCoords.y.toFixed(2)}<br><br>Drag to move, right-click to delete`;
          ui.tooltip.show(tooltipText, evt.evt.clientX, evt.evt.clientY);
        });

        pointCircle.on('mouseleave', () => {
          if (this._prop.mode !== 'freeform')
            return;

          this.stage!.container().style.cursor = 'default';
          ui.tooltip.hide();
        });

        this.pointsGroup!.add(pointCircle);
      });

      this.konvaLine!.points(konvaPoints);
      this.layer!.batchDraw();

      // --- Add special handle (Gaussian peak / Sigmoid inflection) ---
      this.addSpecialHandle(width, height);

      if (notify && this._prop.mode === 'freeform')
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

      const mapper = new CoordMapper(minX, maxX, this.stage!.width(), this.stage!.height());
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
        if (this.isInPlotArea(pos, width, height))
          this.stage.container().style.cursor = 'grab';
        else
          this.stage.container().style.cursor = 'default';
      }
    });

    this.stage.on('mouseout', () => ui.tooltip.hide());

    // Enable curve drag (smooth)
    this.enableCurveDrag(width, height);

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

    const minX = this.getMinX();
    const maxX = this.getMaxX();
    const n = 60;
    const line: DesirabilityLine = [];

    for (let i = 0; i <= n; i++) {
      const x = minX + (maxX - minX) * (i / n);
      let y = 0;

      if (this._prop.mode === 'gaussian') {
        this._prop.mean ??= (minX + maxX) / 2;
        this._prop.sigma ??= (maxX - minX) / 6;
        const mean = this._prop.mean;
        const sigma = this._prop.sigma;
        const z = (x - mean) / sigma;
        y = Math.exp(-0.5 * z * z);
      }

      if (this._prop.mode === 'sigmoid') {
        this._prop.x0 ??= (minX + maxX) / 2;
        this._prop.k ??= 10;
        const x0 = this._prop.x0;
        const k = this._prop.k;
        y = 1 / (1 + Math.exp(-k * (x - x0)));
      }

      line.push([x, y]);
    }

    return line;
  }

  private drawAxes(minX: number, maxX: number, width: number, height: number) {
    this.layer!.add(
      new Konva.Line({
        points: [EDITOR_PADDING.left, height - EDITOR_PADDING.bottom, width - EDITOR_PADDING.right, height - EDITOR_PADDING.bottom],
        stroke: 'grey',
        strokeWidth: 1,
      }),
      new Konva.Line({
        points: [EDITOR_PADDING.left, EDITOR_PADDING.top, EDITOR_PADDING.left, height - EDITOR_PADDING.bottom],
        stroke: 'grey',
        strokeWidth: 1,
      }),
      new Konva.Text({x: EDITOR_PADDING.left, y: height - EDITOR_PADDING.bottom + 3, text: minX.toFixed(1), fontSize: 9, fill: 'grey'}),
      new Konva.Text({x: width - EDITOR_PADDING.right - 15, y: height - EDITOR_PADDING.bottom + 3, text: maxX.toFixed(1), fontSize: 9, fill: 'grey'}),
    );
  }

  getMinX(): number {
    return this._prop.min ?? Math.min(...this._prop.line.map((p) => p[0])) ?? 0;
  }

  getMaxX(): number {
    return this._prop.max ?? Math.max(...this._prop.line.map((p) => p[0])) ?? 1;
  }

  getDefaultMean(): number {
    return (this.getMinX() + this.getMaxX()) / 2;
  }

  getDefaultSigma(): number {
    return Math.max(0.01, (this.getMaxX() - this.getMinX()) / 6);
  }

  getDefaultX0(): number {
    return (this.getMinX() + this.getMaxX()) / 2;
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

    this.layer.destroyChildren();

    this.drawAxes(minX, maxX, width, height);

    this.layer.add(this.konvaLine!, this.pointsGroup!);

    this.redrawFn(notify);

    if (this.barValues)
      this.drawBars();
  }

  drawBars(values?: number[]) {
    if (!this.barsLayer)
      return;
    this.barsLayer.destroyChildren();

    if (values)
      this.barValues = values;

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

    const numBins = 20;
    const plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
    const plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;

    const binWidth = Math.max(1e-9, (maxX - minX) / numBins);
    const bins = new Array(numBins).fill(0);

    (this.barValues ?? []).forEach((v) => {
      const idx = Math.min(Math.floor((v - minX) / binWidth), numBins - 1);
      bins[idx]++;
    });

    const maxCount = Math.max(...bins) || 1;

    bins.forEach((count, i) => {
      const barX = EDITOR_PADDING.left + (i * binWidth / (maxX - minX)) * plotWidth;
      const barW = (binWidth / (maxX - minX)) * plotWidth - 1;
      const barH = (count / maxCount) * plotHeight;

      const rect = new Konva.Rect({
        x: barX,
        y: EDITOR_PADDING.top + plotHeight - barH,
        width: barW,
        height: barH,
        fill: COLORS.barFill,
        stroke: COLORS.line,
        strokeWidth: 0.5,
      });

      this.barsLayer.add(rect);
    });

    this.barsLayer.batchDraw();
  }

  private enableCurveDrag(width: number, height: number) {
    if (!this.stage)
      return;

    let dragging = false;
    let startPointer: {x: number, y: number} | null = null;
    let startMean = 0;
    let startSigma = 0;
    let startX0 = 0;
    let startK = 0;

    const plotWidth = width - EDITOR_PADDING.left - EDITOR_PADDING.right;
    const plotHeight = height - EDITOR_PADDING.top - EDITOR_PADDING.bottom;

    const scaleX = (this.getMaxX() - this.getMinX()) / plotWidth;
    const scaleY = 1 / plotHeight;

    this.stage.on('mousedown touchstart', (evt) => {
      if (this._prop.mode === 'freeform')
        return;
      if (!evt.evt)
        return;

      const pos = this.stage!.getPointerPosition();
      if (!pos)
        return;

      if (!this.isInPlotArea(pos, width, height))
        return;

      dragging = true;
      startPointer = pos;

      startMean = this._prop.mean ?? (this.getMinX() + this.getMaxX()) / 2;
      startSigma = this._prop.sigma ?? (this.getMaxX() - this.getMinX()) / 6;
      startX0 = this._prop.x0 ?? (this.getMinX() + this.getMaxX()) / 2;
      startK = this._prop.k ?? 10;

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

      if (this._prop.mode === 'gaussian') {
        this._prop.mean = startMean + dx * scaleX;
        this._prop.sigma = Math.max(0.01, startSigma + (-dy) * scaleY * (this.getMaxX() - this.getMinX()));
      }

      if (this._prop.mode === 'sigmoid') {
        this._prop.x0 = startX0 + dx * scaleX;
        this._prop.k = Math.max(0.1, startK + (-dy) * scaleY * 50);
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
    const mapper = new CoordMapper(minX, maxX, width, height);

    let x = (minX + maxX) / 2;
    let y = 0.5;

    if (this._prop.mode === 'gaussian') {
      x = this._prop.mean ?? (minX + maxX) / 2;
      y = 1;
    } else if (this._prop.mode === 'sigmoid') {
      x = this._prop.x0 ?? (minX + maxX) / 2;
      y = 0.5;
    }

    const coords = mapper.toCanvasCoords([x, y]);

    if (!this.specialHandle) {
      this.specialHandle = new Konva.Circle({
        x: coords.x,
        y: coords.y,
        radius: 7,
        fill: 'orange',
        stroke: 'black',
        strokeWidth: 1,
        draggable: true,
        hitStrokeWidth: 15,
      });

      this.specialHandle.on('dragmove', (evt) => {
        const pos = evt.target.position();
        const data = mapper.toDataCoords(pos.x, pos.y);

        if (this._prop.mode === 'gaussian') {
          this._prop.mean = data.x;
          this._prop.sigma = Math.max(0.01, Math.abs(data.y - 1));
        } else if (this._prop.mode === 'sigmoid') {
          this._prop.x0 = data.x;
          this._prop.k = clamp(Math.abs(data.y - 0.5) * 30, 0.1, 30);
        }

        this._prop.line = this.computeLine();
        this.redrawFn?.();
        this.onParamsChanged?.(this._prop);
      });

      this.layer!.add(this.specialHandle);
    } else {
      if (!this.specialHandle.getLayer())
        this.layer!.add(this.specialHandle);
      this.specialHandle.position(coords);
    }
    this.layer!.batchDraw();
  }

  setColumn(col: DG.Column | null): void {
    if (!col)
      return;

    const values = col.toList();
    this.drawBars(values);
  }
}
