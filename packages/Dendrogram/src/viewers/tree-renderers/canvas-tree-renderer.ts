import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import * as rxjs from 'rxjs';
import {TreeRendererBase} from './tree-renderer-base';
import {ITreeStyler, MarkupNodeType, renderNode, TraceTargetType, TreeStylerBase} from './markup';
import {RectangleTreeHoverType, RectangleTreePlacer} from './rectangle-tree-placer';
import {TreeHelper} from '../../utils/tree-helper';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

function canvasToTreePoint<TNode extends MarkupNodeType>(
  canvasPoint: DG.Point, canvas: HTMLCanvasElement, placer: RectangleTreePlacer<TNode>
): DG.Point {
  const canvasWidth = canvas.clientWidth - placer.padding.left - placer.padding.right;
  const res: DG.Point = new DG.Point(
    (canvasPoint.x - placer.padding.left) * placer.totalLength / canvasWidth,
    placer.top + canvasPoint.y * placer.height / canvas.clientHeight);
  return res;
}

function treeToCanvasPoint<TNode extends MarkupNodeType>(
  treePoint: DG.Point, canvas: HTMLCanvasElement, placer: RectangleTreePlacer<TNode>
): DG.Point {
  const canvasWidth = canvas.clientWidth - placer.padding.left - placer.padding.right;
  const res: DG.Point = new DG.Point(
    placer.padding.left + canvasWidth * treePoint.x / placer.totalLength,
    canvas.clientHeight * (treePoint.y - placer.top) / placer.height);
  return res;
}

export class CanvasTreeRenderer<TNode extends MarkupNodeType>
  extends TreeRendererBase<TNode, RectangleTreeHoverType<TNode>> {
  protected canvas?: HTMLCanvasElement;

  protected readonly placer: RectangleTreePlacer<TNode>;

  protected readonly lightStyler: ITreeStyler<TNode>;
  protected readonly currentStyler: ITreeStyler<TNode>;
  protected readonly mouseOverStyler: ITreeStyler<TNode>;
  protected readonly selectionStyler: ITreeStyler<TNode>;

  protected _mainStyler: ITreeStyler<TNode>;
  protected _mainStylerOnChangedSub!: rxjs.Unsubscribable;

  get mainStyler(): ITreeStyler<TNode> { return this._mainStyler; }

  set mainStyler(value: ITreeStyler<TNode>) {
    if (this.view)
      this._mainStylerOnChangedSub!.unsubscribe();

    this._mainStyler = value;

    if (this.view) {
      this._mainStylerOnChangedSub = this._mainStyler.onStylingChanged.subscribe(this.stylerOnChanged.bind(this));
      this.render('mainStyler');
    }
  }

  constructor(
    treeRoot: TNode, placer: RectangleTreePlacer<TNode>,
    mainStyler: ITreeStyler<TNode>, lightStyler: ITreeStyler<TNode>,
    currentStyler: ITreeStyler<TNode>, mouseOverStyler: ITreeStyler<TNode>, selectionStyler: ITreeStyler<TNode>
  ) {
    super(treeRoot);

    this.placer = placer;
    this.placer.onPlacingChanged.subscribe(this.placerOnChanged.bind(this));

    this._mainStyler = mainStyler;
    this.lightStyler = lightStyler;
    this.currentStyler = currentStyler;
    this.mouseOverStyler = mouseOverStyler;
    this.selectionStyler = selectionStyler;

    // this.view onSizeChanged event is handled with ancestor TreeRenderBase calls this.render()
  }

  private renderCounter: number = 0;

  /**
   * Leaves along axis has step 1 for every leaf.
   * Draw tree along all canvas height in terms of leaves' axis placement.
   */
  render(purpose: string): void {
    if (!this.view || !this.canvas) return;

    const ctx = this.canvas.getContext('2d')!;
    const dpr: number = window.devicePixelRatio;

    // Here we will render range of leaves
    const plotWidth: number = this.canvas.width - dpr * (this.placer.padding.left + this.placer.padding.right);
    const plotHeight: number = this.canvas.height;
    const placerHeight: number = this.placer.bottom - this.placer.top;
    const lengthRatio: number = plotWidth / this.placer.totalLength; // px/[length unit]
    const stepRatio: number = plotHeight / placerHeight; // px/[step unit, row]

    const t1: number = Date.now();
    ctx.save();
    try {
      (function clearNodeDesc(node: TNode) {
        node.desc = '';
        for (const childNode of (node.children ?? []))
          clearNodeDesc(childNode as TNode);
      })(this.treeRoot);

      ctx.fillStyle = '#FFFFFF';
      ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);

      // // large red diagonal cross
      // ctx.beginPath();
      // ctx.strokeStyle = '#800000';
      // ctx.moveTo(0, 0);
      // ctx.lineTo(ctx.canvas.width, ctx.canvas.height);
      // ctx.moveTo(0, ctx.canvas.height);
      // ctx.lineTo(ctx.canvas.width, 0);
      // ctx.stroke();

      const invisibleStyler = new TreeStylerBase<TNode>('invisible', 0, 0, false, '#00000000', '#00000000');

      this.renderCounter++;
      console.debug(`*** ${this.renderCounter} Dendrogram: CanvasTreeRenderer.render(), ` +
        `main & light, traceback hover & selection, ` +
        `purpose '${purpose}'.`);
      const styler: ITreeStyler<TNode> = !this.mouseOver ? this._mainStyler : this.lightStyler;
      const selectionTraceList: TraceTargetType<TNode>[] = this.selections.map(
        (sel) => { return {target: sel.node, styler: this.selectionStyler}; });
      renderNode(
        {
          ctx: ctx, firstRowIndex: this.placer.top, lastRowIndex: this.placer.bottom,
          leftPadding: this.placer.padding.left, lengthRatio: lengthRatio, stepRatio: stepRatio,
          totalLength: this.placer.totalLength, styler: styler,
        },
        this.treeRoot, 0, [...selectionTraceList]);

      for (const selection of this.selections) {
        renderNode(
          {
            ctx: ctx, firstRowIndex: this.placer.top, lastRowIndex: this.placer.bottom,
            leftPadding: this.placer.padding.left, lengthRatio: lengthRatio, stepRatio: stepRatio,
            styler: this.selectionStyler, totalLength: this.placer.totalLength
          },
          selection.node, selection.nodeHeight, []);
      }

      if (this.current) {
        const currentTraceList: TraceTargetType<TNode>[] = [{target: this.current.node, styler: this.currentStyler}];
        renderNode(
          {
            ctx: ctx, firstRowIndex: this.placer.top, lastRowIndex: this.placer.bottom,
            leftPadding: this.placer.padding.left, lengthRatio: lengthRatio, stepRatio: stepRatio,
            totalLength: this.placer.totalLength, styler: invisibleStyler,
          },
          this.treeRoot, 0, [...currentTraceList]);

        // children
        // renderNode(ctx, this.current.node,
        //   this.placer.top, this.placer.bottom,
        //   this.placer.padding.left, lengthRatio, stepRatio, this.currentStyler,
        //   this.placer.totalLength, this.current.nodeHeight,
        //   []);
      }

      if (this.mouseOver) {
        const mouseOverTraceList: TraceTargetType<TNode>[] = [
          {target: this.mouseOver.node, styler: this.mouseOverStyler}];
        renderNode(
          {
            ctx: ctx, firstRowIndex: this.placer.top, lastRowIndex: this.placer.bottom,
            leftPadding: this.placer.padding.left, lengthRatio: lengthRatio, stepRatio: stepRatio,
            totalLength: this.placer.totalLength, styler: invisibleStyler,
          },
          this.treeRoot, 0, [...mouseOverTraceList]);

        // children
        renderNode(
          {
            ctx: ctx, firstRowIndex: this.placer.top, lastRowIndex: this.placer.bottom,
            leftPadding: this.placer.padding.left, lengthRatio: lengthRatio, stepRatio: stepRatio,
            totalLength: this.placer.totalLength, styler: this.mouseOverStyler,
          },
          this.mouseOver.node, this.mouseOver.nodeHeight,
          []);
      }

      console.debug('');
      console.debug('');
    } catch (err: any) {
      errorToConsole(err);
      throw err;
    } finally {
      ctx.restore();
      const t2: number = Date.now();
      console.debug('Dendrogram: CanvasTreeRenderer.render(), ' + `ET: ${((t2 - t1) / 1000).toFixed(3)}`);
      this._onAfterRender.next({target: this, context: ctx, lengthRatio});
    }
  }

  // -- View --

  public override attach(view: HTMLElement): void {
    super.attach(view);
    this.canvas = ui.canvas();
    this.view!.appendChild(this.canvas);

    this._mainStylerOnChangedSub = this._mainStyler.onStylingChanged.subscribe(this.stylerOnChanged.bind(this));
    this.subs.push(this.lightStyler.onStylingChanged.subscribe(this.stylerOnChanged.bind(this)));
    this.subs.push(this.currentStyler.onStylingChanged.subscribe(this.stylerOnChanged.bind(this)));
    this.subs.push(this.mouseOverStyler.onStylingChanged.subscribe(this.stylerOnChanged.bind(this)));
    this.subs.push(this.selectionStyler.onStylingChanged.subscribe(this.stylerOnChanged.bind(this)));

    this.subs.push(rxjs.fromEvent<WheelEvent>(this.canvas, 'wheel').subscribe(this.canvasOnWheel.bind(this)));
    this.subs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'mousedown').subscribe(this.canvasOnMouseDown.bind(this)));
    this.subs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'mouseup').subscribe(this.canvasOnMouseUp.bind(this)));
    this.subs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'mousemove').subscribe(this.canvasOnMouseMove.bind(this)));
    this.subs.push(rxjs.fromEvent<MouseEvent>(this.canvas, 'click').subscribe(this.canvasOnClick.bind(this)));
  }

  public override detach(): void {
    this.canvas!.remove();
    delete this.canvas;
    this._mainStylerOnChangedSub.unsubscribe();
    super.detach();
  }

  // -- Handle events --

  /** Changes canvas size and canvas outer/style size according to {@link view} size changed */
  protected override viewOnSizeChanged() {
    if (!this.view || !this.canvas) return;
    const dpr: number = window.devicePixelRatio;

    const cw: number = this.view.clientWidth;
    const ch: number = this.view.clientHeight;

    this.canvas.width = dpr * cw;
    this.canvas.height = dpr * ch;

    this.canvas.style.left = `{0}px`;
    this.canvas.style.top = `{0}px`;
    this.canvas.style.width = `${cw}px`;
    this.canvas.style.height = `${ch}px`;

    super.viewOnSizeChanged(); // calls this.render()
  }

  private placerOnChanged() {
    this.render('placerOnChanged()');
  }

  private stylerOnChanged() {
    this.render('stylerOnChanged()');
  }

  protected canvasOnWheel(e: WheelEvent): void {
    if (!this.canvas || this.mouseDragging) return;
    e.preventDefault();

    const pos = canvasToTreePoint(new DG.Point(e.offsetX, e.offsetY), this.canvas, this.placer);

    // @ts-ignore
    const delta: number = e.wheelDelta / -168;
    const newTop = pos.y - (pos.y - this.placer.top) * (1 + 0.2 * delta);
    const newBottom = pos.y + (this.placer.bottom - pos.y) * (1 + 0.2 * delta);

    this.placer.update({top: newTop, bottom: newBottom});

    // e.stopPropagation();
  }

  protected mouseDragging: { pos: DG.Point, top: number, bottom: number } | null = null;

  protected canvasOnMouseDown(e: MouseEvent): void {
    if (!this.view || !this.canvas) return;

    const pos = canvasToTreePoint(new DG.Point(e.offsetX, e.offsetY), this.canvas, this.placer);
    this.mouseDragging = {pos: pos, top: this.placer.top, bottom: this.placer.bottom};
  }

  protected canvasOnMouseUp(e: MouseEvent): void {
    this.mouseDragging = null;
  }

  protected canvasOnMouseMove(e: MouseEvent): void {
    if (!this.view || !this.canvas) return;

    const canvasPoint = new DG.Point(e.offsetX, e.offsetY);

    const md = this.mouseDragging;
    if (md) {
      const mousePosY: number = md.top + e.offsetY * this.placer.height / this.canvas.clientHeight;
      const deltaPosY = md.pos.y - mousePosY;

      this.placer.update({
        top: md.top + deltaPosY,
        bottom: md.bottom + deltaPosY
      });
    } else {
      // console.debug('CanvasTreeRender.onMouseMove() --- getNode() ---');
      this.mouseOver = this.placer.getNode(
        this.treeRoot, canvasPoint, this._mainStyler.lineWidth, this._mainStyler.nodeSize,
        (canvasP: DG.Point): DG.Point => { return treeToCanvasPoint(canvasP, this.canvas!, this.placer); });

      this._mainStyler.fireTooltipShow(this.mouseOver ? this.mouseOver.node : null, e);
    }
  }

  protected canvasOnClick(e: MouseEvent): void {
    if (e.button == 0) {
      if (e.ctrlKey) {
        if (this.mouseOver) {
          const selections = [...this.selections];
          let deselectedIdx = -1;
          let replaceTo = -1;
          const th = new TreeHelper();
          for (let selI = 0; selI < this.selections.length; selI++) {
            const selection = this.selections[selI];
            if (selection.node == this.mouseOver!.node)
              deselectedIdx = selI;

            if (th.includes(selection.node, this.mouseOver.node)) {
              // Do nothing, because the clicked node (this.mouseOver.node) is sub of already selected
            }
            if (th.includes(this.mouseOver.node, selection.node))
              replaceTo = selI;
          }
          if (deselectedIdx != -1)
            selections.splice(deselectedIdx, 1);
          else if (replaceTo != -1)
            selections[replaceTo] = this.mouseOver;
          else
            selections.push(this.mouseOver);

          this.selections = selections;
        }
      } else {
        this.current = this.mouseOver;
      }
    }
  }
}
