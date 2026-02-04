import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {TreeHelper} from '../../utils/tree-helper';
import {CanvasTreeRenderer} from './canvas-tree-renderer';
import {ITreeStyler, MarkupNodeType} from './markup';
import {GridTreePlacer} from './grid-tree-placer';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';

export abstract class GridTreeRendererBase<TNode extends MarkupNodeType> extends CanvasTreeRenderer<TNode> {
  protected _leftPadding: number = 6;
  protected _rightPadding: number = 6;

  private readonly th: ITreeHelper;

  protected readonly grid: DG.Grid;

  get leftPadding(): number { return this._leftPadding; }

  get rightPadding(): number { return this._rightPadding; }

  /** treeRoot can be null in case of the grid.dataFrame.rowCount is zero
   * @param {DG.Grid}grid - grid to render tree for
   * @param {TNode}treeRoot - tree root node
   * @param {GridTreePlacer<TNode>}placer - tree placer
   * @param {ITreeStyler<TNode>}mainStyler - main tree styler
   * @param {ITreeStyler<TNode>}lightStyler - light tree styler
   * @param {ITreeStyler<TNode>}currentStyler - current tree styler
   * @param {ITreeStyler<TNode>}mouseOverStyler - mouse over tree styler
   * @param {ITreeStyler<TNode>}selectionStyler - selection tree styler*/
  protected constructor(
    grid: DG.Grid, treeRoot: TNode | null, placer: GridTreePlacer<TNode>,
    mainStyler: ITreeStyler<TNode>, lightStyler: ITreeStyler<TNode>,
    currentStyler: ITreeStyler<TNode>, mouseOverStyler: ITreeStyler<TNode>, selectionStyler: ITreeStyler<TNode>,
  ) {
    super(treeRoot, placer,
      mainStyler, lightStyler,
      currentStyler, mouseOverStyler, selectionStyler);

    this.th = new TreeHelper();
    this.grid = grid;

    this.subs.push(this.grid.onBeforeDrawContent.subscribe(this.gridOnChanged.bind(this)));
    this.subs.push(ui.onSizeChanged(this.grid.root).subscribe(this.gridOnChanged.bind(this)));

    this.gridOnChanged();
  }

  // -- Handle events --

  /** Override to prevent change canvas size at {@link CanvasTreeRenderer.viewOnSizeChanged() }*/
  protected override viewOnSizeChanged() {
    // super.viewOnSizeChanged();
  }

  protected gridOnChanged() {
    if (!this.view || !this.canvas) return;

    // view is GridNeighbor
    const cw: number = this.view.clientWidth;
    const ch: number = this.view.clientHeight - this.grid.colHeaderHeight;

    this.canvas.width = cw * window.devicePixelRatio;
    this.canvas.height = ch * window.devicePixelRatio;

    // const gridRowCount: number = this.grid.vertScroll.max - this.grid.vertScroll.min;
    // const gridRowHeight: number = this.grid.props.rowHeight;

    this.canvas.style.position = 'absolute';
    this.canvas.style.left = `${0}px`;
    this.canvas.style.top = `${this.grid.colHeaderHeight}px`;
    this.canvas.style.height = `${ch}px`;
    this.canvas.style.width = `${cw}px`;

    this.render('gridOnChanged');
  }

  override canvasOnWheel(e: WheelEvent): void {
    // console.debug('Dendrogram: GridTreeRendererBase.canvasOnWheel()');
    // Intercept wheel event handling to prevent CanvasTreeRender.canvasOnWheel() handler zooming the tree
    //super.canvasOnWheel(e);
    e.preventDefault();

    if (!this.canvas) return;

    // @ts-ignore // for wheelDelta property
    const delta = 5 * e.wheelDelta / -168;

    if (e.ctrlKey || e.metaKey) {
      // Zooming
      // const oldZoomFactor = this.zoomFactor;
      this.xZoomFactor -= delta * 0.2;
      // console.log(delta);
      this.xZoomFactor = Math.min(Math.max(1, this.xZoomFactor), 100);
      this.gridOnChanged();
      return;
    }
    //TODO: Use RangeSlider.scrollBy() method
    //this.grid.vertScroll.scrollBy(delta);
    this.grid.vertScroll.scrollTo(this.grid.vertScroll.min + delta);
  }


  override canvasOnMouseDown(e: MouseEvent): void {
    // console.debug('Dendrogram: GridTreeRendererBase.canvasOnMouseDown()');
    // Intercept to prevent handling drag mode
    // super.canvasOnMouseDown();
    e.preventDefault();
  }

  override canvasOnMouseUp(e: MouseEvent): void {
    // console.debug('Dendrogram: GridTreeRendererBase.canvasOnMouseUp()');
    // Intercept to prevent handling drag mode
    // super.canvasOnMouseDown();
    e.preventDefault();
  }

  override canvasOnMouseMove(e: MouseEvent): void {
    // console.debug('Dendrogram: GridTreeRendererBase.canvasOnMouseMove()');

    super.canvasOnMouseMove(e); // super handler is required to handle mouse over on nodes
  }
}

