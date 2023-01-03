import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {ITreeHelper, NodeType} from '@datagrok-libraries/bio';

import {TreeHelper} from '../../utils/tree-helper';
import {CanvasTreeRenderer} from './canvas-tree-renderer';
import {ITreeStyler, MarkupNodeType} from './markup';
import {GridTreePlacer} from './grid-tree-placer';
import {Mouse} from 'puppeteer';

export abstract class GridTreeRendererBase<TNode extends MarkupNodeType> extends CanvasTreeRenderer<TNode> {
  protected _leftPadding: number = 6;
  protected _rightPadding: number = 6;

  private readonly th: ITreeHelper;

  protected readonly grid: DG.Grid;

  get leftPadding(): number { return this._leftPadding; }

  get rightPadding(): number { return this._rightPadding; }

  protected constructor(
    grid: DG.Grid, treeRoot: TNode, placer: GridTreePlacer<TNode>,
    mainStyler: ITreeStyler<TNode>, lightStyler: ITreeStyler<TNode>,
    currentStyler: ITreeStyler<TNode>, mouseOverStyler: ITreeStyler<TNode>, selectionStyler: ITreeStyler<TNode>
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

  canvasOnWheel(e: WheelEvent): void {
    e.preventDefault();
  }

  protected override canvasOnMouseDown(e: MouseEvent): void {
    e.preventDefault();
  }

  protected override canvasOnMouseUp(e: MouseEvent): void {
    e.preventDefault();
  }
}

