import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable, Subject, Unsubscribable} from 'rxjs';

import {ITreeHelper, NodeType} from '@datagrok-libraries/bio';

import {TreeHelper} from '../utils/tree-helper';
import {TreeRendererBase} from './tree-renderers/tree-renderer-base';
import {CanvasTreeRenderer} from './tree-renderers/canvas-tree-renderer';
import {MarkupNodeType} from './tree-renderers/markup';

// export type GridTreeRendererEventArgsType<TNode extends NodeType> = {
//   target: GridTreeRendererBase<TNode>,
//   context: CanvasRenderingContext2D,
//   lengthRatio: number,
// };

export abstract class GridTreeRendererBase<TNode extends MarkupNodeType>
  extends CanvasTreeRenderer<TNode> {
  protected _leftPadding: number = 6;
  protected _rightPadding: number = 6;

  private readonly th: ITreeHelper;

  protected readonly _grid: DG.Grid;

  get grid(): DG.Grid { return this._grid; }

  get leftPadding(): number { return this._leftPadding; }

  get rightPadding(): number { return this._rightPadding; }

  protected constructor(treeRoot: TNode, totalLength: number, view: HTMLElement, grid: DG.Grid) {
    super(treeRoot, totalLength, view);

    this.th = new TreeHelper();
    this._grid = grid;
    this._totalLength = totalLength;
    this.subs.push(this.grid.onBeforeDrawContent.subscribe(this.render.bind(this)));
    this.subs.push(ui.onSizeChanged(this.grid.root).subscribe(this.render.bind(this)));
  }
}

