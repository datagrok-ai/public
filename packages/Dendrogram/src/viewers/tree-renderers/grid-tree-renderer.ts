import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {NodeType} from '@datagrok-libraries/bio';

import {ITreeStyler, markupNode, MarkupNodeType, TreeStylerBase} from './markup';
import {toRgba, setAlpha} from '@datagrok-libraries/utils/src/color';
import {GridTreeRendererBase} from './grid-tree-renderer-base';
import {LINE_WIDTH, NODE_SIZE, TreeColorNames, TreeDefaultPalette} from '../dendrogram';
import {GridTreePlacer} from './grid-tree-placer';

const TRANS_ALPHA = 0.7;

/** Draws only nodes/leaves visible in leaf range */
export class LeafRangeGridTreeRenderer extends GridTreeRendererBase<MarkupNodeType> {
  constructor(
    grid: DG.Grid, tree: MarkupNodeType, placer: GridTreePlacer<MarkupNodeType>,
    mainStyler: ITreeStyler<MarkupNodeType>, lightStyler: ITreeStyler<MarkupNodeType>,
    currentStyler: ITreeStyler<MarkupNodeType>, mouseOverStyler: ITreeStyler<MarkupNodeType>,
    selectionStyler: ITreeStyler<MarkupNodeType>
  ) {
    super(grid, tree, placer, mainStyler, lightStyler, currentStyler, mouseOverStyler, selectionStyler);
  }

  // -- View --

  public override attach(view: HTMLElement) {
    super.attach(view);

    this.view!.style.setProperty('overflow-y', 'hidden', 'important');
    this.canvas!.style.position = 'absolute';
  }

  // --

  public static create(
    grid: DG.Grid, tree: NodeType, placer: GridTreePlacer<MarkupNodeType>
  ): GridTreeRendererBase<MarkupNodeType> {
    const mainStyler = new TreeStylerBase<MarkupNodeType>('main',
      LINE_WIDTH, NODE_SIZE, true,
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Main], TRANS_ALPHA)),
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Main], TRANS_ALPHA)));
    mainStyler.onTooltipShow.subscribe(({node, e}) => {
      if (node) {
        const tooltip = ui.divV([
          ui.div(`${node.name}`)]);
        ui.tooltip.show(tooltip, e.clientX + 16, e.clientY + 16);
      } else {
        ui.tooltip.hide();
      }
    });

    const lightStyler = new TreeStylerBase<MarkupNodeType>('light',
      LINE_WIDTH, NODE_SIZE, false,
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Light], TRANS_ALPHA)),
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Light], TRANS_ALPHA)));

    const currentStyler = new TreeStylerBase<MarkupNodeType>('current',
      LINE_WIDTH, NODE_SIZE, false,
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Current], TRANS_ALPHA)),
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Current], TRANS_ALPHA)));

    const mouseOverStyler = new TreeStylerBase<MarkupNodeType>('mouseOver',
      LINE_WIDTH, NODE_SIZE, false,
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.MouseOver], TRANS_ALPHA)),
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.MouseOver], TRANS_ALPHA)));

    const selectionStyler = new TreeStylerBase<MarkupNodeType>('selection',
      LINE_WIDTH, NODE_SIZE, false,
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Selection], TRANS_ALPHA)),
      toRgba(setAlpha(TreeDefaultPalette[TreeColorNames.Selection], TRANS_ALPHA)));

    // TODO: Attach Dendrogram properties to grid (type or object?)

    return new LeafRangeGridTreeRenderer(grid, tree as MarkupNodeType, placer,
      mainStyler, lightStyler, currentStyler, mouseOverStyler, selectionStyler);
  }
}
