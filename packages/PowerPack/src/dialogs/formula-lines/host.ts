import * as DG from 'datagrok-api/dg';

import {ItemList} from './item-list';
import {ITEM_CAPTION} from './types';
import {axisProps, isSingleAxis} from './viewer-axes';

/**
 * The table or viewer the dialog edits: its formula lines and annotation regions,
 * and what its viewer type allows and defaults to.
 */
export class Host {
  /** Items of the viewer; absent when the dialog is opened for a table. */
  public readonly viewerItems?: ItemList;
  public readonly dframeItems: ItemList;

  /** Captions of the lines and bands the viewer can show; `undefined` means all of them. */
  public readonly allowedItems?: string[];
  /** Default z-index of new bands and regions: above markers for histogram / bar chart / density plot. */
  public readonly defaultItemZIndex: number;
  /** Color of new bands; `undefined` leaves the platform default. */
  public readonly defaultBandColor?: string;
  /** The value column of a single-axis viewer: the only column its formulas can refer to. */
  public readonly singleAxisColumnName?: string;

  constructor(public readonly src: DG.DataFrame | DG.Viewer) {
    if (!src)
      throw new Error('Source table/viewer not found.');

    const viewer = src instanceof DG.Viewer ? src : null;
    const df = viewer ? viewer.dataFrame : src as DG.DataFrame;
    if (!df)
      throw new Error('Viewer not attached to table.');

    if (viewer)
      this.viewerItems = new ItemList(viewer.meta);
    this.dframeItems = new ItemList(df.meta);

    this.allowedItems = viewer ? allowedItems(viewer) : undefined;
    this.defaultItemZIndex = viewer?.type === DG.VIEWER.HISTOGRAM || viewer?.type === DG.VIEWER.BAR_CHART ||
      viewer?.type === DG.VIEWER.DENSITY_PLOT ? 100 : -100;
    this.defaultBandColor = viewer?.type === DG.VIEWER.BAR_CHART ? '#838383' : undefined;
    this.singleAxisColumnName = viewer && isSingleAxis(viewer) ? viewer.getOptions().look.valueColumnName : undefined;
  }

  /** The items the dialog shows first. */
  public get initialItems(): ItemList {return this.viewerItems ?? this.dframeItems;}

  /** Formula lines of the source other than the one of [items]. */
  public otherFormulaLines(items: ItemList): DG.FormulaLine[] {
    return (items === this.viewerItems ? this.dframeItems : this.viewerItems)?.formulaLines ?? [];
  }

  public save(): void {
    this.dframeItems.save();
    this.viewerItems?.save();
  }
}

/** Single-axis viewers show only the lines and bands that cross their value axis. */
function allowedItems(viewer: DG.Viewer): string[] | undefined {
  if (!isSingleAxis(viewer))
    return undefined;
  return axisProps(viewer).y ?
    [ITEM_CAPTION.HORZ_LINE, ITEM_CAPTION.HORZ_BAND] :
    [ITEM_CAPTION.VERT_LINE, ITEM_CAPTION.VERT_BAND];
}
