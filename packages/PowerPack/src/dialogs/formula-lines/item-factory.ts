import * as DG from 'datagrok-api/dg';

import {wrapCol} from './formulas';
import {AxisColumns, BTN_CAPTION, EditorItem, isAnnotationRegionType, ITEM_CAPTION, ITEM_ORIENTATION,
  ITEM_TYPE} from './types';

/** Host-specific defaults of the new items. */
export interface NewItemDefaults {
  /** Hex color override for newly-created bands. `undefined` = use `setDefaults`' value. */
  defaultBandColor?: string,
  /** Default `zIndex` stamped on freshly-created bands and annotation regions. */
  defaultItemZIndex: number,
}

// Two-axis items require both columns; a single-axis host (box / histogram / bar) only exposes one.
const NEEDS_X: string[] = [ITEM_CAPTION.LINE, ITEM_CAPTION.VERT_LINE, ITEM_CAPTION.VERT_BAND,
  ITEM_CAPTION.FORMULA_REGION];
const NEEDS_Y: string[] = [ITEM_CAPTION.LINE, ITEM_CAPTION.HORZ_LINE, ITEM_CAPTION.HORZ_BAND,
  ITEM_CAPTION.FORMULA_REGION];

const hasRequiredColumns = (caption: string, cols: AxisColumns): boolean =>
  !(NEEDS_X.includes(caption) && !cols.x || NEEDS_Y.includes(caption) && !cols.y);

/** Returns Formula Line type by its user-friendly [caption]. */
function getItemTypeByCaption(caption: string): string {
  switch (caption) {
  case ITEM_CAPTION.LINE:
  case ITEM_CAPTION.HORZ_LINE:
  case ITEM_CAPTION.VERT_LINE: return ITEM_TYPE.LINE;

  case ITEM_CAPTION.BAND:
  case ITEM_CAPTION.HORZ_BAND:
  case ITEM_CAPTION.VERT_BAND: return ITEM_TYPE.BAND;

  default: throw new Error('Unknown item caption.');
  }
}

/** Creates the item of the [caption] menu item on the [cols] axes, or a clone of [currentItem];
 *  null when the viewer lacks an axis column the item needs.
 *  [valX] and [valY] are the point of the viewer a constant line has to pass through. */
export function createItem(caption: string, cols: AxisColumns, defaults: NewItemDefaults,
  currentItem: EditorItem | null, valY?: number, valX?: number): EditorItem | null {
  const colY = cols.y;
  const colX = cols.x;
  // Stats override (used by aggregated hosts like bar chart) — fall through
  // to the column's own stats when not provided.
  const yStats = cols.yStats ?? colY?.stats;
  const xStats = cols.xStats ?? colX?.stats;
  if (!hasRequiredColumns(caption, cols))
    return null;

  const region = (formula1: string, formula2: string): DG.FormulaAnnotationRegion =>
    ({type: ITEM_TYPE.FORMULA_REGION_ANNOTATION, formula1, formula2, zIndex: defaults.defaultItemZIndex});

  let item: DG.FormulaLine = {};

  switch (caption) {
  case ITEM_CAPTION.FORMULA_REGION: {
    const line = `${wrapCol(colY!.name)} = ${wrapCol(colX!.name)}`;
    const offset = yStats!.q2.toFixed(1);
    return region(`${line} + ${offset}`, `${line} - ${offset}`);
  }

  // Axis-bounds formula region — `${col} = q1` / `${col} = q3` for the chosen axis.
  // Mirrors `addAxisAnnotationsMenu`'s `Add Region` in `core/client/d4/lib/src/utils/utils.dart`.
  case ITEM_CAPTION.VERT_RANGE_REGION:
  case ITEM_CAPTION.HORZ_RANGE_REGION: {
    const col = caption === ITEM_CAPTION.VERT_RANGE_REGION ? colX : colY;
    const stats = caption === ITEM_CAPTION.VERT_RANGE_REGION ? xStats : yStats;
    if (!col || !stats)
      return null;
    const lhs = wrapCol(col.name);
    return region(`${lhs} = ${stats.q1.toFixed(1)}`, `${lhs} = ${stats.q3.toFixed(1)}`);
  }

  case ITEM_CAPTION.LINE:
    item.formula = `${wrapCol(colY!.name)} = ${wrapCol(colX!.name)}`;
    break;

  case ITEM_CAPTION.VERT_LINE:
    item.orientation = ITEM_ORIENTATION.VERTICAL;
    item.formula = `${wrapCol(colX!.name)} = ${(valX ?? xStats!.q2).toFixed(1)}`;
    break;

  case ITEM_CAPTION.HORZ_LINE:
    item.orientation = ITEM_ORIENTATION.HORIZONTAL;
    item.formula = `${wrapCol(colY!.name)} = ${(valY ?? yStats!.q2).toFixed(1)}`;
    break;

  case ITEM_CAPTION.VERT_BAND:
    item.formula = `${wrapCol(colX!.name)} in(${xStats!.q1.toFixed(1)}, ${xStats!.q3.toFixed(1)})`;
    item.orientation = ITEM_ORIENTATION.VERTICAL;
    // For single-axis hosts (no Y), fall back to the X column for `column2`.
    item.column2 = (colY ?? colX!).name;
    break;

  case ITEM_CAPTION.HORZ_BAND:
    item.formula = `${wrapCol(colY!.name)} in(${yStats!.q1.toFixed(1)}, ${yStats!.q3.toFixed(1)})`;
    item.orientation = ITEM_ORIENTATION.HORIZONTAL;
    // For single-axis hosts (no X), fall back to the Y column for `column2`.
    item.column2 = (colX ?? colY!).name;
    break;

  case BTN_CAPTION.CLONE:
    // A region takes no formula line defaults
    if (isAnnotationRegionType(currentItem!.type))
      return structuredClone(currentItem!);
    item = currentItem!;
    break;
  }

  item.type ??= getItemTypeByCaption(caption);

  item = DG.FormulaLinesHelper.setDefaults(item);

  // Host-specific overrides for new bands (color + above/below markers).
  if (item.type === ITEM_TYPE.BAND) {
    if (defaults.defaultBandColor)
      item.color = defaults.defaultBandColor;
    item.zIndex = defaults.defaultItemZIndex;
  }

  return item;
}
