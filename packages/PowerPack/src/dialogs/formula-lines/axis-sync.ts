import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {Editor} from './editor';
import {getAreaColumn, isAxisToken, replaceLeadingCol, setAreaColumn, wrapCol} from './formulas';
import {Preview} from './preview';
import {Table} from './table';
import {AxisColumns, isAnnotationRegionType, ITEM_ORIENTATION, ITEM_TYPE} from './types';

const AXIS_PROPERTIES = ['xColumnName', 'yColumnName', 'yColumnNames', 'valueColumnName'];
const Y_AXIS_PROPERTIES = ['yColumnName', 'yColumnNames'];

/** True if both columns share a non-empty `group` tag (GROK-19736). */
const sameGroup = (a: DG.Column | null, b: DG.Column | null): boolean => {
  const g = a?.getTag('group');
  return !!g && g === b?.getTag('group');
};

/**
 * Makes the current item follow the columns the user picks on the axes of the preview:
 * a band takes the other axis as its adjacent column, a constant line moves to the new column
 * of its axis, and a column of any item gives way to a column of the same group.
 * [getTable] returns the table of the current tab, undefined until the tab is shown.
 */
export function syncItemWithAxes(preview: Preview, editor: Editor, getTable: () => Table | undefined): Subscription {
  return preview.viewer.onPropertyValueChanged.subscribe((typeArgs) => {
    const propertyName = (typeArgs.args as unknown as {property: DG.Property}).property.name;
    const table = getTable();
    const item = table?.currentItem;
    if (editor.applyingChange || !table || !item || !AXIS_PROPERTIES.includes(propertyName))
      return;

    // Snapshot once: the getter re-reads viewer props (and scans columns for line charts) on every access.
    const axisCols = preview.axisCols;
    const df = preview.dataFrame;
    const rowIdx = table.currentItemIdx;
    // Grouped siblings render identically via d4 substitution, so the preview already shows the new axis.
    // Only the editor is refreshed (a preview.update would re-derive the axes and snap
    // 1D value-axis viewers back).
    const refreshEditor = () => editor.applyChange(() => editor.update(rowIdx));

    if (isAnnotationRegionType(item.type)) {
      if (syncRegion(item, axisCols, df))
        refreshEditor();
    } else if (item.type === ITEM_TYPE.BAND) {
      if (syncBandColumn2(item, propertyName, axisCols))
        editor.update(rowIdx);
      if (syncLeadingColumn(item, axisCols, df))
        refreshEditor();
    } else if (item.type === ITEM_TYPE.LINE) {
      const colName = syncConstantLine(item, propertyName, axisCols);
      if (colName) {
        editor.applyChange(() => {
          const col = df.col(colName);
          if (editor.columnInput && col)
            editor.columnInput.value = col;
          preview.update(rowIdx);
        });
      }
      if (syncLeadingColumn(item, axisCols, df))
        refreshEditor();
    }
  });
}

/** The adjacent column of a band follows the axis the band does not cross. Returns true if it has changed. */
function syncBandColumn2(item: DG.FormulaLine, propertyName: string, axisCols: AxisColumns): boolean {
  const isHorz = item.orientation === ITEM_ORIENTATION.HORIZONTAL;
  let newCol2: string | undefined;
  if (isHorz && propertyName === 'xColumnName' && axisCols.x)
    newCol2 = axisCols.x.name;
  else if (!isHorz && Y_AXIS_PROPERTIES.includes(propertyName) && axisCols.y)
    newCol2 = axisCols.y.name;
  // Only an actual change counts: the editor is rebuilt on it, and our own preview.update (triggered
  // by editing title/description/etc.) re-fires onPropertyValueChanged for the same value, so the
  // rebuild would steal focus on every keystroke.
  if (newCol2 === undefined || item.column2 === newCol2)
    return false;
  item.column2 = newCol2;
  return true;
}

/** Constant line follows its own axis directly (keeps the value on the newly selected column).
 *  An axis token already does that on its own, so rebinding it would only pin it to one column.
 *  Returns the column the line has moved to, if any. */
function syncConstantLine(item: DG.FormulaLine, propertyName: string, axisCols: AxisColumns): string | undefined {
  const meta = tryGetMeta(item);
  if (!meta || meta.argName || isAxisToken(meta.funcName))
    return undefined;

  let newCol: string | undefined;
  if (item.orientation === ITEM_ORIENTATION.HORIZONTAL) {
    if (Y_AXIS_PROPERTIES.includes(propertyName))
      newCol = axisCols.y?.name;
  } else if (propertyName === 'xColumnName')
    newCol = axisCols.x?.name;

  if (!newCol || newCol === meta.funcName)
    return undefined;
  item.formula = `${wrapCol(newCol)} = ${meta.expression ?? '0'}`;
  return newCol;
}

/** The parsed formula of a line or band; undefined while it does not parse (an edit in progress). */
function tryGetMeta(item: DG.FormulaLine): DG.FormulaLineMeta | undefined {
  try {
    return DG.FormulaLinesHelper.getMeta(item);
  } catch {
    return undefined;
  }
}

/** Line/band: the leading ${column} follows its same-group sibling onto the axis.
 *  Returns true if the formula has changed. */
function syncLeadingColumn(item: DG.FormulaLine, axisCols: AxisColumns, df: DG.DataFrame): boolean {
  const lineMeta = tryGetMeta(item);
  if (!lineMeta)
    return false;
  // Band → range column (argName). Line → left side (funcName), but only for two-variable
  // lines: a constant line has no argName and is already axis-followed by syncConstantLine.
  const leadingCol = item.type === ITEM_TYPE.BAND ?
    lineMeta.argName :
    (lineMeta.argName ? lineMeta.funcName : undefined);
  const candCol = siblingAxisCol(leadingCol, axisCols, df);
  if (!candCol)
    return false;
  item.formula = replaceLeadingCol(item.formula!, candCol.name);
  return true;
}

/** Annotation regions: area x/y fields and formula1/2 leading columns follow same-group siblings.
 *  Returns true if the region has changed. */
function syncRegion(item: DG.AnnotationRegion, axisCols: AxisColumns, df: DG.DataFrame): boolean {
  let changed = false;

  if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
    const region = item as DG.AreaAnnotationRegion;
    for (const axis of ['x', 'y'] as const) {
      const candCol = siblingAxisCol(getAreaColumn(region, axis), axisCols, df);
      if (candCol) {
        setAreaColumn(region, axis, candCol.name);
        changed = true;
      }
    }
  } else if (item.type === ITEM_TYPE.FORMULA_REGION_ANNOTATION) {
    const region = item as DG.FormulaAnnotationRegion;
    for (const key of ['formula1', 'formula2'] as const) {
      const f = region[key];
      if (!f)
        continue;
      let candCol: DG.Column | null = null;
      try {
        candCol = siblingAxisCol(DG.FormulaLinesHelper.getMetaByFormula(f, ITEM_TYPE.LINE).funcName, axisCols, df);
      } catch {
        continue;
      }
      if (candCol) {
        region[key] = replaceLeadingCol(f, candCol.name);
        changed = true;
      }
    }
  }

  return changed;
}

/** The column on an axis of the preview that is in the same group as the [colName] column of [df]. */
function siblingAxisCol(colName: string | undefined, axisCols: AxisColumns, df: DG.DataFrame): DG.Column | null {
  if (!colName)
    return null;
  const target = df.col(colName);
  for (const c of [axisCols.x, axisCols.y]) {
    if (c && c.name !== colName && sameGroup(c, target))
      return c;
  }
  return null;
}
