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
 */
export class AxisSync {
  constructor(
    private preview: Preview,
    private editor: Editor,
    /** The table of the current tab; undefined until the tab is shown */
    private getTable: () => Table | undefined,
  ) {}

  public subscribe(): Subscription {
    return this.preview.viewer.onPropertyValueChanged.subscribe((typeArgs) => {
      const {property} = typeArgs.args as unknown as {property: DG.Property};
      this.onPropertyChanged(property.name);
    });
  }

  private onPropertyChanged(propertyName: string): void {
    const table = this.getTable();
    const currentItem = table?.currentItem;
    if (this.editor.applyingChange || !table || !currentItem || !AXIS_PROPERTIES.includes(propertyName))
      return;

    // Snapshot once: the getter re-reads viewer props (and scans columns for line charts) on every access.
    const axisCols = this.preview.axisCols;
    const rowIdx = table.currentItemIdx;

    if (isAnnotationRegionType(currentItem.type)) {
      // Region renders on the current axis via d4 substitution; refresh the editor inputs only
      // (a preview.update here would re-derive axes and snap 1D value-axis viewers back).
      if (this.syncRegion(currentItem, axisCols))
        this.editor.applyChange(() => this.editor.update(rowIdx));
    } else if (currentItem.type === ITEM_TYPE.BAND) {
      this.syncBandColumn2(currentItem, propertyName, axisCols, rowIdx);
      this.syncLeadingColumn(currentItem, axisCols, rowIdx);
    } else if (currentItem.type === ITEM_TYPE.LINE) {
      this.syncConstantLine(currentItem, propertyName, axisCols, rowIdx);
      this.syncLeadingColumn(currentItem, axisCols, rowIdx);
    }
  }

  /** The adjacent column of a band follows the axis the band does not cross. */
  private syncBandColumn2(item: DG.FormulaLine, propertyName: string, axisCols: AxisColumns, rowIdx: number): void {
    const isHorz = item.orientation === ITEM_ORIENTATION.HORIZONTAL;
    // Only rebuild the editor when the axis column actually changed; otherwise our own
    // preview.update (triggered by editing title/description/etc.) would steal focus on
    // every keystroke, since setOptions re-fires onPropertyValueChanged for the same value.
    let newCol2: string | undefined;
    if (isHorz && propertyName === 'xColumnName' && axisCols.x)
      newCol2 = axisCols.x.name;
    else if (!isHorz && Y_AXIS_PROPERTIES.includes(propertyName) && axisCols.y)
      newCol2 = axisCols.y.name;
    if (newCol2 !== undefined && item.column2 !== newCol2) {
      item.column2 = newCol2;
      this.editor.update(rowIdx);
    }
  }

  /** Constant line follows its own axis directly (keeps the value on the newly selected column).
   *  An axis token already does that on its own, so rebinding it would only pin it to one column. */
  private syncConstantLine(item: DG.FormulaLine, propertyName: string, axisCols: AxisColumns, rowIdx: number): void {
    const meta = DG.FormulaLinesHelper.getMeta(item);
    if (meta.argName || isAxisToken(meta.funcName))
      return;

    let newCol: string | undefined;
    if (item.orientation === ITEM_ORIENTATION.HORIZONTAL) {
      if (Y_AXIS_PROPERTIES.includes(propertyName))
        newCol = axisCols.y?.name;
    } else if (propertyName === 'xColumnName')
      newCol = axisCols.x?.name;

    if (newCol && newCol !== meta.funcName) {
      const colName = newCol;
      this.editor.applyChange(() => {
        item.formula = `${wrapCol(colName)} = ${meta.expression ?? '0'}`;
        const col = this.preview.dataFrame.col(colName);
        if (this.editor.columnInput && col)
          this.editor.columnInput.value = col;
        this.preview.update(rowIdx);
      });
    }
  }

  /** Line/band: the leading ${column} follows its same-group sibling onto the axis. */
  private syncLeadingColumn(item: DG.FormulaLine, axisCols: AxisColumns, rowIdx: number): void {
    // getMeta parses the formula and throws on in-progress/invalid input — bail rather than crash the handler.
    let lineMeta: DG.FormulaLineMeta;
    try {
      lineMeta = DG.FormulaLinesHelper.getMeta(item);
    } catch {
      return;
    }
    // Band → range column (argName). Line → left side (funcName), but only for two-variable
    // lines: a constant line has no argName and is already axis-followed by syncConstantLine.
    const leadingCol = item.type === ITEM_TYPE.BAND ?
      lineMeta.argName :
      (lineMeta.argName ? lineMeta.funcName : undefined);
    const candCol = this.siblingAxisCol(leadingCol, axisCols);
    if (candCol) {
      // Grouped siblings render identically via d4 substitution, so the preview already
      // reflects the user's axis change — only the editor's displayed formula needs refreshing.
      this.editor.applyChange(() => {
        item.formula = replaceLeadingCol(item.formula!, candCol.name);
        this.editor.update(rowIdx);
      });
    }
  }

  /** Annotation regions: area x/y fields and formula1/2 leading columns follow same-group siblings.
   *  Returns true if the region has changed. */
  private syncRegion(item: DG.AnnotationRegion, axisCols: AxisColumns): boolean {
    let changed = false;

    if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
      const region = item as DG.AreaAnnotationRegion;
      for (const axis of ['x', 'y'] as const) {
        const candCol = this.siblingAxisCol(getAreaColumn(region, axis), axisCols);
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
          candCol = this.siblingAxisCol(DG.FormulaLinesHelper.getMetaByFormula(f, ITEM_TYPE.LINE).funcName, axisCols);
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

  /** The column on an axis of the preview that is in the same group as the [colName] column. */
  private siblingAxisCol(colName: string | undefined, axisCols: AxisColumns): DG.Column | null {
    if (!colName)
      return null;
    const target = this.preview.dataFrame.col(colName);
    for (const c of [axisCols.x, axisCols.y]) {
      if (c && c.name !== colName && sameGroup(c, target))
        return c;
    }
    return null;
  }
}
