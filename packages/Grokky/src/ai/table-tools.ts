import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {undoable} from '../claude/undoable';
import type {ViewToolDef} from './view-tools';

const UNDO = ' Recorded as one undo step (Edit | Undo, Ctrl+Z) when the result has `undoable: true`.';
const ROW_CODE = 'JS expression evaluated for every row: `i` is the row index, `t` the table, `v` the current value ' +
  'where there is one (null when missing).';
const column = {type: 'string', description: 'Column name.'};

export const TABLE_TOOL_DEFS: ViewToolDef[] = [
  {
    name: 'table_remove_columns',
    description: 'Remove columns from the table of the current table view.' + UNDO,
    inputSchema: {type: 'object', required: ['columns'], properties: {
      columns: {type: 'array', items: {type: 'string'}, description: 'Names of the columns to remove.'}}},
  },
  {
    name: 'table_remove_rows',
    description: 'Remove rows from the table of the current table view. The user\'s row selection is kept.' + UNDO,
    inputSchema: {type: 'object', required: ['where'], properties: {
      where: {type: 'string',
        description: `True for a row to remove, e.g. "t.getCol('HEIGHT').isNone(i)". ${ROW_CODE}`}}},
  },
  {
    name: 'table_add_column',
    description: 'Add a column to the table of the current table view, filled from `value` or calculated by ' +
      '`formula`.' + UNDO,
    inputSchema: {type: 'object', required: ['name'], properties: {
      name: {type: 'string', description: 'Name of the new column; made unique when taken.'},
      type: {type: 'string', enum: ['double', 'int', 'string', 'bool', 'datetime'],
        description: 'Type of a `value` column.'},
      value: {type: 'string', description: `Value of the row, null for a missing one. ${ROW_CODE}`},
      formula: {type: 'string', description: 'Platform formula of a calculated column, e.g. "${WEIGHT} / 2".'}}},
  },
  {
    name: 'table_set_values',
    description: 'Overwrite the values of a column of the current table view.' + UNDO,
    inputSchema: {type: 'object', required: ['column', 'value'], properties: {
      column, value: {type: 'string', description: `New value of the row, null for a missing one. ${ROW_CODE}`}}},
  },
  {
    name: 'table_edit_column',
    description: 'Rename a column of the current table view or set its metadata.' + UNDO,
    inputSchema: {type: 'object', required: ['column'], properties: {
      column,
      name: {type: 'string', description: 'New name.'},
      format: {type: 'string', description: 'Display format, e.g. "0.00" or "yyyy-MM-dd"; "" clears it.'},
      units: {type: 'string', description: 'Units; "" clears them.'},
      semType: {type: 'string', description: 'Semantic type, e.g. "Molecule".'},
      colorCoding: {type: 'string', enum: ['Linear', 'Categorical', 'Conditional', 'Off'],
        description: 'Color coding with the platform\'s default scheme.'},
      tags: {type: 'object', description: 'Other column tags to set, by name.'}}},
  },
  {
    name: 'table_add_viewer',
    description: 'Add a viewer to the current table view.' + UNDO,
    inputSchema: {type: 'object', required: ['type'], properties: {
      type: {type: 'string', description: 'Viewer type, e.g. "Scatter plot", "Histogram", "Bar chart".'},
      options: {type: 'object',
        description: 'Viewer properties, e.g. {"xColumnName": "HEIGHT", "yColumnName": "WEIGHT"}.'}}},
  },
];

const COLOR_CODING: {[type: string]: (col: DG.Column) => void} = {
  Linear: (col) => col.meta.colors.setLinear(),
  Categorical: (col) => col.meta.colors.setCategorical(),
  Conditional: (col) => col.meta.colors.setConditional(),
  Off: (col) => col.meta.colors.setDisabled(),
};

function tableView(): DG.TableView {
  const view = grok.shell.v;
  if (view?.type !== DG.VIEW_TYPE.TABLE_VIEW)
    throw new Error('The current view is not a table view');
  return view as DG.TableView;
}

function rowCode(code: string, t: DG.DataFrame, col?: DG.Column): (i: number) => any {
  const f = new Function('i', 't', 'v', 'grok', 'DG', `return (${code}\n);`);
  return (i) => f(i, t, col ? (col.isNone(i) ? null : col.get(i)) : undefined, grok, DG);
}

type Runner = (args: any) => any;

export const TABLE_TOOL_RUNNERS = new Map<string, Runner>([
  ['table_remove_columns', (args) => {
    const t = tableView().dataFrame;
    const recorded = undoable.removeColumns(t, args.columns);
    return {success: true, undoable: recorded, columns: t.columns.names()};
  }],
  ['table_remove_rows', async (args) => {
    const t = tableView().dataFrame;
    const before = t.rowCount;
    const recorded = await undoable.removeRows(t, rowCode(args.where, t));
    return {success: true, undoable: recorded, removed: before - t.rowCount, rowCount: t.rowCount};
  }],
  ['table_add_column', async (args) => {
    const t = tableView().dataFrame;
    const name = t.columns.getUnusedName(args.name);
    const value = args.formula ? null : rowCode(args.value, t);
    const col = value ? t.columns.addNew(name, args.type ?? DG.TYPE.FLOAT) :
      await t.columns.addNewCalculated(name, args.formula);
    try {
      if (value)
        col.init(value);
    } catch (e) {
      t.columns.remove(col.name);
      throw e;
    }
    return {success: true, undoable: undoable.columnAdded(col), name: col.name, type: col.type,
      missing: col.stats.missingValueCount};
  }],
  ['table_set_values', async (args) => {
    const t = tableView().dataFrame;
    const col = t.getCol(args.column);
    const value = rowCode(args.value, t, col);
    const recorded = await undoable.changeValues(col, () => col.init(value));
    return {success: true, undoable: recorded, column: col.name, missing: col.stats.missingValueCount};
  }],
  ['table_edit_column', async (args) => {
    const t = tableView().dataFrame;
    const col = t.getCol(args.column);
    if (args.name != null && t.col(args.name) && t.col(args.name) !== col)
      throw new Error(`The table already has a column named "${args.name}"`);
    const recorded = await undoable.editColumn(col, () => {
      if (args.name != null)
        col.name = args.name;
      if (args.format != null)
        col.meta.format = args.format;
      if (args.units != null)
        col.meta.units = args.units;
      if (args.semType != null)
        col.semType = args.semType;
      COLOR_CODING[args.colorCoding]?.(col);
      Object.assign(col.tags, args.tags ?? {});
    });
    return {success: true, undoable: recorded, name: col.name, colorCoding: col.meta.colors.getType(),
      format: col.meta.format, units: col.meta.units, semType: col.semType};
  }],
  ['table_add_viewer', (args) => {
    const viewer = undoable.addViewer(tableView(), args.type, args.options);
    return {success: true, undoable: DG.UndoService != null, type: viewer.type};
  }],
]);
