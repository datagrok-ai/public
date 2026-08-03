/** Column picker for func-node column / column-list inputs.
 *
 *  A column input on a node (e.g. AddNewColumn's `column`, or JoinTables'
 *  `keys1` / `keys2` / `values1` / `values2`) resolves against one of the
 *  node's dataframe inputs. Typing column names from memory is brittle — this
 *  drops a real column selector *menu* right under the picker icon, seeded by
 *  the *actual* upstream table, so the user chooses from a list.
 *
 *  Three cases, keyed off the table this column refers to:
 *   - the table input is not connected  → tell the user to connect a table;
 *   - it's connected and the upstream node has already run → pick immediately
 *     from its captured output table;
 *   - it's connected but not yet computed → offer to run the flow up to that
 *     point, then pick from the produced table.
 *
 *  The picker used to open a modal dialog after the run confirmation — two
 *  dialogs back to back. Now the run confirmation stays modal, but the column
 *  choice itself is a `DG.Menu` popped next to the icon (`singleColumnSelector`
 *  / `multiColumnSelector`) — one click, no second dialog.
 *
 *  When the input is a DG func param carrying a `semType` and/or
 *  `columnTypeFilter` (`numerical` | `categorical` | `int` | `double` |
 *  `string`), the menu is filtered to matching columns only.
 *
 *  Multi-table funcs (JoinTables) carry a per-column `columnTables`
 *  association choosing which dataframe input each column resolves against;
 *  the request already carries the resolved `tableParam`. */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {ExecutionController} from '../execution/execution-controller';
import {ScriptSettings} from '../compiler/script-emitter';
import {detectSemanticTypes} from './func-editor-launcher';

export interface ColumnPickRequest {
  /** Node whose column input is being edited. */
  nodeId: string;
  /** The column / column-list input name. */
  paramName: string;
  /** True for `column_list` (multi-select), false for a single `column`. */
  isList: boolean;
  /** The dataframe input name this column resolves against. */
  tableParam: string;
  /** Current field value (comma-separated for a list). */
  current: string;
  /** Write the chosen name(s) back into the field. */
  apply: (value: string) => void;
  /** The picker icon — the menu is popped up next to it. */
  anchor?: HTMLElement;
}

/** columnTypeFilter values `Column.matches` understands (anything else → skip). */
const COLUMN_TYPE_FILTERS = ['numerical', 'categorical', 'int', 'double', 'string'];

/** Column-menu filter for a param's `semType` / `columnTypeFilter`, or
 *  `undefined` when neither constrains the choice (menu shows every column).
 *  A recognised `columnTypeFilter` (numerical | categorical | int | double |
 *  string) narrows by type; a non-empty `semType` narrows by semantic type;
 *  both together require both. Exported for unit testing. */
export function buildColumnMatchFilter(
  semType: string | null | undefined,
  columnTypeFilter: string | null | undefined,
): ((c: DG.Column) => boolean) | undefined {
  const tests: Array<(c: DG.Column) => boolean> = [];
  if (typeof semType === 'string' && semType.length > 0)
    tests.push((c) => c.semType === semType);
  if (typeof columnTypeFilter === 'string' && COLUMN_TYPE_FILTERS.includes(columnTypeFilter))
    tests.push((c) => c.matches(columnTypeFilter as DG.ColumnType | 'numerical' | 'categorical'));
  if (tests.length === 0) return undefined;
  return (c) => tests.every((t) => t(c));
}

export class ColumnPicker {
  constructor(
    private flow: FlowEditor,
    private exec: ExecutionController,
    private getSettings: () => ScriptSettings,
  ) {}

  async pick(req: ColumnPickRequest): Promise<void> {
    const src = this.flow.getInputSource(req.nodeId, req.tableParam);
    if (!src) {
      grok.shell.info(`Connect a table to the “${req.tableParam}” input first, then choose columns.`);
      return;
    }

    let table = this.exec.cloneForNode(src.node.id);
    if (!table) {
      const shouldRun = await this.confirmRun(src.node.label);
      if (!shouldRun) return;
      table = await this.exec.produceTableForNode(src.node.id, this.getSettings());
      if (!table) {
        grok.shell.error('The flow ran but no table was produced for that input.');
        return;
      }
    }
    // Semtype-filtered column inputs (Molecule, …) need the table's semantic
    // types detected before the picker opens — a captured clone may not have
    // been through detection yet.
    await detectSemanticTypes([table]);
    this.openMenu(table, req);
  }

  /** Build a column filter from the func param's `semType` / `columnTypeFilter`.
   *  Non-func nodes (Select Column utilities) and params without either
   *  attribute get no filter — the menu shows every column. */
  private buildColumnFilter(req: ColumnPickRequest): ((c: DG.Column) => boolean) | undefined {
    const prop = this.flow.getNodeById(req.nodeId)?.dgFunc?.inputs
      .find((p) => p.name === req.paramName);
    if (!prop) return undefined;
    return buildColumnMatchFilter(prop.semType, prop.columnTypeFilter);
  }

  /** Modal confirm before running a slice to materialize the upstream table. */
  private confirmRun(sourceLabel: string): Promise<boolean> {
    return new Promise((resolve) => {
      let decided = false;
      const settle = (v: boolean): void => {if (!decided) {decided = true; resolve(v);}};
      ui.dialog('Run to load columns')
        .add(ui.divText(`The table from “${sourceLabel}” hasn’t been computed yet. ` +
          'Run the flow up to that point now to load its columns?'))
        .onOK(() => settle(true))
        .onCancel(() => settle(false))
        .show();
    });
  }

  /** Pop a column selector menu next to the picker icon. Single inputs use
   *  `singleColumnSelector` (click → set → close); lists use
   *  `multiColumnSelector` writing back the checked set live on each toggle. */
  private openMenu(table: DG.DataFrame, req: ColumnPickRequest): void {
    const names = table.columns.names();
    const columnFilter = this.buildColumnFilter(req);
    const show = (menu: DG.Menu): void => {
      const bb = req.anchor?.getBoundingClientRect();
      menu.show(bb ? {x: bb.x + 30, y: bb.y + 30, element: document.body} : {});
    };

    if (req.isList) {
      const checked = req.current.split(',').map((s) => s.trim()).filter((n) => n && names.includes(n));
      const menu = DG.Menu.popup().multiColumnSelector(table, {
        initialValue: checked,
        columnFilter,
        editable: true,
        onChange: (grid) => req.apply(grid.getCheckedColumnNames().join(', ')),
      });
      menu.closeOnClick = false;
      show(menu);
    } else {
      const cur = req.current.trim();
      show(DG.Menu.popup().singleColumnSelector(table, {
        initialValue: names.includes(cur) ? cur : undefined,
        columnFilter,
        changeOnHover: false,
        closeOnClick: true,
        onChange: (_grid, column) => {if (column) req.apply(column.name);},
      }));
    }
  }
}
