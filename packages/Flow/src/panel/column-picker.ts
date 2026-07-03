/** Column picker for func-node column / column-list inputs.
 *
 *  A column input on a node (e.g. AddNewColumn's `column`, or JoinTables'
 *  `keys1` / `keys2` / `values1` / `values2`) resolves against one of the
 *  node's dataframe inputs. Typing column names from memory is brittle — this
 *  opens a real column / columns picker dialog seeded by the *actual* upstream
 *  table, so the user chooses from a list.
 *
 *  Three cases, keyed off the table this column refers to:
 *   - the table input is not connected  → tell the user to connect a table;
 *   - it's connected and the upstream node has already run → pick immediately
 *     from its captured output table;
 *   - it's connected but not yet computed → offer to run the flow up to that
 *     point, then pick from the produced table.
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
    this.openDialog(table, req);
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

  private openDialog(table: DG.DataFrame, req: ColumnPickRequest): void {
    const names = table.columns.names();
    if (req.isList) {
      const checked = req.current.split(',').map((s) => s.trim()).filter((n) => n && names.includes(n));
      const input = ui.input.columns(req.paramName, {table, value: table.columns.byNames(checked)});
      ui.dialog('Select columns')
        .add(input.root)
        .onOK(() => req.apply((input.value ?? []).map((c) => c.name).join(', ')))
        .show();
    } else {
      const cur = req.current.trim();
      const input = ui.input.column(req.paramName, {
        table, value: names.includes(cur) ? table.col(cur) ?? undefined : undefined,
      });
      ui.dialog('Select column')
        .add(input.root)
        .onOK(() => {if (input.value) req.apply(input.value.name);})
        .show();
    }
  }
}
