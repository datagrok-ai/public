/** Column picker for func-node column / column-list inputs — a DG.Menu column selector next to the
 *  picker icon, seeded by the actual upstream table (run up to that point if needed). */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {ExecutionController} from '../execution/execution-controller';
import {ScriptSettings} from '../compiler/script-emitter';
import {detectSemanticTypes} from './func-editor-launcher';

export interface ColumnPickRequest {
  nodeId: string;
  paramName: string;
  isList: boolean;
  /** The dataframe input name this column resolves against. */
  tableParam: string;
  /** Current field value (comma-separated for a list). */
  current: string;
  apply: (value: string) => void;
  anchor?: HTMLElement;
}

/** columnTypeFilter values `Column.matches` understands (anything else → skip). */
const COLUMN_TYPE_FILTERS = ['numerical', 'categorical', 'int', 'double', 'string'];

/** Column filter from a param's `semType` / `columnTypeFilter`; `undefined` when neither constrains. */
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
    await detectSemanticTypes([table]);
    this.openMenu(table, req);
  }

  private buildColumnFilter(req: ColumnPickRequest): ((c: DG.Column) => boolean) | undefined {
    const prop = this.flow.getNodeById(req.nodeId)?.dgFunc?.inputs
      .find((p) => p.name === req.paramName);
    if (!prop) return undefined;
    return buildColumnMatchFilter(prop.semType, prop.columnTypeFilter);
  }

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
