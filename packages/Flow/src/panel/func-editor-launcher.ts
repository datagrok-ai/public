/** Launches a function's own custom editor (declared via `editor:` meta, or on
 *  the explicit allowlist in func-editor-utils.ts) for a func node, seeded with
 *  its real upstream tables.
 *
 *  Mirrors the column picker's table-resolution ladder, but for *every*
 *  dataframe input of the node at once:
 *   - a table input is not connected      → balloon: connect a table first;
 *   - connected and the upstream has run  → reuse its captured output table;
 *   - connected but not yet computed      → offer to run the flow up to that
 *     point, then use the produced table.
 *
 *  With all tables in hand, a `FuncCall` is prepared from the connected values
 *  (live registry) + the panel-edited `inputValues` (column names resolved to
 *  real columns of the seeded table), handed to the platform editor dialog
 *  (`createFuncCallEditor`), and the edited values are written back into
 *  `node.inputValues` — connected inputs are never overridden, and columns come
 *  back as name strings (the panel edits names, not live columns). */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {defaultTableParam} from '../rete/nodes/func-node';
import {ExecutionController} from '../execution/execution-controller';
import {ScriptSettings} from '../compiler/script-emitter';
import {createFuncCallEditor} from '../utils/func-editor-utils';

/** Run semantic-type detection on the resolved upstream tables before they
 *  seed a picker/editor dialog — semtype-filtered column inputs (Molecule, …)
 *  are empty otherwise, because a captured clone may not have been through
 *  detection yet. Guarded per table: a detection failure never blocks the
 *  dialog. Shared by the column picker and the function-editor launcher. */
export async function detectSemanticTypes(tables: Iterable<DG.DataFrame>): Promise<void> {
  for (const t of tables) {
    try {
      await t.meta.detectSemanticTypes();
    } catch {/* detection is best-effort — open the dialog regardless */}
  }
}

/** The dataframe input a column/column_list param resolves against: the node's
 *  explicit `columnTables` association, else the shared numeric-suffix pairing
 *  (`keys2`→`table2`), else the first dataframe input — the same ladder the
 *  compiler uses (`tableExprForColumnParam`). */
export function tableParamForColumn(node: FlowNode, paramName: string, dataframeParams: string[]): string | undefined {
  if (dataframeParams.length === 0) return undefined;
  const assoc = node.properties['columnTables'] as Record<string, string> | undefined;
  const explicit = assoc?.[paramName];
  if (explicit && dataframeParams.includes(explicit)) return explicit;
  return defaultTableParam(paramName, dataframeParams);
}

/** Convert one edited FuncCall input back into its panel representation:
 *  columns → their name, column/string lists → a comma-separated name string
 *  (the panel edits names), plain `list` values stay arrays, primitives pass
 *  through. Returns `undefined` for values the panel can't hold (dataframes). */
export function editorValueToPanelValue(v: unknown, propertyType: string): unknown {
  if (v instanceof DG.DataFrame) return undefined;
  if (v instanceof DG.Column) return v.name;
  if (Array.isArray(v)) {
    const isColumns = v.some((x) => x instanceof DG.Column);
    const names = v.map((x) => (x instanceof DG.Column ? x.name : String(x)));
    if (propertyType === 'column_list' || propertyType === 'string_list' || isColumns)
      return names.join(', ');
    return v; // a plain `list` input holds a JS array in the panel
  }
  return v;
}

/** Write the editor-configured inputs back into the node's panel-editable
 *  values. Connected inputs are never overridden (the wire wins), dataframe
 *  inputs are skipped (always connection-fed), and only slots the panel edits
 *  (`name in inputValues`) are touched. Returns the updated param names. */
export function applyEditorResult(
  node: FlowNode, func: DG.Func, fc: DG.FuncCall, isConnected: (name: string) => boolean,
): string[] {
  const applied: string[] = [];
  for (const inp of func.inputs) {
    const name = inp.name;
    const pt = String(inp.propertyType);
    if (pt === 'dataframe') continue;
    if (isConnected(name)) continue;
    if (!(name in node.inputValues)) continue;
    let v: unknown;
    try {
      v = fc.inputs[name];
    } catch {
      continue;
    }
    if (v === undefined || v === null) continue;
    const panelValue = editorValueToPanelValue(v, pt);
    if (panelValue === undefined) continue;
    node.inputValues[name] = panelValue;
    applied.push(name);
  }
  return applied;
}

export class FuncEditorLauncher {
  constructor(
    private flow: FlowEditor,
    private exec: ExecutionController,
    private getSettings: () => ScriptSettings,
  ) {}

  /** Open the function's own editor dialog for this node; resolves `true` when
   *  the dialog round-trip completed and values were written back. */
  async open(node: FlowNode): Promise<boolean> {
    const func = node.dgFunc;
    if (!func) return false;
    const dfParams = func.inputs.filter((p) => String(p.propertyType) === 'dataframe').map((p) => p.name);

    // Every table input must be wired — the editor needs real tables to seed
    // its column pickers (same rule as choosing columns).
    const unconnected = dfParams.filter((p) => !this.flow.isInputConnected(node.id, p));
    if (unconnected.length > 0) {
      const list = unconnected.map((p) => `“${p}”`).join(', ');
      grok.shell.info(`Connect a table to the ${list} input${unconnected.length > 1 ? 's' : ''} first, ` +
        'then edit the parameters in the editor.');
      return false;
    }

    // Resolve each table: a captured (completed, non-stale) upstream result is
    // reused; anything not yet computed is produced by running its slice —
    // after one confirm covering all of them.
    const tables = new Map<string, DG.DataFrame>();
    const missing: Array<{param: string; srcId: string; srcLabel: string}> = [];
    for (const p of dfParams) {
      const src = this.flow.getInputSource(node.id, p)!;
      const table = this.exec.cloneForNode(src.node.id);
      if (table) tables.set(p, table);
      else missing.push({param: p, srcId: src.node.id, srcLabel: String(src.node.label ?? '')});
    }
    if (missing.length > 0) {
      const ok = await this.confirmRun([...new Set(missing.map((m) => m.srcLabel))]);
      if (!ok) return false;
      for (const m of missing) {
        const table = await this.exec.produceTableForNode(m.srcId, this.getSettings());
        if (!table) {
          grok.shell.error(`The flow ran but no table was produced for “${m.param}”.`);
          return false;
        }
        tables.set(m.param, table);
      }
    }

    // The editor's column pickers filter by semantic type (Molecule, …) — make
    // sure the seeded tables carry their semtypes before the dialog opens (a
    // captured clone may not have been through detection yet).
    await detectSemanticTypes(tables.values());

    const fc = func.prepare(this.buildParams(node, func, tables, dfParams));
    const edited = await createFuncCallEditor(fc);
    applyEditorResult(node, func, edited, (name) => this.flow.isInputConnected(node.id, name));
    void this.flow.updateNode(node.id);
    return true;
  }

  /** Modal confirm before running slices to materialize the upstream tables. */
  private confirmRun(sourceLabels: string[]): Promise<boolean> {
    return new Promise((resolve) => {
      let decided = false;
      const settle = (v: boolean): void => {
        if (!decided) {
          decided = true;
          resolve(v);
        }
      };
      const names = sourceLabels.map((l) => `“${l}”`).join(', ');
      ui.dialog('Run to load tables')
        .add(ui.divText(`The table${sourceLabels.length > 1 ? 's' : ''} from ${names} ` +
          `hasn’t been computed yet. Run the flow up to that point now so the editor can show real data?`))
        .onOK(() => settle(true))
        .onCancel(() => settle(false))
        .show();
    });
  }

  /** Seed the FuncCall parameters: resolved tables for dataframe inputs, live
   *  captured values for other connected inputs (when a prior run stashed
   *  them), and the panel-edited `inputValues` for the rest — column names
   *  resolved to real columns of their seeded table. */
  private buildParams(
    node: FlowNode, func: DG.Func, tables: Map<string, DG.DataFrame>, dfParams: string[],
  ): Record<string, unknown> {
    const params: Record<string, unknown> = {};
    for (const inp of func.inputs) {
      const name = inp.name;
      const pt = String(inp.propertyType);
      if (pt === 'dataframe') {
        const t = tables.get(name);
        if (t) params[name] = t;
        continue;
      }
      if (this.flow.isInputConnected(node.id, name)) {
        const src = this.flow.getInputSource(node.id, name)!;
        const live = this.exec.liveValue(src.node.id, src.outputKey);
        if (live !== undefined && live !== null) params[name] = live;
        continue;
      }
      if (!(name in node.inputValues)) continue;
      const v = node.inputValues[name];
      if (v === undefined || v === null) continue;
      if (pt === 'column' || pt === 'column_list') {
        const tp = tableParamForColumn(node, name, dfParams);
        const t = tp ? tables.get(tp) : undefined;
        if (!t) continue;
        const names = String(v).split(',').map((s) => s.trim()).filter(Boolean);
        if (pt === 'column') {
          const col = names.length > 0 ? t.col(names[0]) : null;
          if (col) params[name] = col;
        } else {
          const cols = names.map((n) => t.col(n)).filter((c): c is DG.Column => !!c);
          if (cols.length > 0) params[name] = cols;
        }
        continue;
      }
      if (pt === 'string_list') {
        // The panel edits a comma-separated string; the funccall wants an array.
        const items = String(v).split(',').map((s) => s.trim()).filter(Boolean);
        if (items.length > 0) params[name] = items;
        continue;
      }
      if (String(v) === '' && pt !== 'string') continue;
      params[name] = v;
    }
    return params;
  }
}
