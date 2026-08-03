/** Context-aware suggestion engine — powers the toolbox Suggestions pane.
 *
 *  Reads EVERYTHING the canvas knows right now — the selection, the focused
 *  (last-clicked) node, the *data* captured inside nodes (live dataframes →
 *  columns → semantic types), the cell last clicked in the output preview,
 *  and what is already on the canvas — and ranks the ≤10 next steps a user
 *  most likely wants.
 *
 *  Two halves:
 *  - `collectSuggestContext` (async) — reads the live editor + execution
 *    controller into a plain {@link SuggestContext};
 *  - `computeSuggestions` (pure) — rule pipeline over the context and the
 *    function catalog. Testable with hand-built contexts.
 *
 *  The rules are DATA-DRIVEN from the catalog signatures, not per-package
 *  lists: a function declaring `column {semType: Molecule}` is a Molecule
 *  suggestion wherever a Molecule column is detected; `string {semType:
 *  Molecule}` functions fire on a clicked molecule cell; functions with two
 *  dataframe inputs fire on a two-table selection (and get auto-wired).
 *  See docs/func-catalog-snapshot.md + the audit behind this design. */

import * as DG from 'datagrok-api/dg';

import {FlowNode, isExecKey} from '../rete/scheme';
import {FlowEditor} from '../rete/flow-editor';
import {ExecutionController} from '../execution/execution-controller';
import {FuncInfo, getRegisteredFuncs, isCommonNextFunc} from '../rete/node-factory';
import {TypedSocket} from '../rete/sockets';
import {domainSection} from '../types/type-map';
import {getPackageName} from '../utils/dart-proxy-utils';
import {funcWrapperOf} from '../utils/func-input-overrides';

// ---------- context model ----------

export interface ColumnSignal {
  name: string;
  type: string;
  semType: string | null;
}

/** A dataframe-carrying output of a node on the canvas — the raw material of
 *  most suggestions. `columns` is empty until the value is captured by a run. */
export interface TableSignal {
  nodeId: string;
  nodeLabel: string;
  outputKey: string;
  /** True when the output is a `__pt` pass-through (used only when the node
   *  has no real dataframe output). */
  passthrough: boolean;
  /** From the currently selected node set (vs a canvas-wide fallback scan). */
  selected: boolean;
  columns: ColumnSignal[];
}

/** The cell last clicked in the output preview grid. */
export interface CellSignal {
  semType: string | null;
  column: string;
  value: unknown;
}

export interface ScalarSignal {
  nodeId: string;
  outputKey: string;
  dgType: string;
}

export interface SuggestContext {
  nodeCount: number;
  selectedCount: number;
  /** Focus-first: tables of the focused node, then other selected nodes, then
   *  (when nothing is selected) every canvas node with a dataframe output. */
  tables: TableSignal[];
  scalars: ScalarSignal[];
  cell: CellSignal | null;
  /** Lower-cased simple func names already on the canvas. */
  canvasFuncNames: Set<string>;
  /** Domain sections present on the canvas ('Cheminformatics' / 'Bioinformatics'). */
  canvasDomains: Set<string>;
  /** `${nodeId}|${typeName}` pairs that are ALREADY wired — a suggestion whose
   *  every source already feeds a node of that type is dropped (done that). */
  wiredTargets: Set<string>;
}

// ---------- suggestion model ----------

export interface SuggestionWire {
  fromNodeId: string;
  fromOutputKey: string;
  /** Input name on the new node; resolved to the Nth dataframe input at apply
   *  time when omitted. */
  toInput?: string;
}

export interface Suggestion {
  typeName: string;
  label: string;
  /** Short human explanation shown under the label ("Molecule column 'smiles'"). */
  reason: string;
  score: number;
  wire: SuggestionWire[];
  /** `inputValues` to prefill on the created node (column names, a clicked
   *  molecule string, …). Reported via `notifyNodeParamsChanged` on apply. */
  prefill?: Record<string, unknown>;
}

// ---------- catalog signature index ----------

interface FuncSlots {
  info: FuncInfo;
  /** Dataframe input names, declaration order. */
  dfInputs: string[];
  /** Column / column_list inputs with their semType qualifier. */
  colInputs: Array<{name: string; semType: string | null; isList: boolean}>;
  /** Scalar (non-column) inputs carrying a semType (e.g. `string {semType: Molecule}`). */
  semScalarInputs: Array<{name: string; type: string; semType: string}>;
  inputCount: number;
}

const _slotsCache = new Map<string, FuncSlots>();

function slotsFor(info: FuncInfo): FuncSlots {
  let cached = _slotsCache.get(info.nodeTypeName);
  if (cached) return cached;
  const slots: FuncSlots = {info, dfInputs: [], colInputs: [], semScalarInputs: [], inputCount: 0};
  try {
    // A wrapped func (FUNC_WRAPPERS) is matched and auto-wired by what its
    // NODE exposes, not its raw signature — AppendTables reads as a two-table
    // combiner even though the function itself takes one dataframe_list.
    const wrapper = funcWrapperOf(info.func);
    const params: Array<{name: string; type: string; semType: string | null}> = wrapper ?
      wrapper.inputs.map((s) => ({name: s.name, type: String(s.type), semType: null})) :
      info.func.inputs.map((p) => {
        let sem: string | null = null;
        try {
          sem = p.semType ? String(p.semType) : null;
        } catch {/* Dart proxy */}
        return {name: p.name, type: String(p.propertyType), semType: sem};
      });
    for (const p of params) {
      slots.inputCount++;
      if (p.type === 'dataframe')
        slots.dfInputs.push(p.name);
      else if (p.type === 'column' || p.type === 'column_list')
        slots.colInputs.push({name: p.name, semType: p.semType, isList: p.type === 'column_list'});
      else if (p.semType)
        slots.semScalarInputs.push({name: p.name, type: p.type, semType: p.semType});
    }
  } catch {/* introspection failure — empty slots */}
  cached = slots;
  _slotsCache.set(info.nodeTypeName, cached);
  return cached;
}

/** Featured "everyone reaches for these" operations, floated above their
 *  semType siblings (simple func name, lower-cased). Editable. */
export const FEATURED_FUNCS = new Set([
  'descriptors', 'addchempropertiescolumns', 'addchemriskscolumns', 'getinchis', 'curate',
  'chemicalspaceusingumap', 'runelementalanalysis',
  'toatomiclevel', 'pepseamsa', 'splittomonomerstopmenu', 'getregiontopmenu',
]);

function simpleName(info: FuncInfo): string {
  try {
    return String(info.func.name || '').toLowerCase();
  } catch {
    return '';
  }
}

// ---------- scores (exported for tests / tuning) ----------

export const SCORE = {
  semTypeColumn: 100,
  /** A clicked preview cell is a more direct signal than the column it sits
   *  in — but still below {@link SCORE.twoTables} even fully boosted (105 + 6 + 2). */
  cellValue: 105,
  /** Selecting two tables is the most explicit intent signal there is — the
   *  combiners (Join, Append, Compare, …) must outrank every per-table match,
   *  including a semType hit with all its bonuses (100 + 6 + 2 + 3 = 111). */
  twoTables: 115,
  commonNext: 60,
  viewer: 55,
  tableOutput: 42,
  valueOutput: 35,
  emptyCanvas: 50,
  featuredBonus: 6,
  domainBonus: 2,
  selectedBonus: 3,
  /** Per-table decay: suggestions for the 2nd/3rd… context table rank lower. */
  tableDecay: 4,
} as const;

/** Cap of same-rule items per context table, so one Molecule column doesn't
 *  fill all ten slots. */
const PER_TABLE_SEMTYPE_CAP = 6;
export const MAX_SUGGESTIONS = 10;

/** The single-table next steps offered, most-likely first; only the first
 *  {@link COMMON_NEXT_CAP} found in the catalog make it in. */
const COMMON_NEXT_ORDER = [
  'addnewcolumn', 'aggregate', 'filterrows', 'jointables', 'extractcolumns', 'pivot', 'unpivot',
];
const COMMON_NEXT_CAP = 4;

// ---------- the engine ----------

export function computeSuggestions(
  ctx: SuggestContext,
  catalog: FuncInfo[] = getRegisteredFuncs(),
  limit: number = MAX_SUGGESTIONS,
): Suggestion[] {
  const out = new Map<string, Suggestion>();
  const add = (s: Suggestion): void => {
    // A suggestion whose every wire source already feeds a node of this type
    // was already taken — don't nag.
    if (s.wire.length > 0 && s.wire.every((w) => ctx.wiredTargets.has(`${w.fromNodeId}|${s.typeName}`)))
      return;
    const prev = out.get(s.typeName);
    if (!prev || prev.score < s.score) out.set(s.typeName, s);
  };

  const domainBoost = (info: FuncInfo): number => {
    const d = domainSection(getPackageName(info.func));
    return d && ctx.canvasDomains.has(d) ? SCORE.domainBonus : 0;
  };

  if (ctx.nodeCount === 0) {
    add({typeName: 'Inputs/Table Input', label: 'Table Input', score: SCORE.emptyCanvas,
      reason: 'Start with a table', wire: []});
    const openFile = catalog.find((f) => simpleName(f) === 'openfile');
    if (openFile) {
      add({typeName: openFile.nodeTypeName, label: openFile.name, score: SCORE.emptyCanvas + 1,
        reason: 'Start by opening a file', wire: []});
    }
    return rank(out, limit);
  }

  // R1 — semantic-type column rules (Molecule → Chem ops, Macromolecule → Bio
  // ops, PDB_ID → structure fetch, …), fully data-driven from the catalog.
  ctx.tables.forEach((table, ti) => {
    const semCols = table.columns.filter((c) => c.semType);
    if (semCols.length === 0) return;
    const candidates: Array<{slots: FuncSlots; col: ColumnSignal; colInput: string}> = [];
    for (const info of catalog) {
      const slots = slotsFor(info);
      if (slots.dfInputs.length === 0 || slots.colInputs.length === 0) continue;
      for (const ci of slots.colInputs) {
        if (!ci.semType) continue;
        const col = semCols.find((c) => c.semType === ci.semType);
        if (col) {
          candidates.push({slots, col, colInput: ci.name});
          break;
        }
      }
    }
    candidates.sort((a, b) =>
      Number(FEATURED_FUNCS.has(simpleName(b.slots.info))) - Number(FEATURED_FUNCS.has(simpleName(a.slots.info))) ||
      a.slots.info.name.localeCompare(b.slots.info.name));
    for (const {slots, col, colInput} of candidates.slice(0, PER_TABLE_SEMTYPE_CAP)) {
      const featured = FEATURED_FUNCS.has(simpleName(slots.info)) ? SCORE.featuredBonus : 0;
      add({
        typeName: slots.info.nodeTypeName,
        label: slots.info.name,
        reason: `${col.semType} column "${col.name}"`,
        score: SCORE.semTypeColumn - ti * SCORE.tableDecay + featured +
          domainBoost(slots.info) + (table.selected ? SCORE.selectedBonus : 0),
        wire: [{fromNodeId: table.nodeId, fromOutputKey: table.outputKey, toInput: slots.dfInputs[0]}],
        prefill: {[colInput]: col.name},
      });
    }
  });

  // R2 — the cell clicked in the output preview: a Molecule/CHEMBL_ID/… value
  // feeds functions declaring a semType-qualified scalar input, prefilled.
  if (ctx.cell?.semType) {
    for (const info of catalog) {
      const slots = slotsFor(info);
      const si = slots.semScalarInputs.find((s) => s.semType === ctx.cell!.semType);
      if (!si) continue;
      add({
        typeName: info.nodeTypeName,
        label: info.name,
        reason: `Clicked ${ctx.cell.semType} value`,
        score: SCORE.cellValue + domainBoost(info) +
          (FEATURED_FUNCS.has(simpleName(info)) ? SCORE.featuredBonus : 0),
        wire: [],
        prefill: {[si.name]: String(ctx.cell.value ?? '')},
      });
    }
  }

  // R3 — two (or more) selected tables → multi-dataframe functions, auto-wired.
  const selectedTables = distinctByNode(ctx.tables.filter((t) => t.selected));
  if (selectedTables.length >= 2) {
    for (const info of catalog) {
      const slots = slotsFor(info);
      if (slots.dfInputs.length < 2) continue;
      const wires = slots.dfInputs.slice(0, selectedTables.length).map((input, i) => ({
        fromNodeId: selectedTables[i].nodeId,
        fromOutputKey: selectedTables[i].outputKey,
        toInput: input,
      }));
      add({
        typeName: info.nodeTypeName,
        label: info.name,
        reason: `${Math.min(slots.dfInputs.length, selectedTables.length)} selected tables`,
        score: SCORE.twoTables + (isCommonNextFunc(info.nodeTypeName) ? SCORE.featuredBonus : 0) +
          domainBoost(info),
        wire: wires,
      });
    }
  }

  // R4 — a table in context → the common next pipeline steps + viewers.
  const table = ctx.tables[0];
  if (table) {
    const wire = [{fromNodeId: table.nodeId, fromOutputKey: table.outputKey}];
    // The full COMMON_NEXT set would swamp the 10-slot cap with equal-scored
    // plumbing — take the top few in likelihood order and leave room for
    // viewers, outputs, and the semType rules.
    let commonBudget = COMMON_NEXT_CAP;
    for (const name of COMMON_NEXT_ORDER) {
      if (commonBudget === 0) break;
      const info = catalog.find((f) => simpleName(f) === name);
      if (!info) continue;
      const slots = slotsFor(info);
      if (slots.dfInputs.length === 0) continue;
      commonBudget--;
      add({
        typeName: info.nodeTypeName,
        label: info.name,
        reason: `Table from "${table.nodeLabel}"`,
        score: SCORE.commonNext - (COMMON_NEXT_CAP - commonBudget) + (table.selected ? SCORE.selectedBonus : 0),
        wire: [{...wire[0], toInput: slots.dfInputs[0]}],
      });
    }

    // Viewers chosen by what the columns can show. Without captured columns,
    // Grid is always a sane default.
    const numeric = table.columns.filter((c) => c.type === 'int' || c.type === 'double' || c.type === 'float').length;
    const categorical = table.columns.filter((c) => c.type === 'string' && !c.semType).length;
    const viewers: Array<[string, number, string]> = [];
    if (numeric >= 2) viewers.push(['Viewers/Scatter Plot', SCORE.viewer + 3, `${numeric} numeric columns`]);
    if (numeric >= 1) viewers.push(['Viewers/Histogram', SCORE.viewer, 'Numeric column']);
    if (categorical >= 1) viewers.push(['Viewers/Bar Chart', SCORE.viewer - 1, `Category column`]);
    viewers.push(['Viewers/Grid', table.columns.length ? SCORE.viewer - 5 : SCORE.viewer, `Table from "${table.nodeLabel}"`]);
    for (const [typeName, score, reason] of viewers) {
      add({typeName, label: typeName.split('/')[1], reason,
        score: score + (table.selected ? SCORE.selectedBonus : 0),
        wire: [{...wire[0], toInput: 'table'}]});
    }

    add({typeName: 'Outputs/Table Output', label: 'Table Output', reason: 'Publish this table as a flow output',
      score: SCORE.tableOutput, wire: [{...wire[0], toInput: 'table'}]});
  }

  // R5 — scalar outputs on the focus node → Value Output.
  if (!table && ctx.scalars.length > 0) {
    const s = ctx.scalars[0];
    add({typeName: 'Outputs/Value Output', label: 'Value Output',
      reason: `Publish the ${s.dgType} result`, score: SCORE.valueOutput,
      wire: [{fromNodeId: s.nodeId, fromOutputKey: s.outputKey, toInput: 'value'}]});
  }

  return rank(out, limit);
}

function distinctByNode(tables: TableSignal[]): TableSignal[] {
  const seen = new Set<string>();
  const out: TableSignal[] = [];
  for (const t of tables) {
    if (seen.has(t.nodeId)) continue;
    seen.add(t.nodeId);
    out.push(t);
  }
  return out;
}

function rank(out: Map<string, Suggestion>, limit: number): Suggestion[] {
  return Array.from(out.values())
    .sort((a, b) => b.score - a.score || a.label.localeCompare(b.label))
    .slice(0, limit);
}

// ---------- context collection (live editor + execution state) ----------

/** Dataframe-carrying output keys of a node: real dataframe outputs, else (the
 *  user rule) its dataframe pass-throughs. */
export function dataframeOutputKeys(node: FlowNode): {key: string; passthrough: boolean} | null {
  let firstPt: string | null = null;
  for (const [key, out] of Object.entries(node.outputs as Record<string, {socket: TypedSocket} | undefined>)) {
    if (!out || isExecKey(key) || out.socket.dgType !== 'dataframe') continue;
    if (key.endsWith('__pt')) {
      firstPt = firstPt ?? key;
      continue;
    }
    return {key, passthrough: false};
  }
  return firstPt ? {key: firstPt, passthrough: true} : null;
}

function scalarOutputKeys(node: FlowNode): ScalarSignal[] {
  const out: ScalarSignal[] = [];
  for (const [key, o] of Object.entries(node.outputs as Record<string, {socket: TypedSocket} | undefined>)) {
    if (!o || isExecKey(key) || key.endsWith('__pt')) continue;
    const t = o.socket.dgType;
    if (t === 'string' || t === 'int' || t === 'double' || t === 'bool' || t === 'datetime')
      out.push({nodeId: node.id, outputKey: key, dgType: t});
  }
  return out;
}

async function columnsOf(df: DG.DataFrame): Promise<ColumnSignal[]> {
  try {
    // Captured clones may predate semantic-type detection — run it once, so a
    // freshly opened SDF/CSV suggests chem ops without any manual step.
    if (!df.columns.toList().some((c) => c.semType))
      await df.meta.detectSemanticTypes();
    return df.columns.toList().map((c) => ({
      name: c.name, type: String(c.type), semType: c.semType ? String(c.semType) : null,
    }));
  } catch {
    return [];
  }
}

/** Read the live canvas into a plain context. `focusNodeId` — the node the
 *  user last clicked (property-panel node); `cell` — the preview cell last
 *  clicked. Selection is read fresh from the editor. */
export async function collectSuggestContext(
  flow: FlowEditor,
  exec: ExecutionController | null,
  focusNodeId: string | null,
  cell: CellSignal | null,
): Promise<SuggestContext> {
  const nodes = flow.getNodes();
  const selectedIds = new Set(flow.getSelectedNodeIds());
  if (focusNodeId && !nodes.some((n) => n.id === focusNodeId)) focusNodeId = null;

  // Focus-first candidate order: focus node, other selected, then — only when
  // nothing is selected — every other node (canvas-wide fallback scan).
  const ordered: FlowNode[] = [];
  const push = (n: FlowNode): void => {
    if (!ordered.includes(n)) ordered.push(n);
  };
  const focus = focusNodeId ? nodes.find((n) => n.id === focusNodeId) : undefined;
  if (focus && (selectedIds.size === 0 || selectedIds.has(focus.id))) push(focus);
  for (const n of nodes) if (selectedIds.has(n.id)) push(n);
  if (ordered.length === 0) for (const n of nodes) push(n);

  const tables: TableSignal[] = [];
  const scalars: ScalarSignal[] = [];
  for (const node of ordered) {
    const df = dataframeOutputKeys(node);
    if (df) {
      const live = exec?.liveValue(node.id, df.key);
      tables.push({
        nodeId: node.id,
        nodeLabel: node.label,
        outputKey: df.key,
        passthrough: df.passthrough,
        selected: selectedIds.has(node.id),
        columns: live instanceof DG.DataFrame ? await columnsOf(live) : [],
      });
    }
    if (scalars.length === 0 && (selectedIds.has(node.id) || node.id === focusNodeId))
      scalars.push(...scalarOutputKeys(node));
  }

  const canvasFuncNames = new Set<string>();
  const canvasDomains = new Set<string>();
  for (const n of nodes) {
    try {
      const fn = n.dgFunc?.name;
      if (fn) canvasFuncNames.add(String(fn).toLowerCase());
      const pkg = n.dgFunc ? getPackageName(n.dgFunc) : '';
      const d = pkg ? domainSection(pkg) : null;
      if (d) canvasDomains.add(d);
    } catch {/* Dart proxy */}
  }

  const wiredTargets = new Set<string>();
  for (const c of flow.getConnections()) {
    const target = flow.getNodeById(c.target);
    if (target?.dgTypeName) wiredTargets.add(`${c.source}|${target.dgTypeName}`);
  }

  return {
    nodeCount: nodes.length,
    selectedCount: selectedIds.size,
    tables,
    scalars,
    cell,
    canvasFuncNames,
    canvasDomains,
    wiredTargets,
  };
}

