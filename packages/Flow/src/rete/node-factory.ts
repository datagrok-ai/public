/** Node-type registry: maps type names like `Inputs/Table Input` to zero-arg FlowNode factories. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ClassicPreset} from 'rete';
import {FlowNode, EXEC_IN_KEY, EXEC_OUT_KEY, ORDER_SOCKET_TYPE} from './scheme';
import {FuncNode} from './nodes/func-node';
import {TypedSocket, getSocket} from './sockets';
import {areTypesCompatible, categorizeBySignature, domainCategory, domainSection} from '../types/type-map';

import {
  TableInputNode, ColumnInputNode, ColumnListInputNode, StringInputNode, MoleculeInputNode,
  HelmInputNode, NumberInputNode, IntInputNode, BooleanInputNode, DateTimeInputNode,
  FileInputNode, MapInputNode, DynamicInputNode, StringListInputNode, BlobInputNode,
} from './nodes/input-nodes';
import {TableOutputNode, ValueOutputNode} from './nodes/output-nodes';
import {
  SelectColumnNode, SelectColumnsNode, SelectTableNode, AddTableViewNode, LogNode, InfoNode,
  WarningNode, ToStringNode, FromJsonNode, ToJsonNode,
  ConstStringNode, ConstIntNode, ConstDoubleNode, ConstBoolNode, ConstListNode,
} from './nodes/utility-nodes';
import {
  EqualsNode, NotEqualsNode, GreaterThanNode, GreaterOrEqualNode,
  LessThanNode, LessOrEqualNode, ContainsNode, StartsWithNode, EndsWithNode, IsNullNode,
} from './nodes/comparison-nodes';
import {BreakpointNode} from './nodes/breakpoint-node';
import {ViewerNode, CORE_VIEWER_SPECS, genericViewerSpec, VIEWER_TYPE_PREFIX, ViewerSpec} from './nodes/viewer-node';
import {
  getRole, getTags, getPackageName, getFuncDisplayName, getFuncQualifiedName, safeGet,
} from '../utils/dart-proxy-utils';
import {propertyNameToFriendly} from '../utils/naming';
import {INCLUDED_FUNC_NQNAMES} from './included-funcs';
import {builtinNodeDesc} from './builtin-catalog';

export interface FuncInfo {
  func: DG.Func;
  name: string;
  role: string | null;
  tags: string[];
  packageName: string;
  nodeTypeName: string;
  /** Lazily computed by {@link funcCategory}. */
  category?: string;
}

/** A saved Flow — a `DG.Script` whose language is `flow`. */
export function isWorkflowFunc(func: DG.Func): boolean {
  try {
    // `language` is typed as the core-language union, which doesn't know `flow` — compare the raw value.
    return func instanceof DG.Script && String((func as DG.Script).language) === 'flow';
  } catch {
    return false;
  }
}

/** The "what it does" category for a catalog function, cached on the FuncInfo (Dart-proxy reads aren't free). */
export function funcCategory(info: FuncInfo): string {
  if (info.category) return info.category;
  let cat = 'Other';
  try {
    if (isWorkflowFunc(info.func))
      cat = 'Workflows';
    else {
      const ins = info.func.inputs.map((p) => String(p.propertyType));
      const outs = info.func.outputs.map((p) => String(p.propertyType));
      cat = domainCategory(info.packageName, ins) ?? categorizeBySignature(ins, outs, info.role);
    }
  } catch {/* Dart proxy edge cases — keep 'Other' */}
  info.category = cat;
  return cat;
}

type Factory = () => FlowNode;

const FACTORIES = new Map<string, Factory>();
let funcRegistry: FuncInfo[] = [];
let registered = false;

function register(name: string, factory: Factory): void {
  FACTORIES.set(name, factory);
}

/** Dev/test packages whose queries would otherwise flood the Queries pane. */
const DEV_TEST_PACKAGES = new Set<string>([
  'Dbtests', 'ApiTests', 'UiTests', 'DevTools', 'Tutorials', 'ApiSamples', 'UsageAnalysis',
]);

export function shouldIncludeFunc(func: DG.Func): boolean {
  try {
    const include = safeGet(func.options, 'includeInFlow');
    if (include === false || String(include).toLowerCase() === 'false') return false;
    if (include === true || String(include).toLowerCase() === 'true') return true;
  } catch {/* options can throw on odd Dart proxies — fall through */}

  if (isWorkflowFunc(func)) return true;
  if (func instanceof DG.DataQuery)
    return !DEV_TEST_PACKAGES.has(getPackageName(func));

  try {
    return INCLUDED_FUNC_NQNAMES.has(func.nqName);
  } catch {/* nqName can throw on odd Dart proxies */}
  return false;
}

/** Viewer node types for the Viewers toolbox pane. */
export interface ViewerNodeType {label: string; nodeTypeName: string; core: boolean;}
export const VIEWER_NODE_TYPES: ViewerNodeType[] = [];

function registerViewerSpec(spec: ViewerSpec, core: boolean): void {
  const typeName = `${VIEWER_TYPE_PREFIX}${spec.label}`;
  if (FACTORIES.has(typeName)) return;
  register(typeName, () => new ViewerNode(spec));
  VIEWER_NODE_TYPES.push({label: spec.label, nodeTypeName: typeName, core});
}

export function registerBuiltinNodes(): void {
  if (registered) return;
  registered = true;

  register('Inputs/Table Input', () => new TableInputNode());
  register('Inputs/Column Input', () => new ColumnInputNode());
  register('Inputs/Column List Input', () => new ColumnListInputNode());
  register('Inputs/String Input', () => new StringInputNode());
  register('Inputs/Sketcher Input', () => new MoleculeInputNode());
  register('Inputs/Helm Input', () => new HelmInputNode());
  register('Inputs/Number Input', () => new NumberInputNode());
  register('Inputs/Int Input', () => new IntInputNode());
  register('Inputs/Boolean Input', () => new BooleanInputNode());
  register('Inputs/DateTime Input', () => new DateTimeInputNode());
  register('Inputs/File Input', () => new FileInputNode());
  register('Inputs/Map Input', () => new MapInputNode());
  register('Inputs/Dynamic Input', () => new DynamicInputNode());
  register('Inputs/String List Input', () => new StringListInputNode());
  register('Inputs/Blob Input', () => new BlobInputNode());

  register('Outputs/Table Output', () => new TableOutputNode());
  register('Outputs/Value Output', () => new ValueOutputNode());

  register('Utilities/Select Column', () => new SelectColumnNode());
  register('Utilities/Select Columns', () => new SelectColumnsNode());
  register('Utilities/Select Table', () => new SelectTableNode());
  register('Utilities/Add Table View', () => new AddTableViewNode());
  register('Utilities/Log', () => new LogNode());
  register('Utilities/Info', () => new InfoNode());
  register('Utilities/Warning', () => new WarningNode());
  register('Utilities/ToString', () => new ToStringNode());
  register('Utilities/FromJSON', () => new FromJsonNode());
  register('Utilities/ToJSON', () => new ToJsonNode());

  register('Constants/String', () => new ConstStringNode());
  register('Constants/Int', () => new ConstIntNode());
  register('Constants/Double', () => new ConstDoubleNode());
  register('Constants/Boolean', () => new ConstBoolNode());
  register('Constants/List', () => new ConstListNode());

  register('Comparisons/Equals (==)', () => new EqualsNode());
  register('Comparisons/Not Equals (!=)', () => new NotEqualsNode());
  register('Comparisons/Greater Than (>)', () => new GreaterThanNode());
  register('Comparisons/Greater Or Equal (>=)', () => new GreaterOrEqualNode());
  register('Comparisons/Less Than (<)', () => new LessThanNode());
  register('Comparisons/Less Or Equal (<=)', () => new LessOrEqualNode());
  register('Comparisons/Contains', () => new ContainsNode());
  register('Comparisons/Starts With', () => new StartsWithNode());
  register('Comparisons/Ends With', () => new EndsWithNode());
  register('Comparisons/Is Null', () => new IsNullNode());

  register('Debug/Breakpoint', () => new BreakpointNode());

  // Registered synchronously so saved flows with viewer nodes deserialize before registerAllFunctions runs.
  for (const spec of CORE_VIEWER_SPECS) registerViewerSpec(spec, true);
}

function registerDiscoveredViewers(): void {
  try {
    const taken = new Set(VIEWER_NODE_TYPES.map((v) => v.label.toLowerCase()));
    const found: string[] = [];
    for (const f of DG.Func.find({meta: {role: 'viewer'}})) {
      const name = String((f as DG.Func).friendlyName ?? f.name ?? '').trim();
      if (!name || taken.has(name.toLowerCase())) continue;
      taken.add(name.toLowerCase());
      found.push(name);
    }
    for (const name of found.sort((a, b) => a.localeCompare(b)))
      registerViewerSpec(genericViewerSpec(name, name), false);
  } catch (e) {
    console.warn('FuncFlow: viewer discovery failed', e);
  }
}

let funcsRegistered = false;

export function registerAllFunctions(): FuncInfo[] {
  if (funcsRegistered) return funcRegistry;
  funcsRegistered = true;

  registerDiscoveredViewers();

  const allFuncs = DG.Func.find({});
  funcRegistry = [];
  const seen = new Set<string>();

  for (const func of allFuncs) {
    try {
      if (func.inputs.length === 0 && func.outputs.length === 0) continue;
      if (!shouldIncludeFunc(func)) continue;
      const role = getRole(func);
      const tags = getTags(func);
      const pkgName = getPackageName(func);

      const category = role || 'Uncategorized';

      let typeName = `DG Functions/${category}/${func.name}`;
      if (seen.has(typeName))
        typeName = `DG Functions/${category}/${pkgName}:${func.name}`;
      if (seen.has(typeName)) continue;
      seen.add(typeName);

      const capturedFunc = func;
      register(typeName, () => new FuncNode(capturedFunc));

      funcRegistry.push({
        func, name: getFuncDisplayName(func) || func.name,
        role, tags, packageName: pkgName, nodeTypeName: typeName,
      });
    } catch {
      // Dart proxy introspection can throw.
    }
  }

  registerVariableFuncs();

  return funcRegistry;
}

/** SetVar/GetVar are excluded by the catalog rules, yet saved flows depend on them — register unconditionally. */
function registerVariableFuncs(): void {
  for (const name of ['SetVar', 'GetVar']) {
    try {
      const func = DG.Func.find({name})[0];
      if (func) ensureFuncNodeType(func);
    } catch {
      // No live backend.
    }
  }
}

export function getRegisteredFuncs(): FuncInfo[] {
  return funcRegistry;
}

export function ensureFunctionsRegistered(): void {
  registerBuiltinNodes();
  registerAllFunctions();
}

/** Node type name for a DG.Func, registering a factory on the fly for off-catalog funcs;
 *  matches by qualified name so a freshly parsed Dart instance finds its registered twin. */
export function ensureFuncNodeType(func: DG.Func): string {
  registerAllFunctions();

  const qName = getFuncQualifiedName(func).toLowerCase();
  for (const info of funcRegistry) {
    if (getFuncQualifiedName(info.func).toLowerCase() === qName)
      return info.nodeTypeName;
  }

  const role = getRole(func);
  const pkgName = getPackageName(func);
  const category = role || 'Uncategorized';
  let typeName = `DG Functions/${category}/${func.name}`;
  if (FACTORIES.has(typeName))
    typeName = `DG Functions/${category}/${pkgName}:${func.name}`;
  for (let i = 2; FACTORIES.has(typeName); i++)
    typeName = `DG Functions/${category}/${pkgName}:${func.name}#${i}`;

  const capturedFunc = func;
  register(typeName, () => new FuncNode(capturedFunc));
  funcRegistry.push({
    func, name: getFuncDisplayName(func) || func.name,
    role, tags: getTags(func), packageName: pkgName, nodeTypeName: typeName,
  });
  return typeName;
}

let queryFuncsPromise: Promise<FuncInfo[]> | null = null;

/** Queries-pane catalog: `DG.Func.find({})` misses queries, so load the authoritative
 *  server list (cached; a failed load clears the cache so a later render retries). */
export function loadQueryFuncs(): Promise<FuncInfo[]> {
  if (queryFuncsPromise) return queryFuncsPromise;
  queryFuncsPromise = (async () => {
    const queries = await grok.dapi.queries.list();
    const result: FuncInfo[] = [];
    for (const q of queries) {
      try {
        if (!q.connection) continue;
        const pkgName = getPackageName(q);
        if (DEV_TEST_PACKAGES.has(pkgName)) continue;
        const typeName = ensureFuncNodeType(q);
        result.push({
          func: q, name: getFuncDisplayName(q) || q.name,
          role: getRole(q), tags: getTags(q), packageName: pkgName, nodeTypeName: typeName,
        });
      } catch {
        // Dart proxy introspection can throw.
      }
    }
    return result;
  })().catch((e) => {
    queryFuncsPromise = null;
    throw e;
  });
  return queryFuncsPromise;
}

/** Instantiates a registered node type; stamps `dgTypeName` for the serializer. */
export function createNode(typeName: string): FlowNode | null {
  const factory = FACTORIES.get(typeName);
  if (!factory) return null;
  const node = factory();
  node.dgTypeName = typeName;
  humanizeSlotLabels(node);
  addExecPorts(node);
  return node;
}

/** A slot label still equal to its raw key is a fallback, not a caption — humanize it;
 *  declared captions and the `→` pass-through markers differ from their key and stay untouched. */
function humanizeSlotLabels(node: FlowNode): void {
  const ports = [
    ...Object.entries(node.inputs as Record<string, {label?: string} | undefined>),
    ...Object.entries(node.outputs as Record<string, {label?: string} | undefined>),
  ];
  for (const [key, port] of ports) {
    if (port && port.label === key)
      port.label = propertyNameToFriendly(key);
  }
}

/** Exec ports go after the node's own ports so `passthroughCount` indexing and the
 *  factory type probes (which must never see an `order` slot) are unaffected. */
function addExecPorts(node: FlowNode): void {
  node.addInput(EXEC_IN_KEY, new ClassicPreset.Input(getSocket(ORDER_SOCKET_TYPE), 'before', true));
  node.addOutput(EXEC_OUT_KEY, new ClassicPreset.Output(getSocket(ORDER_SOCKET_TYPE), 'after', true));
}

export function getRegisteredTypeNames(): string[] {
  return Array.from(FACTORIES.keys());
}

/** One sample node per type is constructed to read its input types, then cached —
 *  the suggestion menu probes hundreds of factories. */
const _sampleInputTypesCache = new Map<string, string[]>();

export function getInputTypesForType(typeName: string): string[] {
  let cached = _sampleInputTypesCache.get(typeName);
  if (cached !== undefined) return cached;
  const factory = FACTORIES.get(typeName);
  if (!factory) return [];
  try {
    const sample = factory();
    cached = (Object.values(sample.inputs) as Array<{socket: TypedSocket} | undefined>)
      .map((i) => i?.socket.dgType ?? 'dynamic');
  } catch {
    cached = [];
  }
  _sampleInputTypesCache.set(typeName, cached);
  return cached;
}

/** Output slot dgTypes split into real vs pass-through (`__pt`) — real producers rank above threaders. */
const _sampleOutputTypesCache = new Map<string, {real: string[]; passthrough: string[]}>();

export function getOutputTypesForType(typeName: string): {real: string[]; passthrough: string[]} {
  let cached = _sampleOutputTypesCache.get(typeName);
  if (cached !== undefined) return cached;
  const factory = FACTORIES.get(typeName);
  if (!factory) return {real: [], passthrough: []};
  try {
    const sample = factory();
    cached = {real: [], passthrough: []};
    for (const [key, out] of Object.entries(sample.outputs) as Array<[string, {socket: TypedSocket} | undefined]>) {
      if (!out) continue;
      (key.endsWith('__pt') ? cached.passthrough : cached.real).push(out.socket.dgType);
    }
  } catch {
    cached = {real: [], passthrough: []};
  }
  _sampleOutputTypesCache.set(typeName, cached);
  return cached;
}

/** A semType-tagged input node (Sketcher = String Input + `semType: Molecule`) is a specialization —
 *  a bare-type drag must still land on the general one. */
const _specializedInputCache = new Map<string, boolean>();

function isSpecializedInput(typeName: string): boolean {
  let cached = _specializedInputCache.get(typeName);
  if (cached !== undefined) return cached;
  try {
    cached = String(FACTORIES.get(typeName)?.().properties['semType'] ?? '').trim().length > 0;
  } catch {
    cached = false;
  }
  _specializedInputCache.set(typeName, cached);
  return cached;
}

export interface CompatibleNodeType {
  typeName: string;
  label: string;
  isBuiltin: boolean;
  description?: string;
  /** Exact type match (vs a `dynamic`/`object` wildcard) — ranks first in its tier. */
  exact?: boolean;
  /** Reverse menu only: a real output matches (vs pass-through only). */
  realOutput?: boolean;
  /** Lower-cased search haystack — the same texts the toolbox search matches. */
  searchText?: string;
  /** Set when the toolbox suggestion engine also recommends this type — such items lead the menu. */
  reason?: string;
}

function candidateSearchText(typeName: string, label: string, info?: FuncInfo): string {
  const parts = [label, typeName];
  if (info) {
    parts.push(info.name, info.role ?? '', info.packageName, ...info.tags);
    try {
      parts.push(String(info.func.name || ''), String(info.func.friendlyName || ''),
        String(info.func.description || ''));
    } catch {/* Dart proxy */}
  }
  else
    parts.push(builtinNodeDesc(typeName));
  return parts.filter(Boolean).join(' ').toLowerCase();
}

export function candidateMatchesQuery(c: CompatibleNodeType, query: string): boolean {
  const q = query.toLowerCase().trim();
  if (!q) return true;
  const hay = c.searchText ?? `${c.label} ${c.typeName}`.toLowerCase();
  return hay.includes(q) || hay.replace(/\s+/g, '').includes(q.replace(/\s+/g, ''));
}

/** One toolbox-suggestion-engine pick, projected for the drag-out menu. */
export interface SocketSuggestion {
  typeName: string;
  reason: string;
  prefill?: Record<string, unknown>;
}

/** Suggestion-engine picks lead in engine order (with reasons); a suggested type not in `candidates` is dropped. */
export function prioritizeCandidates(
  candidates: CompatibleNodeType[], suggested: SocketSuggestion[],
): CompatibleNodeType[] {
  if (suggested.length === 0) return candidates;
  const order = new Map<string, {idx: number; s: SocketSuggestion}>();
  suggested.forEach((s, idx) => {
    if (!order.has(s.typeName)) order.set(s.typeName, {idx, s});
  });
  const lead: CompatibleNodeType[] = [];
  const rest: CompatibleNodeType[] = [];
  for (const c of candidates) {
    const hit = order.get(c.typeName);
    if (hit) lead.push({...c, reason: hit.s.reason});
    else rest.push(c);
  }
  lead.sort((a, b) => order.get(a.typeName)!.idx - order.get(b.typeName)!.idx);
  return [...lead, ...rest];
}

/** Built-ins show their trailing typeName segment; DG funcs the friendly name + category
 *  (never the raw role segment, which reads `AddNewColumn (Uncategorized)` for most funcs). */
function labelForTypeName(typeName: string, info?: FuncInfo): string {
  if (info) return `${info.name}  (${funcCategory(info)})`;
  const parts = typeName.split('/');
  if (parts[0] === 'DG Functions' && parts.length >= 3)
    return `${parts[parts.length - 1]}  (${parts[1]})`;
  return parts[parts.length - 1];
}

/** Usage-proxy for common next-steps (no client-side telemetry); matched on the lower-cased simple name. */
const COMMON_NEXT_FUNCS = new Set([
  'jointables', 'linktables', 'addnewcolumn', 'addnewcolumnlist', 'aggregate',
  'filterrows', 'extractrows', 'extractcolumns', 'pivot', 'unpivot', 'splitcolumn',
  'addtableview', 'clonetable', 'renamecolumn', 'changecolumnstype',
]);

function simpleFuncName(typeName: string): string {
  const last = typeName.split('/').pop() ?? typeName;
  return (last.split(':').pop() ?? last).toLowerCase();
}

export function isCommonNextFunc(typeName: string): boolean {
  return COMMON_NEXT_FUNCS.has(simpleFuncName(typeName));
}

/** Canvas context the drag-out suggestion menus rank against; all optional. */
export interface SuggestionContext {
  /** Package of the node the drag started from ('' / undefined for built-ins). */
  sourcePackageName?: string | null;
  /** Source packages of every node already on the canvas. */
  graphPackageNames?: Iterable<string>;
  /** Simple function names already used on the canvas (any case). */
  graphFuncNames?: Iterable<string>;
}

/** Node types with an input compatible with `sourceType`, ranked for the drag-out suggestion menu. */
export function findNodeTypesAcceptingInput(
  sourceType: string, context?: SuggestionContext,
): CompatibleNodeType[] {
  const matches: Array<CompatibleNodeType & {exact: boolean}> = [];
  const infoByTypeName = new Map(funcRegistry.map((f) => [f.nodeTypeName, f]));
  for (const typeName of FACTORIES.keys()) {
    const inputTypes = getInputTypesForType(typeName);
    if (inputTypes.length === 0) continue;
    if (!inputTypes.some((t) => areTypesCompatible(sourceType, t))) continue;
    const info = infoByTypeName.get(typeName);
    const label = labelForTypeName(typeName, info);
    matches.push({
      typeName,
      label,
      isBuiltin: !typeName.startsWith('DG Functions/'),
      exact: inputTypes.includes(sourceType),
      searchText: candidateSearchText(typeName, label, info),
    });
  }

  const {preferredDomains, usedFuncs} = contextBoosts(context);

  const rank = (t: CompatibleNodeType): number => {
    if (t.typeName === 'Outputs/Value Output') return 0;
    const info = infoByTypeName.get(t.typeName);
    if (info && preferredDomains.size > 0 && preferredDomains.has(funcCategory(info))) return 1;
    if (COMMON_NEXT_FUNCS.has(simpleFuncName(t.typeName))) return 2;
    return t.isBuiltin ? 3 : 4;
  };
  const used = (t: CompatibleNodeType): number => (usedFuncs.has(simpleFuncName(t.typeName)) ? 0 : 1);
  matches.sort((a, b) =>
    rank(a) - rank(b) ||
    Number(b.exact) - Number(a.exact) ||
    used(a) - used(b) ||
    a.label.localeCompare(b.label));
  return matches;
}

/** Ranking boosts from the canvas context: the domain in play and the funcs already used. */
function contextBoosts(context?: SuggestionContext): {preferredDomains: Set<string>; usedFuncs: Set<string>} {
  const preferredDomains = new Set<string>();
  const srcDomain = domainSection(context?.sourcePackageName || undefined);
  if (srcDomain) preferredDomains.add(srcDomain);
  else {
    for (const p of context?.graphPackageNames ?? []) {
      const d = domainSection(p);
      if (d) preferredDomains.add(d);
    }
  }
  const usedFuncs = new Set<string>();
  for (const n of context?.graphFuncNames ?? []) usedFuncs.add(n.toLowerCase());
  return {preferredDomains, usedFuncs};
}

/** Node types with an output (or pass-through) compatible with `targetType` —
 *  the reverse "what produces this?" suggestion menu. */
export function findNodeTypesProducingOutput(
  targetType: string, context?: SuggestionContext,
): CompatibleNodeType[] {
  const matches: Array<CompatibleNodeType & {exact: boolean; realOutput: boolean}> = [];
  const infoByTypeName = new Map(funcRegistry.map((f) => [f.nodeTypeName, f]));
  for (const typeName of FACTORIES.keys()) {
    const {real, passthrough} = getOutputTypesForType(typeName);
    const realCompat = real.some((t) => areTypesCompatible(t, targetType));
    if (!realCompat && !passthrough.some((t) => areTypesCompatible(t, targetType))) continue;
    const info = infoByTypeName.get(typeName);
    const label = labelForTypeName(typeName, info);
    matches.push({
      typeName,
      label,
      isBuiltin: !typeName.startsWith('DG Functions/'),
      realOutput: realCompat,
      exact: (realCompat ? real : passthrough).includes(targetType),
      searchText: candidateSearchText(typeName, label, info),
    });
  }

  const {preferredDomains, usedFuncs} = contextBoosts(context);

  const rank = (t: CompatibleNodeType & {exact: boolean; realOutput: boolean}): number => {
    if (t.typeName.startsWith('Inputs/') && t.realOutput && t.exact) return 0;
    const info = infoByTypeName.get(t.typeName);
    if (info && preferredDomains.size > 0 && preferredDomains.has(funcCategory(info))) return 1;
    if (COMMON_NEXT_FUNCS.has(simpleFuncName(t.typeName))) return 2;
    if (info && funcCategory(info) === 'Data Sources') return 2;
    return t.isBuiltin ? 3 : 4;
  };
  const used = (t: CompatibleNodeType): number => (usedFuncs.has(simpleFuncName(t.typeName)) ? 0 : 1);
  const specialized = (t: CompatibleNodeType): number => Number(isSpecializedInput(t.typeName));
  matches.sort((a, b) =>
    rank(a) - rank(b) ||
    Number(b.realOutput) - Number(a.realOutput) ||
    Number(b.exact) - Number(a.exact) ||
    specialized(a) - specialized(b) ||
    used(a) - used(b) ||
    a.label.localeCompare(b.label));
  return matches;
}
