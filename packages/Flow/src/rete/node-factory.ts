/** Node-type registry — replaces LiteGraph's `registerNodeType` /
 *  `createNode` pair with a plain factory map.
 *
 *  Each node type has a string name like `Inputs/Table Input` (kept for
 *  parity with the old serialization keys, function-browser categories, etc.)
 *  and a zero-arg factory that returns a fresh `FlowNode` instance. */

import * as DG from 'datagrok-api/dg';
import {ClassicPreset} from 'rete';
import {FlowNode, EXEC_IN_KEY, EXEC_OUT_KEY, ORDER_SOCKET_TYPE} from './scheme';
import {FuncNode} from './nodes/func-node';
import {TypedSocket, getSocket} from './sockets';
import {areTypesCompatible, categorizeBySignature, domainCategory, domainSection} from '../types/type-map';

import {
  TableInputNode, ColumnInputNode, ColumnListInputNode, StringInputNode,
  NumberInputNode, IntInputNode, BooleanInputNode, DateTimeInputNode,
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

export interface FuncInfo {
  func: DG.Func;
  name: string;
  role: string | null;
  tags: string[];
  packageName: string;
  nodeTypeName: string;
  /** Lazily-computed "what it does" category (see {@link funcCategory}). */
  category?: string;
}

/** A saved Flow — a `DG.Script` whose language is `flow`. Runnable like any
 *  function, but grouped under its own 'Workflows' section in every browser
 *  grouping (its signature would otherwise misfile it under Data Sources). */
export function isWorkflowFunc(func: DG.Func): boolean {
  try {
    // `language` is typed as the core-language union, which doesn't know
    // about `flow` — compare the raw value.
    return func instanceof DG.Script && String((func as DG.Script).language) === 'flow';
  } catch {
    return false;
  }
}

/** The "what it does" category for a catalog function — workflows (saved flows)
 *  first, then domain (chem/bio) for operations on data, else the signature
 *  category; the same routing as the browser's `categorizeFunc` and the node
 *  title-bar coloring. Cached on the FuncInfo (Dart-proxy input/output reads
 *  aren't free). */
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

/** Dev/test packages whose *queries* would otherwise flood the Queries pane —
 *  queries are included by kind, and Dbtests alone ships ~290 test queries. */
const DEV_TEST_PACKAGES = new Set<string>([
  'Dbtests', 'ApiTests', 'UiTests', 'DevTools', 'Tutorials', 'ApiSamples', 'UsageAnalysis',
]);

/** Allowlist-based inclusion — the opposite of the old deny pipeline. A
 *  function enters the Flow catalog only when:
 *   - it declares **`meta.includeInFlow: true`** (the author opt-in — meta.*
 *     surfaces as `func.options`; the value may arrive as a string), or
 *   - it is a **saved flow** (`isWorkflowFunc`) or a **query**
 *     (`DG.DataQuery`, except from a dev/test package) — user artifacts that
 *     can't be enumerated statically (own toolbox panes: Workflows section /
 *     per-connection Queries pane), or
 *   - its `nqName` is in [`INCLUDED_FUNC_NQNAMES`](./included-funcs.ts) — the
 *     curated include list, frozen from the catalog the old filters produced.
 *  `meta.includeInFlow: false` opts out unconditionally (wins over the list).
 *  Everything else — dev/test packages, panels, sketchers, dialog wrappers,
 *  scalar helpers — is simply *not on the list*. */
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

/** Viewer node types for the Viewers toolbox pane (core first, then discovered
 *  package viewers). Populated by `registerBuiltinNodes` + `registerAllFunctions`. */
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

  // Inputs
  register('Inputs/Table Input', () => new TableInputNode());
  register('Inputs/Column Input', () => new ColumnInputNode());
  register('Inputs/Column List Input', () => new ColumnListInputNode());
  register('Inputs/String Input', () => new StringInputNode());
  register('Inputs/Number Input', () => new NumberInputNode());
  register('Inputs/Int Input', () => new IntInputNode());
  register('Inputs/Boolean Input', () => new BooleanInputNode());
  register('Inputs/DateTime Input', () => new DateTimeInputNode());
  register('Inputs/File Input', () => new FileInputNode());
  register('Inputs/Map Input', () => new MapInputNode());
  register('Inputs/Dynamic Input', () => new DynamicInputNode());
  register('Inputs/String List Input', () => new StringListInputNode());
  register('Inputs/Blob Input', () => new BlobInputNode());

  // Outputs
  register('Outputs/Table Output', () => new TableOutputNode());
  register('Outputs/Value Output', () => new ValueOutputNode());

  // Utilities
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

  // Constants
  register('Constants/String', () => new ConstStringNode());
  register('Constants/Int', () => new ConstIntNode());
  register('Constants/Double', () => new ConstDoubleNode());
  register('Constants/Boolean', () => new ConstBoolNode());
  register('Constants/List', () => new ConstListNode());

  // Comparisons
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

  // Debug
  register('Debug/Breakpoint', () => new BreakpointNode());

  // Core viewers (sync — no catalog lookup — so saved flows with viewer nodes
  // always deserialize, even before `registerAllFunctions` runs).
  for (const spec of CORE_VIEWER_SPECS) registerViewerSpec(spec, true);
}

/** Discover non-core package viewers (role=viewer) and register a generic
 *  viewer node per distinct friendly name. Idempotent; needs the live catalog. */
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

/** Discover all DG.Func instances, register a factory per function, return
 *  the catalog used by the function-browser. */
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

      // Capture func by reference for the factory.
      const capturedFunc = func;
      register(typeName, () => new FuncNode(capturedFunc));

      funcRegistry.push({
        func, name: getFuncDisplayName(func) || func.name,
        role, tags, packageName: pkgName, nodeTypeName: typeName,
      });
    } catch {
      // Skip funcs that fail to introspect (Dart proxy edge cases).
    }
  }

  registerVariableFuncs();

  return funcRegistry;
}

/** SetVar / GetVar fall to the primitive-only exclusion rule, yet saved flows
 *  depend on them (every imported creation script terminates in SetVar nodes,
 *  and Flow treats SetVar as an output). Register them unconditionally —
 *  otherwise a saved .ffjson with SetVar/GetVar nodes deserializes only in a
 *  session where a creation-script import happened to register them first. */
function registerVariableFuncs(): void {
  for (const name of ['SetVar', 'GetVar']) {
    try {
      const func = DG.Func.find({name})[0];
      if (func) ensureFuncNodeType(func);
    } catch {
      // No live backend / lookup failure — nothing to register.
    }
  }
}

export function getRegisteredFuncs(): FuncInfo[] {
  return funcRegistry;
}

/** Deterministic one-shot registration of everything `createNode` may need:
 *  built-ins plus the DG.Func catalog. Load paths (deserializer, entity open,
 *  headless compile) call this instead of racing the view's deferred timer —
 *  both underlying calls are synchronous and idempotent. */
export function ensureFunctionsRegistered(): void {
  registerBuiltinNodes();
  registerAllFunctions();
}

/** Node type name for a DG.Func, registering a factory on the fly when the
 *  function is not in the catalog (i.e. not on the include list). Used by the
 *  creation-script importer, where parsed funcs must always yield a node —
 *  and a `dgTypeName` the serializer can persist. Matching is by qualified
 *  name, so a freshly parsed Dart instance finds its registered twin. */
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

/** Instantiate a registered node type by name. Returns null if unknown.
 *  Stamps `dgTypeName` on the new node so the serializer can persist it. */
export function createNode(typeName: string): FlowNode | null {
  const factory = FACTORIES.get(typeName);
  if (!factory) return null;
  const node = factory();
  node.dgTypeName = typeName;
  humanizeSlotLabels(node);
  addExecPorts(node);
  return node;
}

/** A slot label that still equals its raw key is a fallback, not a caption —
 *  humanize it (core's `propertyNameToFriendly`) so 'table' reads 'Table' and
 *  a caption-less func output 'maxNumOfSomething' reads 'Max Num Of Something',
 *  mirroring `ui.input.forProperty`. Runs before the exec ports are added
 *  (their labels are 'before'/'after' anyway); declared captions and the `→`
 *  pass-through markers differ from their key and stay untouched. Keys —
 *  connections, `inputValues`, compilation — are never affected. */
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

/** Add the execution-ordering port pair to a node: an exec-in (accepts many
 *  predecessors) and an exec-out. Added after the node's own ports so func
 *  `passthroughCount` indexing and the data-port rows are unaffected. The
 *  factory builders are left untouched (so the suggestion-menu type probe in
 *  `getInputTypesForType` never sees an `order` slot). */
function addExecPorts(node: FlowNode): void {
  node.addInput(EXEC_IN_KEY, new ClassicPreset.Input(getSocket(ORDER_SOCKET_TYPE), 'before', true));
  node.addOutput(EXEC_OUT_KEY, new ClassicPreset.Output(getSocket(ORDER_SOCKET_TYPE), 'after', true));
}

/** All registered type names, mostly for debugging / completeness checks. */
export function getRegisteredTypeNames(): string[] {
  return Array.from(FACTORIES.keys());
}

// ---------- input-type cache + suggestion helper ----------

/** Cache of (typeName → input slot dgTypes), populated lazily. We construct
 *  one sample of each node to read its input socket types and never again —
 *  important for the suggestion menu where we may probe hundreds of factories. */
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

/** Cache of (typeName → output slot dgTypes), split into real outputs vs
 *  pass-throughs (`__pt` keys) — the reverse suggestion menu ("what produces
 *  this?") ranks real producers above threaders. Same one-sample probe as
 *  `getInputTypesForType`. */
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

export interface CompatibleNodeType {
  typeName: string;
  label: string;
  isBuiltin: boolean;
  description?: string;
  /** Whether a slot matches the dragged type exactly (vs a `dynamic`/`object`
   *  wildcard) — exact matches rank first in their tier. */
  exact?: boolean;
  /** Reverse menu only: whether a **real** output matches (vs pass-through
   *  only) — real producers rank above threaders in their tier. */
  realOutput?: boolean;
}

/** Display label for a suggestion-menu candidate. Built-ins show their trailing
 *  typeName segment; DG functions show the same **friendly name** the toolbox
 *  uses, with the "what it does" category in parentheses for orientation —
 *  NOT the raw func name / role segment baked into the typeName (which reads
 *  as `AddNewColumn (Uncategorized)` for the role-less majority). */
function labelForTypeName(typeName: string, info?: FuncInfo): string {
  if (info) return `${info.name}  (${funcCategory(info)})`;
  const parts = typeName.split('/');
  if (parts[0] === 'DG Functions' && parts.length >= 3)
    return `${parts[parts.length - 1]}  (${parts[1]})`;
  return parts[parts.length - 1];
}

/** The data-pipeline steps a scientist reaches for most. Used as a usage-proxy
 *  to float common next-steps to the top of the drag-out suggestion menu (real
 *  telemetry isn't available client-side). Matched on the simple function name,
 *  lower-cased. */
const COMMON_NEXT_FUNCS = new Set([
  'jointables', 'linktables', 'addnewcolumn', 'addnewcolumnlist', 'aggregate',
  'filterrows', 'extractrows', 'extractcolumns', 'pivot', 'unpivot', 'splitcolumn',
  'addtableview', 'clonetable', 'renamecolumn', 'changecolumnstype',
]);

/** Simple (trailing) function name of a typeName, lower-cased — `DG Functions/
 *  Transform/Pkg:JoinTables` → `jointables`. */
function simpleFuncName(typeName: string): string {
  const last = typeName.split('/').pop() ?? typeName;
  return (last.split(':').pop() ?? last).toLowerCase();
}

/** Whether a node type is one of the common next pipeline steps (Join,
 *  AddNewColumn, Aggregate, …) — shared by the drag-out menu ranking and the
 *  toolbox suggestion engine. */
export function isCommonNextFunc(typeName: string): boolean {
  return COMMON_NEXT_FUNCS.has(simpleFuncName(typeName));
}

/** Canvas context the drag-out suggestion menu ranks against — what the drag
 *  came from and what the user is already building with. All optional; without
 *  it the ranking degrades to the context-free order. */
export interface SuggestionContext {
  /** Package of the node the drag started from ('' / undefined for built-ins). */
  sourcePackageName?: string | null;
  /** Source packages of every node already on the canvas. */
  graphPackageNames?: Iterable<string>;
  /** Simple function names already used on the canvas (any case). */
  graphFuncNames?: Iterable<string>;
}

/** All registered node types whose inputs include at least one slot
 *  type-compatible with `sourceType`. Used by the drag-output suggestion menu.
 *
 *  Ranked by discovery heuristics (lower tier = higher in the list):
 *  1. **Value Output** — the universal "expose this value" terminal.
 *  2. **The science you're doing** — Cheminformatics / Bioinformatics functions
 *     when the drag *source* is from that domain's package (else when the
 *     canvas already holds nodes of that domain): dragging out of Chemical
 *     Properties should offer chem next-steps first, not bury them
 *     alphabetically.
 *  3. **Common next-step funcs** (`COMMON_NEXT_FUNCS` — Join, AddNewColumn,
 *     Aggregate, Filter, …) — the usage-proxy for core table plumbing.
 *  4. Other built-ins (viewers, outputs, utilities).
 *  5. The remaining DG functions.
 *
 *  Within a tier: an **exact type match** beats a wildcard (`dynamic`/`object`)
 *  acceptor — a `column` drag offers real column consumers before catch-alls;
 *  then functions **already used on the canvas** (pipelines repeat their ops);
 *  then alphabetical. */
export function findNodeTypesAcceptingInput(
  sourceType: string, context?: SuggestionContext,
): CompatibleNodeType[] {
  const matches: Array<CompatibleNodeType & {exact: boolean}> = [];
  // typeName → catalog entry, for friendly names + "what it does" categories.
  const infoByTypeName = new Map(funcRegistry.map((f) => [f.nodeTypeName, f]));
  for (const typeName of FACTORIES.keys()) {
    const inputTypes = getInputTypesForType(typeName);
    if (inputTypes.length === 0) continue;
    if (!inputTypes.some((t) => areTypesCompatible(sourceType, t))) continue;
    matches.push({
      typeName,
      label: labelForTypeName(typeName, infoByTypeName.get(typeName)),
      isBuiltin: !typeName.startsWith('DG Functions/'),
      exact: inputTypes.includes(sourceType),
    });
  }

  const {preferredDomains, usedFuncs} = contextBoosts(context);

  // Lower rank number = higher in the list.
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

/** Shared ranking boosts derived from the canvas context: the science in play
 *  (the drag-origin node's own domain wins; a domain-less origin — OpenFile, a
 *  utility — falls back to any domain already on the canvas) and the functions
 *  the user already reached for. */
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

/** All registered node types with an output — or a pass-through — compatible
 *  with `targetType`. The reverse suggestion menu: dragging an *input* to empty
 *  canvas asks "what produces this?".
 *
 *  Ranked like `findNodeTypesAcceptingInput`, with the terminals swapped for
 *  the upstream direction (lower tier = higher):
 *  1. The **matching Input node** (Table Input for a table drag, …) — the
 *     universal "make this a script parameter" producer.
 *  2. **The science in play** (same domain boost as the forward menu).
 *  3. **Data Sources** functions (OpenFile, queries, generators — the natural
 *     answer to "where does a table come from?") and the common table funcs.
 *  4. Other built-ins; 5. the remaining DG functions.
 *
 *  Within a tier: a **real output** match beats a pass-through-only threader
 *  (`realOutput`), then an **exact type match** beats a wildcard, then
 *  already-used-on-canvas funcs, then alphabetical. */
export function findNodeTypesProducingOutput(
  targetType: string, context?: SuggestionContext,
): CompatibleNodeType[] {
  const matches: Array<CompatibleNodeType & {exact: boolean; realOutput: boolean}> = [];
  const infoByTypeName = new Map(funcRegistry.map((f) => [f.nodeTypeName, f]));
  for (const typeName of FACTORIES.keys()) {
    const {real, passthrough} = getOutputTypesForType(typeName);
    const realCompat = real.some((t) => areTypesCompatible(t, targetType));
    if (!realCompat && !passthrough.some((t) => areTypesCompatible(t, targetType))) continue;
    matches.push({
      typeName,
      label: labelForTypeName(typeName, infoByTypeName.get(typeName)),
      isBuiltin: !typeName.startsWith('DG Functions/'),
      realOutput: realCompat,
      exact: (realCompat ? real : passthrough).includes(targetType),
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
  matches.sort((a, b) =>
    rank(a) - rank(b) ||
    Number(b.realOutput) - Number(a.realOutput) ||
    Number(b.exact) - Number(a.exact) ||
    used(a) - used(b) ||
    a.label.localeCompare(b.label));
  return matches;
}
