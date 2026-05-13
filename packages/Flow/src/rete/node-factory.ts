/** Node-type registry — replaces LiteGraph's `registerNodeType` /
 *  `createNode` pair with a plain factory map.
 *
 *  Each node type has a string name like `Inputs/Table Input` (kept for
 *  parity with the old serialization keys, function-browser categories, etc.)
 *  and a zero-arg factory that returns a fresh `FlowNode` instance. */

import * as DG from 'datagrok-api/dg';
import {FlowNode} from './scheme';
import {FuncNode} from './nodes/func-node';
import {TypedSocket} from './sockets';
import {areTypesCompatible} from '../types/type-map';

import {
  TableInputNode, ColumnInputNode, ColumnListInputNode, StringInputNode,
  NumberInputNode, IntInputNode, BooleanInputNode, DateTimeInputNode,
  FileInputNode, MapInputNode, DynamicInputNode, StringListInputNode, BlobInputNode,
} from './nodes/input-nodes';
import {TableOutputNode, ValueOutputNode} from './nodes/output-nodes';
import {
  SelectColumnNode, SelectColumnsNode, AddTableViewNode, LogNode, InfoNode,
  WarningNode, ToStringNode, FromJsonNode, ToJsonNode,
  ConstStringNode, ConstIntNode, ConstDoubleNode, ConstBoolNode, ConstListNode,
} from './nodes/utility-nodes';
import {
  EqualsNode, NotEqualsNode, GreaterThanNode, GreaterOrEqualNode,
  LessThanNode, LessOrEqualNode, ContainsNode, StartsWithNode, EndsWithNode, IsNullNode,
} from './nodes/comparison-nodes';
import {BreakpointNode} from './nodes/breakpoint-node';
import {getRole, getTags, getPackageName, getFuncDisplayName} from '../utils/dart-proxy-utils';

export interface FuncInfo {
  func: DG.Func;
  name: string;
  role: string | null;
  tags: string[];
  packageName: string;
  nodeTypeName: string;
}

type Factory = () => FlowNode;

const FACTORIES = new Map<string, Factory>();
let funcRegistry: FuncInfo[] = [];
let registered = false;

function register(name: string, factory: Factory): void {
  FACTORIES.set(name, factory);
}

/** Tags that mark a function for exclusion from the catalog. */
export const EXCLUDED_TAGS: string[] = [];

/** Roles that mark a function for exclusion. */
export const EXCLUDED_ROLES: string[] = [
  DG.FUNC_TYPES.APP, 'aiSearchProvider', 'antibodyNumbering', 'appTreeBrowser', 'canonicalizer',
  DG.FUNC_TYPES.CELL_RENDERER, 'dashboard', DG.FUNC_TYPES.FILE_VIEWER, DG.FUNC_TYPES.FILE_IMPORTER,
  DG.FUNC_TYPES.FILE_EXPORTER, DG.FUNC_TYPES.FILTER, DG.FUNC_TYPES.FOLDER_VIEWER, DG.FUNC_TYPES.MOLECULE_SKETCHER,
  'notationProviderConstructor', 'notationRefiner', 'packageSettingsEditor', 'searchProvider', 'semTypeDetector',
  'semValueExtractor', 'valueEditor',
];

function shouldExcludeFunc(role: string | null, tags: string[]): boolean {
  if (role && EXCLUDED_ROLES.includes(role)) return true;
  for (const tag of tags) if (EXCLUDED_TAGS.includes(tag)) return true;
  return false;
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
}

let funcsRegistered = false;

/** Discover all DG.Func instances, register a factory per function, return
 *  the catalog used by the function-browser. */
export function registerAllFunctions(): FuncInfo[] {
  if (funcsRegistered) return funcRegistry;
  funcsRegistered = true;

  const allFuncs = DG.Func.find({});
  funcRegistry = [];
  const seen = new Set<string>();

  for (const func of allFuncs) {
    try {
      if (func.inputs.length === 0 && func.outputs.length === 0) continue;
      const role = getRole(func);
      const tags = getTags(func);
      if (shouldExcludeFunc(role, tags)) continue;

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

  return funcRegistry;
}

export function getRegisteredFuncs(): FuncInfo[] {
  return funcRegistry;
}

/** Instantiate a registered node type by name. Returns null if unknown.
 *  Stamps `dgTypeName` on the new node so the serializer can persist it. */
export function createNode(typeName: string): FlowNode | null {
  const factory = FACTORIES.get(typeName);
  if (!factory) return null;
  const node = factory();
  node.dgTypeName = typeName;
  return node;
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

function getInputTypesForType(typeName: string): string[] {
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

export interface CompatibleNodeType {
  typeName: string;
  label: string;
  isBuiltin: boolean;
  description?: string;
}

/** Display label parsed from a typeName like `Inputs/Table Input` or
 *  `DG Functions/Transform/JoinTables` → the trailing segment, but for DG
 *  functions we also include the role for orientation. */
function labelForTypeName(typeName: string): string {
  const parts = typeName.split('/');
  if (parts[0] === 'DG Functions' && parts.length >= 3)
    return `${parts[parts.length - 1]}  (${parts[1]})`;
  return parts[parts.length - 1];
}

/** All registered node types whose inputs include at least one slot
 *  type-compatible with `sourceType`. Used by the drag-output suggestion
 *  menu — Value Output is sorted first because dynamic accepts everything. */
export function findNodeTypesAcceptingInput(sourceType: string): CompatibleNodeType[] {
  const matches: CompatibleNodeType[] = [];
  for (const typeName of FACTORIES.keys()) {
    const inputTypes = getInputTypesForType(typeName);
    if (inputTypes.length === 0) continue;
    if (!inputTypes.some((t) => areTypesCompatible(sourceType, t))) continue;
    matches.push({
      typeName,
      label: labelForTypeName(typeName),
      isBuiltin: !typeName.startsWith('DG Functions/'),
    });
  }
  // Stable sort: Value Output first, then the rest of built-ins, then DG
  // funcs alphabetically.
  matches.sort((a, b) => {
    if (a.typeName === 'Outputs/Value Output') return -1;
    if (b.typeName === 'Outputs/Value Output') return 1;
    if (a.isBuiltin !== b.isBuiltin) return a.isBuiltin ? -1 : 1;
    return a.label.localeCompare(b.label);
  });
  return matches;
}
