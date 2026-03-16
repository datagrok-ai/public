import {LiteGraph} from 'litegraph.js';
import * as DG from 'datagrok-api/dg';
import {createFuncNodeClass} from './func-node';
import {registerInputNodes} from './input-nodes';
import {registerOutputNodes} from './output-nodes';
import {registerUtilityNodes} from './utility-nodes';
import {registerComparisonNodes} from './comparison-nodes';
import {registerBreakpointNode} from './breakpoint-node';
import {getRole, getTags, getPackageName} from '../utils/dart-proxy-utils';
import {registerSlotColors} from '../types/type-map';

export interface FuncInfo {
  func: DG.Func;
  name: string;
  role: string | null;
  tags: string[];
  packageName: string;
  nodeTypeName: string;
}

/** Tags that cause a function to be excluded from the node catalog.
 * Add tag strings here to filter out unwanted functions. */
export const EXCLUDED_TAGS: string[] = [
  // Add tags to exclude, e.g.: 'internal', 'deprecated'
];

/** Roles that cause a function to be excluded from the node catalog.
 * Add role strings here to filter out unwanted functions. */
export const EXCLUDED_ROLES: string[] = [
  DG.FUNC_TYPES.APP, 'aiSearchProvider', 'antibodyNumbering', 'appTreeBrowser', 'canonicalizer',
  DG.FUNC_TYPES.CELL_RENDERER, 'dashboard', DG.FUNC_TYPES.FILE_VIEWER, DG.FUNC_TYPES.FILE_IMPORTER,
  DG.FUNC_TYPES.FILE_EXPORTER, DG.FUNC_TYPES.FILTER, DG.FUNC_TYPES.FOLDER_VIEWER, DG.FUNC_TYPES.MOLECULE_SKETCHER,
  'notationProviderConstructor', 'notationRefiner', 'packageSettingsEditor', 'searchProvider', 'semTypeDetector',
  'semValueExtractor', 'valueEditor',
  // Add roles to exclude, e.g.: 'semTypeDetector', 'cellRenderer'
];

let registeredFuncs: FuncInfo[] = [];
let builtinsRegistered = false;
let funcsRegistered = false;

/** Clear all default LiteGraph node types (math, basic, etc.) keeping only ours */
function clearDefaultNodeTypes(): void {
  const registered = LiteGraph.registered_node_types;
  for (const key of Object.keys(registered)) {
    // Keep our custom categories
    if (key.startsWith('Inputs/') || key.startsWith('Outputs/') ||
        key.startsWith('Utilities/') || key.startsWith('Constants/') ||
        key.startsWith('Comparisons/') || key.startsWith('DG Functions/'))
      continue;
    delete registered[key];
  }
}

/** Check if a function should be excluded based on its tags and role */
function shouldExcludeFunc(role: string | null, tags: string[]): boolean {
  if (role && EXCLUDED_ROLES.length > 0 && EXCLUDED_ROLES.includes(role))
    return true;
  if (tags.length > 0 && EXCLUDED_TAGS.length > 0) {
    for (const tag of tags) {
      if (EXCLUDED_TAGS.includes(tag))
        return true;
    }
  }
  return false;
}

/** Register all built-in nodes (inputs, outputs, utilities). No-op on subsequent calls. */
export function registerBuiltinNodes(): void {
  if (builtinsRegistered) return;
  builtinsRegistered = true;

  // Clear LiteGraph defaults first (math, basic, etc.)
  clearDefaultNodeTypes();

  registerInputNodes();
  registerOutputNodes();
  registerUtilityNodes();
  registerComparisonNodes();
  registerBreakpointNode();
}

/** Register all DG.Func as LiteGraph node types and return metadata. No-op on subsequent calls. */
export function registerAllFunctions(): FuncInfo[] {
  if (funcsRegistered) return registeredFuncs;
  funcsRegistered = true;

  registerSlotColors();

  const allFuncs = DG.Func.find({});
  registeredFuncs = [];
  const seenNames = new Set<string>();

  for (const func of allFuncs) {
    try {
      // Skip functions with no inputs and no outputs (not useful in a chain)
      if (func.inputs.length === 0 && func.outputs.length === 0) continue;

      const role = getRole(func);
      const tags = getTags(func);

      // Skip functions with excluded tags or roles
      if (shouldExcludeFunc(role, tags)) continue;

      // Use safe getter to avoid Dart proxy crash on func.package
      const pkgName = getPackageName(func);
      const category = role || 'Uncategorized';

      let typeName = `DG Functions/${category}/${func.name}`;
      if (seenNames.has(typeName))
        typeName = `DG Functions/${category}/${pkgName}:${func.name}`;
      if (seenNames.has(typeName)) continue;
      seenNames.add(typeName);

      const nodeClass = createFuncNodeClass(func);
      LiteGraph.registerNodeType(typeName, nodeClass);

      registeredFuncs.push({
        func,
        name: func.name,
        role,
        tags,
        packageName: pkgName,
        nodeTypeName: typeName,
      });
    } catch (e) {
      // Silently skip functions that fail (Dart proxy issues, etc.)
    }
  }

  return registeredFuncs;
}

/** Get the cached list of registered functions */
export function getRegisteredFuncs(): FuncInfo[] {
  return registeredFuncs;
}
