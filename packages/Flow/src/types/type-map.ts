/** DG type ↔ socket type & color mappings, plus type-compatibility rules
 *  used by `TypedSocket.isCompatibleWith`. */

export const DG_TYPE_MAP: Record<string, {slotType: string; color: string}> = {
  'dataframe': {slotType: 'dataframe', color: '#E67E22'},
  'column': {slotType: 'column', color: '#3498DB'},
  'column_list': {slotType: 'column_list', color: '#5DADE2'},
  'string': {slotType: 'string', color: '#2ECC71'},
  'int': {slotType: 'int', color: '#1ABC9C'},
  'double': {slotType: 'double', color: '#00BCD4'},
  'bool': {slotType: 'bool', color: '#E74C3C'},
  'datetime': {slotType: 'datetime', color: '#9B59B6'},
  'string_list': {slotType: 'string_list', color: '#8BC34A'},
  'object': {slotType: 'object', color: '#95A5A6'},
  'dynamic': {slotType: 'dynamic', color: '#7E8C8D'},
  'map': {slotType: 'map', color: '#FFC107'},
  'funccall': {slotType: 'funccall', color: '#E91E63'},
  'list': {slotType: 'list', color: '#AB47BC'},
  'file': {slotType: 'file', color: '#78909C'},
  'byte_array': {slotType: 'byte_array', color: '#607D8B'},
  'bitset': {slotType: 'bitset', color: '#FF7043'},
  'num': {slotType: 'num', color: '#26C6DA'},
  'viewer': {slotType: 'viewer', color: '#42A5F5'},
  'graphics': {slotType: 'graphics', color: '#66BB6A'},
  'blob': {slotType: 'byte_array', color: '#607D8B'},
  'view': {slotType: 'view', color: '#5C6BC0'},
  // Execution-ordering ports (control flow, not data). Gray, and deliberately
  // isolated from every other type (see areTypesCompatible).
  'order': {slotType: 'order', color: '#9E9E9E'},
};

/** Role → title-bar color (white body). Looked up by `FuncNode` from
 *  `func.options.role`. */
export const ROLE_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'app': {color: '#7986CB', bgcolor: '#ffffff'},
  'panel': {color: '#9575CD', bgcolor: '#ffffff'},
  'viewer': {color: '#64B5F6', bgcolor: '#ffffff'},
  'transform': {color: '#4DB6AC', bgcolor: '#ffffff'},
  'Transform': {color: '#4DB6AC', bgcolor: '#ffffff'},
  'filter': {color: '#FFD54F', bgcolor: '#ffffff'},
  'converter': {color: '#FFB74D', bgcolor: '#ffffff'},
  'widget': {color: '#F06292', bgcolor: '#ffffff'},
  'cellRenderer': {color: '#A1887F', bgcolor: '#ffffff'},
  'semTypeDetector': {color: '#DCE775', bgcolor: '#ffffff'},
  'fileViewer': {color: '#4DD0E1', bgcolor: '#ffffff'},
  'fileExporter': {color: '#4DD0E1', bgcolor: '#ffffff'},
  'editor': {color: '#4FC3F7', bgcolor: '#ffffff'},
  'searchProvider': {color: '#AED581', bgcolor: '#ffffff'},
  'tooltip': {color: '#FF8A65', bgcolor: '#ffffff'},
};

export const DEFAULT_NODE_COLOR = '#BDBDBD';
export const DEFAULT_NODE_BGCOLOR = '#ffffff';

// ---- categorize a function by what it does (shared by the browser + coloring) ----

const VIS_TYPES = ['viewer', 'view', 'widget', 'graphics'];
const SCALAR_TYPES = ['string', 'int', 'double', 'bool', 'datetime', 'num', 'bigint', 'qnum'];
const COL_TYPES = ['column', 'column_list'];

/** Bucket a function by its input/output signature (and viewer role). Pure —
 *  operates on the lists of DG property-type strings, so it's shared by the
 *  function browser's grouping AND the node title-bar coloring. */
export function categorizeBySignature(ins: string[], outs: string[], role: string | null): string {
  const has = (arr: string[], set: string[]): boolean => arr.some((t) => set.includes(t));
  const dfIn = ins.filter((t) => t === 'dataframe').length;
  const outDf = outs.includes('dataframe');
  const noOut = outs.length === 0;
  const roleHasViewer = !!role && role.split(',').some((r) => r.trim() === 'viewer');

  if (has(outs, VIS_TYPES) || roleHasViewer) return 'Visualize';
  if (dfIn >= 2) return 'Combine Tables';
  if (outDf && dfIn === 0) return 'Data Sources';
  if ((outDf && dfIn === 1) || (noOut && dfIn >= 1)) return 'Transform Tables';
  if (has(outs, COL_TYPES)) return 'Column Operations';
  if (has(outs, SCALAR_TYPES)) return 'Compute Values';
  return 'Other';
}

/** Title-bar color per task category — so a function with no role (most of them:
 *  JoinTables, AddNewColumn, chem properties, …) still reads its job from color
 *  instead of all being gray. */
export const CATEGORY_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'Data Sources': {color: '#FF8A65', bgcolor: '#ffffff'},      // orange — bring data in
  'Combine Tables': {color: '#BA68C8', bgcolor: '#ffffff'},    // purple — join/union
  'Transform Tables': {color: '#4DB6AC', bgcolor: '#ffffff'},  // teal — reshape
  'Column Operations': {color: '#5C9DED', bgcolor: '#ffffff'}, // blue — derive columns
  'Compute Values': {color: '#9CCC65', bgcolor: '#ffffff'},    // green — scalars
  'Visualize': {color: '#4DD0E1', bgcolor: '#ffffff'},         // cyan — viewers
  'Other': {color: '#90A4AE', bgcolor: '#ffffff'},             // blue-gray — the rest
};

/** Per-function title-bar colors, keyed by simple function name
 *  (case-insensitive). Checked before role-based coloring, so specific
 *  functions can be visually pinned regardless of their role. Add an entry to
 *  give any function a fixed color. */
export const FUNC_NAME_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'setvar': {color: '#EF5350', bgcolor: '#ffffff'}, // red — variable assignment
  'getvar': {color: '#EF9A9A', bgcolor: '#ffffff'}, // light red — variable read
};

/** Symmetric compat map: an output of type K can connect to an input of any
 *  type listed in `COMPATIBLE_TYPES[K]`, and vice-versa. `'*'` is a wildcard. */
const COMPATIBLE_TYPES: Record<string, string[]> = {
  'double': ['int', 'num'],
  'tableview': ['view'],
  'num': ['int', 'double'],
  'list': ['string_list'],
  'string_list': ['list'],
  'dynamic': ['*'],
  'object': ['*'],
};

export function areTypesCompatible(outputType: string, inputType: string): boolean {
  if (outputType === inputType) return true;
  // Execution-ordering ports connect ONLY to each other — checked before the
  // dynamic/object wildcards so a data port can never plug into an exec port.
  if (outputType === 'order' || inputType === 'order') return false;
  if (outputType === 'dynamic' || inputType === 'dynamic') return true;
  if (outputType === 'object' || inputType === 'object') return true;
  const inCompat = COMPATIBLE_TYPES[inputType];
  if (inCompat && (inCompat.includes('*') || inCompat.includes(outputType))) return true;
  const outCompat = COMPATIBLE_TYPES[outputType];
  if (outCompat && (outCompat.includes('*') || outCompat.includes(inputType))) return true;
  return false;
}

export function dgTypeToSlotType(dgType: string): string {
  const mapped = DG_TYPE_MAP[dgType];
  return mapped ? mapped.slotType : dgType;
}

export function getSlotColor(slotType: string): string {
  const mapped = DG_TYPE_MAP[slotType];
  return mapped ? mapped.color : '#95A5A6';
}

export function getNodeColors(
  role: string | null, funcName?: string, category?: string,
): {color: string; bgcolor: string} {
  if (funcName) {
    const override = FUNC_NAME_COLORS[funcName.toLowerCase()];
    if (override) return override;
  }
  if (role && ROLE_COLORS[role]) return ROLE_COLORS[role];
  // Fall back to the task category — gives role-less functions (the gray
  // majority) a color that reflects what they do.
  if (category && CATEGORY_COLORS[category]) return CATEGORY_COLORS[category];
  return {color: DEFAULT_NODE_COLOR, bgcolor: DEFAULT_NODE_BGCOLOR};
}
