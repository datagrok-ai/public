/** DG type ↔ socket type & color mappings, plus type-compatibility rules
 *  used by `TypedSocket.isCompatibleWith`. */

import * as DG from 'datagrok-api/dg';

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

// ---- node identity colors come from the platform's categorical palette ----

/** Core's Standard palette (`Color.category20` in d4's color.dart) — the
 *  compile-time fallback when `DG.Color.categoricalPalette` is unreachable.
 *  Order matters: `CAT` below names entries by index. */
const STANDARD_PALETTE = [
  '#1f77b4', '#ffbb78', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
  '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#98df8a', '#ff9896',
  '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
];

/** Hue names → index in the Standard palette. Identity colors are picked by
 *  index at runtime, so a user-customized platform palette carries into Flow. */
export const CAT = {
  blue: 0, orange: 1, green: 2, red: 3, purple: 4, brown: 5, pink: 6, gray: 7,
  olive: 8, cyan: 9, lightGreen: 10, lightRed: 11, lightPurple: 12,
  lightBrown: 13, lightPink: 14, lightGray: 15, lightOlive: 16, lightCyan: 17,
};

let _palette: string[] | null = null;

/** The i-th color of the platform's categorical palette
 *  (`DG.Color.categoricalPalette` — what users see on every categorical
 *  coloring across Datagrok), looping like `getCategoricalColor`. */
export function categoricalColor(i: number): string {
  if (_palette == null) {
    try {
      const raw = DG.Color.categoricalPalette;
      _palette = raw?.length ? raw.map((c) => DG.Color.toHtml(c)) : STANDARD_PALETTE;
    }
    catch (_) {
      _palette = STANDARD_PALETTE;
    }
  }
  return _palette[i % _palette.length];
}

const white = (color: string): {color: string; bgcolor: string} => ({color, bgcolor: '#ffffff'});

/** Role → title-bar color (white body). Looked up by `FuncNode` from
 *  `func.options.role`. */
export const ROLE_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'app': white(categoricalColor(CAT.blue)),
  'panel': white(categoricalColor(CAT.purple)),
  'viewer': white(categoricalColor(CAT.lightCyan)),
  'transform': white(categoricalColor(CAT.cyan)),
  'Transform': white(categoricalColor(CAT.cyan)),
  'filter': white(categoricalColor(CAT.olive)),
  'converter': white(categoricalColor(CAT.orange)),
  'widget': white(categoricalColor(CAT.pink)),
  'cellRenderer': white(categoricalColor(CAT.brown)),
  'semTypeDetector': white(categoricalColor(CAT.lightOlive)),
  'fileViewer': white(categoricalColor(CAT.lightCyan)),
  'fileExporter': white(categoricalColor(CAT.lightCyan)),
  'editor': white(categoricalColor(CAT.lightCyan)),
  'searchProvider': white(categoricalColor(CAT.lightGreen)),
  'tooltip': white(categoricalColor(CAT.lightRed)),
};

export const DEFAULT_NODE_COLOR = categoricalColor(CAT.lightGray);
export const DEFAULT_NODE_BGCOLOR = '#ffffff';

/** How much white is mixed into a node's identity color for the on-canvas
 *  title bar (see `pastelize`). One knob for the whole palette. */
export const TITLE_WHITE_RATIO = 0.6;

/** Soften an identity color by mixing it with white — same hue, much lighter,
 *  so title bars blend with the canvas. The vivid original stays canonical
 *  (`node.color`: minimap, tests, future legends); only the rendered title bar
 *  uses the pastel. Non-`#rrggbb` inputs are returned unchanged. */
export function pastelize(hex: string, whiteRatio: number = TITLE_WHITE_RATIO): string {
  const m = /^#([0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return hex;
  const n = parseInt(m[1], 16);
  const mix = (c: number): number => Math.round(c + (255 - c) * whiteRatio);
  const [r, g, b] = [mix((n >> 16) & 0xff), mix((n >> 8) & 0xff), mix(n & 0xff)];
  return '#' + ((1 << 24) | (r << 16) | (g << 8) | b).toString(16).slice(1);
}

// ---- domain sections (cheminformatics / bioinformatics) ----

/** Packages whose functions are grouped under **Cheminformatics** (small
 *  molecules: structures, descriptors, similarity/substructure, ADMET, docking,
 *  reactions, chemical DBs). Membership is by source package — a domain a
 *  scientist recognizes at a glance, orthogonal to the signature-based task
 *  categories used for everything else. */
export const CHEMINFORMATICS_PACKAGES = new Set<string>([
  'Chem', 'Chembl', 'ChemblApi', 'PubchemApi', 'Chemspace', 'Surechembl',
  'Admetica', 'Docking', 'Retrosynthesis', 'Marvin', 'ChemDrawSketcher',
  'KetcherSketcher', 'HitTriage', 'Datagrokdsmf', 'Curves',
]);

/** Packages whose functions are grouped under **Bioinformatics** (sequences,
 *  peptides, macromolecules, structures, oligos). */
export const BIOINFORMATICS_PACKAGES = new Set<string>([
  'Bio', 'SequenceTranslator', 'Helm', 'Proteomics', 'Bionemo', 'Biologics',
  'OligoBatchCalculator', 'Parabilisseq', 'Sequenceutils', 'BiostructureViewer',
  'PhyloTreeViewer', 'Peptides',
]);

/** The domain section a function belongs to based on its source package, or
 *  `null` for general (task-categorized) functions. Chem wins over Bio for the
 *  (currently empty) intersection. */
export function domainSection(packageName: string | null | undefined): 'Cheminformatics' | 'Bioinformatics' | null {
  if (!packageName) return null;
  if (CHEMINFORMATICS_PACKAGES.has(packageName)) return 'Cheminformatics';
  if (BIOINFORMATICS_PACKAGES.has(packageName)) return 'Bioinformatics';
  return null;
}

/** Whether a function *operates on data it is given* — i.e. takes a dataframe or
 *  column input. The domain sections hold only such operations; a chem/bio
 *  function that merely *produces* a table from scalars (a DB query, fetch, or
 *  generator) is a data source, not an operation, and is left to its signature
 *  category (Data Sources) so "Cheminformatics"/"Bioinformatics" stay about
 *  doing something to the scientist's table — never queries. */
export function isDomainOperation(inputTypes: string[]): boolean {
  return inputTypes.some((t) => t === 'dataframe' || t === 'column' || t === 'column_list');
}

/** The domain section for a function only when it's an operation on data (see
 *  `isDomainOperation`); otherwise `null` (fall back to the task category). */
export function domainCategory(
  packageName: string | null | undefined, inputTypes: string[]): 'Cheminformatics' | 'Bioinformatics' | null {
  const section = domainSection(packageName);
  return section && isDomainOperation(inputTypes) ? section : null;
}

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
  'Data Sources': white(categoricalColor(CAT.orange)),         // orange — bring data in
  'Combine Tables': white(categoricalColor(CAT.purple)),       // purple — join/union
  'Transform Tables': white(categoricalColor(CAT.cyan)),       // cyan — reshape
  'Column Operations': white(categoricalColor(CAT.blue)),      // blue — derive columns
  'Compute Values': white(categoricalColor(CAT.lightGreen)),   // green — scalars
  'Visualize': white(categoricalColor(CAT.lightCyan)),         // cyan — viewers
  'Cheminformatics': white(categoricalColor(CAT.pink)),        // pink — small molecules
  'Bioinformatics': white(categoricalColor(CAT.lightPurple)),  // purple — sequences
  'Other': white(categoricalColor(CAT.lightGray)),             // gray — the rest
};

/** Per-function title-bar colors, keyed by simple function name
 *  (case-insensitive). Checked before role-based coloring, so specific
 *  functions can be visually pinned regardless of their role. Add an entry to
 *  give any function a fixed color. */
export const FUNC_NAME_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'setvar': white(categoricalColor(CAT.red)),      // red — variable assignment
  'getvar': white(categoricalColor(CAT.lightRed)), // light red — variable read
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
  // `list<string>` is just the parametrized spelling of `string_list` — fold it
  // so the socket type (and everything keyed off it) is the same.
  if (dgType === 'list<string>') return 'string_list';
  const mapped = DG_TYPE_MAP[dgType];
  return mapped ? mapped.slotType : dgType;
}

/** A comma-separated string-list input — editable inline like a column-list
 *  (the value is a comma-separated string; the compiler turns it into a JS
 *  array of trimmed, non-empty strings). `list<string>` is the same as
 *  `string_list`; plain `list` (which may hold non-strings) is intentionally
 *  excluded. Matches either a raw DG `propertyType` or a resolved slot type. */
export function isStringListType(dgType: string): boolean {
  return dgType === 'string_list' || dgType === 'list<string>';
}

/** Comma-separated string → JS array literal of trimmed, non-empty strings
 *  (`"a, b ,c"` → `["a", "b", "c"]`, empty → `[]`). Shared by the compilers. */
export function stringListToArrayLiteral(value: unknown): string {
  const items = String(value ?? '').split(',').map((s) => s.trim()).filter((s) => s.length > 0);
  return `[${items.map((s) => JSON.stringify(s)).join(', ')}]`;
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
