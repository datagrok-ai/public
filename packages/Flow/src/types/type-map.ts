/** DG type ↔ socket type & color mappings, plus type-compatibility rules
 *  used by `TypedSocket.isCompatibleWith`. */

import * as DG from 'datagrok-api/dg';

/** Core's Standard palette — the fallback when `DG.Color.categoricalPalette` is unreachable.
 *  Order matters: `CAT` names entries by index. */
const STANDARD_PALETTE = [
  '#1f77b4', '#ffbb78', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
  '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#98df8a', '#ff9896',
  '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
];

/** Hue names → index in the Standard palette. */
export const CAT = {
  blue: 0, orange: 1, green: 2, red: 3, purple: 4, brown: 5, pink: 6, gray: 7,
  olive: 8, cyan: 9, lightGreen: 10, lightRed: 11, lightPurple: 12,
  lightBrown: 13, lightPink: 14, lightGray: 15, lightOlive: 16, lightCyan: 17,
};

let _palette: string[] | null = null;

/** The i-th color of the platform's categorical palette, looping. */
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

/** Column-data socket types wear core's `Color.typeColors` so a socket's letter+color
 *  pair matches the Column Manager. */
export const DG_TYPE_MAP: Record<string, {slotType: string; color: string}> = {
  'dataframe': {slotType: 'dataframe', color: '#E67E22'},
  'column': {slotType: 'column', color: '#3498DB'},
  'column_list': {slotType: 'column_list', color: '#5DADE2'},
  'string': {slotType: 'string', color: categoricalColor(CAT.orange)},
  'int': {slotType: 'int', color: categoricalColor(CAT.green)},
  'double': {slotType: 'double', color: categoricalColor(CAT.pink)},
  'bool': {slotType: 'bool', color: categoricalColor(CAT.blue)},
  'datetime': {slotType: 'datetime', color: categoricalColor(CAT.brown)},
  'bigint': {slotType: 'bigint', color: categoricalColor(CAT.red)},
  'qnum': {slotType: 'qnum', color: categoricalColor(CAT.purple)},
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
  // Execution-ordering ports: gray, deliberately isolated from every other type.
  'order': {slotType: 'order', color: '#9E9E9E'},
};

/** Role → title-bar color (white body). */
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

/** White ratio mixed into a node's identity color for the title bar — one knob for the whole palette. */
export const TITLE_WHITE_RATIO = 0.6;

/** Mix an identity color with white. The vivid original stays canonical (`node.color`) —
 *  only the rendered title bar uses the pastel. Non-`#rrggbb` inputs return unchanged. */
export function pastelize(hex: string, whiteRatio: number = TITLE_WHITE_RATIO): string {
  const m = /^#([0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return hex;
  const n = parseInt(m[1], 16);
  const mix = (c: number): number => Math.round(c + (255 - c) * whiteRatio);
  const [r, g, b] = [mix((n >> 16) & 0xff), mix((n >> 8) & 0xff), mix(n & 0xff)];
  return '#' + ((1 << 24) | (r << 16) | (g << 8) | b).toString(16).slice(1);
}

/** Packages whose functions are grouped under Cheminformatics. */
export const CHEMINFORMATICS_PACKAGES = new Set<string>([
  'Chem', 'Chembl', 'ChemblApi', 'PubchemApi', 'Chemspace', 'Surechembl',
  'Admetica', 'Docking', 'Retrosynthesis', 'Marvin', 'ChemDrawSketcher',
  'KetcherSketcher', 'HitTriage', 'Datagrokdsmf', 'Curves',
]);

/** Packages whose functions are grouped under Bioinformatics. */
export const BIOINFORMATICS_PACKAGES = new Set<string>([
  'Bio', 'SequenceTranslator', 'Helm', 'Proteomics', 'Bionemo', 'Biologics',
  'OligoBatchCalculator', 'Parabilisseq', 'Sequenceutils', 'BiostructureViewer',
  'PhyloTreeViewer', 'Peptides',
]);

/** The domain section by source package, or null; Chem wins over Bio for the intersection. */
export function domainSection(packageName: string | null | undefined): 'Cheminformatics' | 'Bioinformatics' | null {
  if (!packageName) return null;
  if (CHEMINFORMATICS_PACKAGES.has(packageName)) return 'Cheminformatics';
  if (BIOINFORMATICS_PACKAGES.has(packageName)) return 'Bioinformatics';
  return null;
}

/** Whether a function operates on data it is given. The domain sections hold only such
 *  operations — a chem/bio source (table from scalars) stays with its signature category. */
export function isDomainOperation(inputTypes: string[]): boolean {
  return inputTypes.some((t) => t === 'dataframe' || t === 'column' || t === 'column_list');
}

/** The domain section only for operations on data, else null. */
export function domainCategory(
  packageName: string | null | undefined, inputTypes: string[]): 'Cheminformatics' | 'Bioinformatics' | null {
  const section = domainSection(packageName);
  return section && isDomainOperation(inputTypes) ? section : null;
}

const VIS_TYPES = ['viewer', 'view', 'widget', 'graphics'];
const SCALAR_TYPES = ['string', 'int', 'double', 'bool', 'datetime', 'num', 'bigint', 'qnum'];
const COL_TYPES = ['column', 'column_list'];

/** Bucket a function by its input/output signature; pure — shared by the browser grouping and node coloring. */
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

/** Title-bar color per task category — the fallback that keeps role-less functions from all being gray. */
export const CATEGORY_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'Data Sources': white(categoricalColor(CAT.orange)),
  'Combine Tables': white(categoricalColor(CAT.purple)),
  'Transform Tables': white(categoricalColor(CAT.cyan)),
  'Column Operations': white(categoricalColor(CAT.blue)),
  'Compute Values': white(categoricalColor(CAT.lightGreen)),
  'Visualize': white(categoricalColor(CAT.lightCyan)),
  'Cheminformatics': white(categoricalColor(CAT.pink)),
  'Bioinformatics': white(categoricalColor(CAT.lightPurple)),
  'Other': white(categoricalColor(CAT.lightGray)),
};

/** Per-function title-bar colors by lower-cased simple name, checked before role-based coloring. */
export const FUNC_NAME_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'setvar': white(categoricalColor(CAT.red)),
  'getvar': white(categoricalColor(CAT.lightRed)),
};

/** Symmetric compat map: an output of type K can connect to an input of any listed type,
 *  and vice-versa; `'*'` is a wildcard. */
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
  // Exec ports connect only to each other — checked before the dynamic/object wildcards.
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
  // `list<string>` is the parametrized spelling of `string_list` — fold it.
  if (dgType === 'list<string>') return 'string_list';
  const mapped = DG_TYPE_MAP[dgType];
  return mapped ? mapped.slotType : dgType;
}

/** Scalar types whose input rows are hidden on the node card by default — edited in the panel. */
const PRIMITIVE_SLOT_TYPES = new Set(['string', 'int', 'double', 'num', 'bool', 'datetime', 'bigint', 'qnum']);

export function isPrimitiveSlotType(dgType: string): boolean {
  return PRIMITIVE_SLOT_TYPES.has(dgTypeToSlotType(dgType));
}

/** A comma-separated string-list input; plain `list` (which may hold non-strings) is
 *  intentionally excluded. */
export function isStringListType(dgType: string): boolean {
  return dgType === 'string_list' || dgType === 'list<string>';
}

/** Comma-separated string → JS array literal of trimmed, non-empty strings. */
export function stringListToArrayLiteral(value: unknown): string {
  const items = String(value ?? '').split(',').map((s) => s.trim()).filter((s) => s.length > 0);
  return `[${items.map((s) => JSON.stringify(s)).join(', ')}]`;
}

export function getSlotColor(slotType: string): string {
  const mapped = DG_TYPE_MAP[slotType];
  return mapped ? mapped.color : '#95A5A6';
}

/** Overrides for slot letters that shouldn't be the plain first character. */
const SLOT_LETTERS: Record<string, string> = {
  'dataframe': 't', // "table" — how users say it
  'dynamic': '?',   // wildcard — a letter would suggest a concrete type
  'object': '?',    // wildcard
};

/** Single-letter type abbreviation in a socket chip — mirrors the Column Manager. */
export function getSlotLetter(dgType: string): string {
  const slot = dgTypeToSlotType(dgType);
  const letter = SLOT_LETTERS[slot] ?? slot.trim().charAt(0).toLowerCase();
  return letter === '' ? '?' : letter;
}

export function getNodeColors(
  role: string | null, funcName?: string, category?: string,
): {color: string; bgcolor: string} {
  if (funcName) {
    const override = FUNC_NAME_COLORS[funcName.toLowerCase()];
    if (override) return override;
  }
  if (role && ROLE_COLORS[role]) return ROLE_COLORS[role];
  if (category && CATEGORY_COLORS[category]) return CATEGORY_COLORS[category];
  return {color: DEFAULT_NODE_COLOR, bgcolor: DEFAULT_NODE_BGCOLOR};
}
