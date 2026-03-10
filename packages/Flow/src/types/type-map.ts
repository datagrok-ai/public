import {LiteGraph} from 'litegraph.js';

/** Datagrok TYPE → LiteGraph slot type string mapping with colors */
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
};

/** Role → node color mapping (light theme: white body, colored title bar) */
export const ROLE_COLORS: Record<string, {color: string; bgcolor: string}> = {
  'app': {color: '#5C6BC0', bgcolor: '#ffffff'},
  'panel': {color: '#7E57C2', bgcolor: '#ffffff'},
  'viewer': {color: '#42A5F5', bgcolor: '#ffffff'},
  'transform': {color: '#26A69A', bgcolor: '#ffffff'},
  'Transform': {color: '#26A69A', bgcolor: '#ffffff'},
  'filter': {color: '#FFCA28', bgcolor: '#ffffff'},
  'converter': {color: '#FFA726', bgcolor: '#ffffff'},
  'widget': {color: '#EC407A', bgcolor: '#ffffff'},
  'cellRenderer': {color: '#8D6E63', bgcolor: '#ffffff'},
  'semTypeDetector': {color: '#C0CA33', bgcolor: '#ffffff'},
  'fileViewer': {color: '#26C6DA', bgcolor: '#ffffff'},
  'fileExporter': {color: '#26C6DA', bgcolor: '#ffffff'},
  'editor': {color: '#29B6F6', bgcolor: '#ffffff'},
  'searchProvider': {color: '#9CCC65', bgcolor: '#ffffff'},
  'tooltip': {color: '#FF7043', bgcolor: '#ffffff'},
};

export const DEFAULT_NODE_COLOR = '#BDBDBD';
export const DEFAULT_NODE_BGCOLOR = '#ffffff';

/** Type compatibility: which types can connect to which */
const COMPATIBLE_TYPES: Record<string, string[]> = {
  'double': ['int', 'num'],
  'num': ['int', 'double'],
  'list': ['string_list'],
  'string_list': ['list'],
  'dynamic': ['*'], // special: connects to everything
  'object': ['*'],
};

export function areTypesCompatible(outputType: string, inputType: string): boolean {
  if (outputType === inputType) return true;
  if (outputType === 'dynamic' || inputType === 'dynamic') return true;
  if (outputType === 'object' || inputType === 'object') return true;
  const compat = COMPATIBLE_TYPES[inputType];
  if (compat && (compat.includes('*') || compat.includes(outputType))) return true;
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

export function getNodeColors(role: string | null): {color: string; bgcolor: string} {
  if (role && ROLE_COLORS[role]) return ROLE_COLORS[role];
  return {color: DEFAULT_NODE_COLOR, bgcolor: DEFAULT_NODE_BGCOLOR};
}

/** Register all DG types with LiteGraph's link_type_colors for colored connections */
export function registerSlotColors(): void {
  for (const [, {slotType, color}] of Object.entries(DG_TYPE_MAP))
    LGraphCanvas.link_type_colors[slotType] = color;
}

// We need to import LGraphCanvas separately for the static property
import {LGraphCanvas} from 'litegraph.js';
