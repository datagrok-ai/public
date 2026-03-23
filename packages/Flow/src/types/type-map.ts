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
  'blob': {slotType: 'byte_array', color: '#607D8B'},
  'view': {slotType: 'view', color: '#5C6BC0'},
};

/** Role → node color mapping (light theme: colored title bar, white body) */
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
