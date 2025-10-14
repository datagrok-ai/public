import {BaseType, Selection} from 'd3-selection';
// Base types for nodes and reactions
export interface BaseNode {
  bigg_id: string
  name: string | null
  data?: number | null | string | any
  data_string?: string | null
  [key: string]: any
}

export interface BaseReaction {
  bigg_id: string
  name: string | null
  data?: any
  data_string?: string | null
  [key: string]: any
}

// COBRA Model types
export interface CobraGene {
  bigg_id: string
  name: string | null
  [key: string]: any
}

export interface CobraMetabolite extends BaseNode {
  formula?: string
  compartment?: string
  // Node properties
  x?: number
  y?: number
  node_type?: string
  connected_segments?: string[]
  label_x?: number
  label_y?: number
  index?: number
  count?: number
  is_primary?: boolean,
  coefficient: number
}

export interface CobraReaction extends BaseReaction {
  metabolites: { [key: string]: number } // metabolite_id -> coefficient
  gene_reaction_rule?: string
  genes?: CobraGene[]
  reversibility: boolean
  segments?: { [key: string]: MapSegment }
  // Reaction properties
  label_x?: number
  label_y?: number
  objective_coefficient?: number
  id: string
  lower_bound?: number
  upper_bound?: number
  subsystem?: string
}

export type ReactionBounds = {
  upper_bound: number
  lower_bound: number
}

export interface CobraModelData {
  id?: string
  name?: string
  description?: string
  notes?: any // TODO: define the structure of notes
  reactions: CobraReaction[]
  metabolites: CobraMetabolite[]
  genes: CobraGene[]
}

// Map types
export interface MapNode extends BaseNode {
  x: number
  y: number
  node_type: string
  node_id?: IDType
  node_is_primary?: boolean
  connected_segments: {segment_id: string,
            reaction_id: string}[]
  label_x: number
  label_y: number,
  [key: string]: any
}

export interface MapSegment {
  id?: string
  from_node_id: string
  to_node_id: string
  b1?: { x: number; y: number } | null
  b2?: { x: number; y: number } | null
  data?: any | null
  reversibility?: boolean,
  from_node_coefficient?: number | null
  to_node_coefficient?: number | null
  segment_id?: string | number | null
  reverse_flux?: boolean,
  unconnected_segment_with_arrow?: boolean,
  reaction_id?: string
}

export type IDType = string;

export interface MapReaction extends BaseReaction {
  gene_reaction_rule: string
  genes: CobraGene[]
  metabolites: {name?: string, bigg_id: string, coefficient: number, index?: number, count?: number, is_primary?: boolean}[]
  segments: {[key: IDType]: MapSegment}
  reversibility: boolean
  reverse_flux?: boolean,
  label_x: number
  label_y: number,
  reaction_id: string,
  gene_string?: { bigg_id: string | null, name: string | null, text: string }[] | null,
  [key: string]: any
}

export interface MapData {
  reactions: { [key: string]: MapReaction }
  nodes: { [key: string]: MapNode }
  text_labels: { [key: string]: TextLabel }
  canvas: {
    x: number
    y: number
    width: number
    height: number
  }
}

export type EscherMapDataType = [{map_name: string, map_id: string, map_description: string}, MapData];

export interface TextLabel {
  text: string
  x: number
  y: number
}

// Data types
export type DataValue = string | number | null
export type DataArray = DataValue[]
export type GeneValues = { [key: string]: DataArray }
export type CompareStyle = 'diff' | 'fold' | 'log2_fold'

// Settings types
export interface SettingsType {
  //todo: add types here
  reaction_data: { [key: string]: (number | string | null) } | null
  metabolite_data: any | null
  gene_data: any | null
  reaction_styles: string[]
  metabolite_styles: string[]
  identifiers_on_map: 'bigg_id' | 'name'
  reaction_compare_style: CompareStyle
  metabolite_compare_style: CompareStyle
  and_method_in_gene_reaction_rule: 'mean' | 'min'
  highlight_missing: boolean
  enable_search: boolean
  enable_editing: boolean
  disabled_buttons: string[]
  zoom_to_element?: {
    type: 'reaction' | 'node'
    id: string
  }
  starting_reaction?: string
  never_ask_before_quit: boolean
  show_gene_reaction_rules: boolean
  hide_secondary_metabolites: boolean
  hide_all_labels: boolean
  allow_building_duplicate_reactions: boolean
  bezier_curve_options?: any
  tooltip_component?: any
  enable_tooltips?: string[]
  enable_keys_with_tooltip?: boolean
  enable_tooltips_global?: boolean
  reaction_scale: [number, number]
  metabolite_scale: [number, number]
  primary_metabolite_radius: number
  secondary_metabolite_radius: number
  scroll_behavior?: 'zoom' | 'pan'
  fill_screen?: boolean
  samplingFunction?: (mp: CobraModelData) => SamplingFunctionResult | Promise<SamplingFunctionResult>
  marker_radius: number
  gene_font_size: number
  canvas_size_and_loc?: {
    x: number
    y: number
    width: number
    height: number
  }
  saveAction?: (() => void) | null
  loadAction?: (() => void) | null
  pathFindingDisabled?: boolean,
  runFBA?: () => Promise<void>
}

export type SamplingFunctionResult = {
  upper_bound: number;
  lower_bound: number;
  data: Map<string, number[]>;
  cancled?: boolean;
}

// DataStyles types
export type Styles = string[]

// Builder types
export interface BuilderData {
  map_data?: MapData
  model_data?: CobraModelData
  settings?: Partial<SettingsType>
}

export type D3Selection<T extends Element = Element> = Selection<T, unknown, BaseType, any>
export type D3SVGSelection = Selection<SVGSVGElement, unknown, BaseType, any>
export type D3DivSelection = Selection<HTMLDivElement, unknown, BaseType, any>
export type D3ElementSelection = Selection<Element, unknown, BaseType, any>

export interface D3DragEvent {
  x: number
  y: number
  dx: number
  dy: number
  sourceEvent: MouseEvent | TouchEvent
}
export type Coord = {x: number, y: number}

export type ReactionSamplingDistribution = {
  lower_bound: number
  upper_bound: number
  data: Map<IDType, number[]> // reaction_id -> [values of distribution counts]
}

export interface TooltipComponentProps {
  display: boolean;
  disableTooltips: () => void;
  biggId: string;
  loc: Coord;
  data: string;
  type: string;
  name: string;
  [key: string]: any;
}
