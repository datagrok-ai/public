export const STORAGE_NAME = 'admet_models';
export const KEY = 'selected';
export const TEMPLATES_FOLDER = 'System:AppData/Admetica/templates';
export const DEFAULT_LOWER_VALUE = 0.8;
export const DEFAULT_UPPER_VALUE = 1.0;
export const DEFAULT_TABLE_NAME = 'table';

export interface AdmeticaResponse {
  success: boolean;
  error: string | null;
  result: string | null;
}

export interface ModelProperty {
  name: string;
  description: string;
  units: string;
  direction?: 'Lower is better' | 'Higher is better';
  ranges?: { [key: string]: string };
  weight: number;
  object: { [key: string]: {
    [key: string]: string
  }};
}

export interface ModelColoring {
  type: 'Linear' | 'Conditional';
  min?: number;
  max?: number;
  colors?: string;
}

export interface Model {
  name: string;
  checked: boolean;
  units: string;
  min: number;
  max: number;
  properties: ModelProperty[];
  coloring: ModelColoring;
  line: any;
  weight: number;
}

export interface Template {
  name: string;
  subgroup: Subgroup[];
}

export interface Subgroup {
  name: string;
  expanded: boolean;
  checked: boolean,
  description: string;
  models: Model[];
}

export const TAGS = {
  SECTOR_COLOR: '.sectorColor',
  LOW: '.low',
  HIGH: '.high',
  WEIGHT: '.weight',
  GROUP_NAME: '.group-name',
};

export const ERROR_MESSAGES = {
  MALFORMED: 'Molecule is possibly malformed',
  EMPTY: 'Molecule is empty',
};

export const colorsDictionary: { [key: string]: string } = {
  'Absorption': '#1f77b4',
  'Distribution': '#cf5a0c',
  'Metabolism': '#2ca02c',
  'Excretion': '#d62728',
  'Toxicity': '#9467bd',
};
