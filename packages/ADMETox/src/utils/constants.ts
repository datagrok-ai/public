export let STORAGE_NAME = 'admet_models';
export let KEY = 'selected';
export let TEMPLATES_FOLDER = 'System:AppData/Admetox/templates';

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
  units: string;
  min: number | string;
  max: number | string;
  properties: ModelProperty[];
  coloring: ModelColoring;
}

export interface Template {
  name: string;
  subgroup: Subgroup[];
}

export interface Subgroup {
  name: string;
  models: Model[];
}