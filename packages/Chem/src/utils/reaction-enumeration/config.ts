import * as yaml from 'js-yaml';

export interface ProductsSpecs {
  exclusion_smarts_products_file: string | null;
  exclusion_smarts_products_file_smarts_col: string;
  max_num_heavy_atoms: number;
  min_num_carbon_atoms: number;
  max_num_carbon_atoms: number;
  max_num_hetero_atoms: number;
  max_num_nitrogen: number;
  max_num_sulfur: number;
  max_num_oxygen: number;
  max_num_metals: number;
  max_num_halogens: number;
  max_num_aromatic_atoms: number;
  max_num_unsaturated_nonaromatic_bonds: number;
  only_these_atoms_allowed: string[];
  remove_radicals: boolean;
  remove_isotope_information: boolean;
  remove_charged_species: boolean;
}

export interface EnumerationSpecs {
  template_file: string;
  smarts_col: string;
  reactant_blocking_groups_per_template_column: string;
  bb_file: string;
  bb_smiles_column: string;
  reagent_file: string | null;
  reagent_smiles_column: string;
  output_file: string;
  reaction_name_col: string;
  delimiter: string;
  depth_first: boolean;
  num_rounds: number;
}

export interface EnumeratorConfig {
  keep_building_blocks_in_final_output: boolean;
  max_num_components: number;
  max_num_routes_per_compound: number;
  max_num_combinations_per_template: number;
  products_specs: ProductsSpecs;
  enumeration: EnumerationSpecs;
}

export const DEFAULT_CONFIG: EnumeratorConfig = {
  keep_building_blocks_in_final_output: false,
  max_num_components: 4,
  max_num_routes_per_compound: -1,
  max_num_combinations_per_template: 50,
  products_specs: {
    exclusion_smarts_products_file: 'ex_smarts.csv',
    exclusion_smarts_products_file_smarts_col: 'SMARTS',
    max_num_heavy_atoms: -1,
    min_num_carbon_atoms: 10,
    max_num_carbon_atoms: 30,
    max_num_hetero_atoms: 10,
    max_num_nitrogen: -1,
    max_num_sulfur: -1,
    max_num_oxygen: -1,
    max_num_metals: 0,
    max_num_halogens: -1,
    max_num_aromatic_atoms: -1,
    max_num_unsaturated_nonaromatic_bonds: 5,
    only_these_atoms_allowed: ['C', 'H', 'O', 'N', 'S', 'P'],
    remove_radicals: true,
    remove_isotope_information: true,
    remove_charged_species: true,
  },
  enumeration: {
    template_file: 'reactions.csv',
    smarts_col: 'reaction_smarts',
    reactant_blocking_groups_per_template_column: 'blocking_fg',
    bb_file: 'bb.csv',
    bb_smiles_column: 'SMILES',
    reagent_file: null,
    reagent_smiles_column: 'SMILES',
    output_file: 'enumeration_output.csv',
    reaction_name_col: 'reaction_name',
    delimiter: ',',
    depth_first: true,
    num_rounds: 2,
  },
};

export function cloneConfig(c: EnumeratorConfig): EnumeratorConfig {
  return JSON.parse(JSON.stringify(c));
}

export function configToYaml(c: EnumeratorConfig): string {
  const dump: any = JSON.parse(JSON.stringify(c));
  if (dump.enumeration?.reagent_file == null)
    dump.enumeration.reagent_file = 'None';
  return yaml.dump(dump, {lineWidth: 120, noRefs: true, sortKeys: false});
}

export function configFromYaml(text: string): EnumeratorConfig {
  const raw = yaml.load(text);
  if (!raw || typeof raw !== 'object')
    throw new Error('YAML did not parse to an object.');
  return mergeWithDefaults(raw as Partial<EnumeratorConfig>);
}

export function mergeWithDefaults(partial: Partial<EnumeratorConfig> | any): EnumeratorConfig {
  const out = cloneConfig(DEFAULT_CONFIG);
  if (!partial) return out;
  for (const k of Object.keys(out) as (keyof EnumeratorConfig)[]) {
    const v = (partial as any)[k];
    if (v == null) continue;
    if (k === 'products_specs' || k === 'enumeration')
      Object.assign(out[k], v);
    else
      (out as any)[k] = v;
  }
  const rf = out.enumeration.reagent_file;
  if (typeof rf === 'string' && (rf === 'None' || rf === 'none' || rf.trim() === ''))
    out.enumeration.reagent_file = null;
  return out;
}
