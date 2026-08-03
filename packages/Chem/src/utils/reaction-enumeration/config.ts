import * as yaml from 'js-yaml';

export interface ProductsSpecs {
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

// No file-path fields here (template/BB/reagent/output file names, CSV delimiter) — the app never
// auto-loads a file from config, the user always picks each one manually, so persisting a file name
// would only ever go stale and mislead (it can't be kept in sync with what's actually loaded).
export interface EnumerationSpecs {
  smarts_col: string;
  reactant_blocking_groups_per_template_column: string;
  bb_smiles_column: string;
  reagent_smiles_column: string;
  reaction_name_col: string;
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
    smarts_col: 'reaction_smarts',
    reactant_blocking_groups_per_template_column: 'blocking_fg',
    bb_smiles_column: 'SMILES',
    reagent_smiles_column: 'SMILES',
    reaction_name_col: 'reaction_name',
    depth_first: true,
    num_rounds: 2,
  },
};

export function cloneConfig(c: EnumeratorConfig): EnumeratorConfig {
  return JSON.parse(JSON.stringify(c));
}

export function configToYaml(c: EnumeratorConfig): string {
  return yaml.dump(cloneConfig(c), {lineWidth: 120, noRefs: true, sortKeys: false});
}

// Walks `partial` against `defaults`' own shape — DEFAULT_CONFIG already has every field with its
// correct type, so this needs no separately-maintained schema. Catches the case a hand-edited or
// corrupted YAML silently carries a wrong-typed value all the way into a numeric comparison at
// filter time (e.g. a string in a product-filter field makes that filter's `>= 0` check
// permanently false, disabling it with no indication anything's wrong).
function validateShape(partial: any, defaults: any, path: string, errors: string[]): void {
  for (const k of Object.keys(defaults)) {
    if (!(k in partial)) continue;
    const expected = defaults[k];
    const actual = partial[k];
    if (Array.isArray(expected)) {
      if (!Array.isArray(actual) || !actual.every((x: unknown) => typeof x === 'string'))
        errors.push(`'${path}${k}' must be a list of strings.`);
    } else if (typeof expected === 'object') {
      // Same reason as configFromYaml's own top-level guard: an array or primitive here would
      // otherwise recurse as if it were the section's object (silently validating nothing, since
      // `k in partial` never matches a named key against an array/primitive) or throw a raw
      // `TypeError: Cannot use 'in' operator...` instead of this function's own clean message.
      if (actual != null && (typeof actual !== 'object' || Array.isArray(actual)))
        errors.push(`'${path}${k}' must be an object.`);
      else
        validateShape(actual ?? {}, expected, `${path}${k}.`, errors);
    } else if (typeof actual !== typeof expected || (typeof expected === 'number' && !Number.isFinite(actual))) {
      errors.push(`'${path}${k}' must be a ${typeof expected}.`);
    }
  }
}

export function configFromYaml(text: string): EnumeratorConfig {
  const raw = yaml.load(text);
  // typeof [] === 'object' and arrays are always truthy, so a top-level YAML list would otherwise
  // sail past this guard — every key lookup below (`k in partial`, `partial[k]`) then misses on an
  // array, so validateShape finds nothing to flag and mergeWithDefaults finds nothing to merge,
  // silently returning pure defaults with no error at all.
  if (!raw || typeof raw !== 'object' || Array.isArray(raw))
    throw new Error('YAML did not parse to an object.');
  const errors: string[] = [];
  validateShape(raw, DEFAULT_CONFIG, '', errors);
  if (errors.length > 0) throw new Error(`Invalid config: ${errors.join('; ')}`);
  return mergeWithDefaults(raw as Partial<EnumeratorConfig>);
}

// File-name fields (template_file, bb_file, reagent_file, output_file, delimiter,
// exclusion_smarts_products_file) were removed from the schema — a YAML saved before that change
// may still carry them. Object.assign below would otherwise let them tag along as inert extra
// properties, so an old file loaded then re-saved keeps leaking them forward.
const LEGACY_ENUMERATION_KEYS = ['template_file', 'bb_file', 'reagent_file', 'output_file', 'delimiter'];
const LEGACY_PRODUCTS_SPECS_KEYS = ['exclusion_smarts_products_file'];

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
  for (const k of LEGACY_ENUMERATION_KEYS) delete (out.enumeration as any)[k];
  for (const k of LEGACY_PRODUCTS_SPECS_KEYS) delete (out.products_specs as any)[k];
  return out;
}
