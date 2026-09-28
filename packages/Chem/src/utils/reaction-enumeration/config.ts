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
  max_num_aromatic_rings: number;
  max_num_unsaturated_nonaromatic_bonds: number;
  only_these_atoms_allowed: string[];
  remove_radicals: boolean;
  remove_isotope_information: boolean;
  remove_charged_species: boolean;
}

export const AGGREGATIONS = ['sum', 'multiply', 'avg', 'min', 'max'] as const;
export type Aggregation = typeof AGGREGATIONS[number];
/** Column name → its aggregations over the route, one result column each; [] keeps the per-step columns only. */
export type PropagatedColumns = Record<string, Aggregation[]>;

// No file-path fields: the app never auto-loads a file from config, so a persisted name could only
// go stale and mislead.
export interface EnumerationSpecs {
  smarts_col: string;
  reactant_blocking_groups_per_template_column: string;
  bb_smiles_column: string;
  reagent_smiles_column: string;
  reaction_name_col: string;
  template_propagated_columns: PropagatedColumns;
  bb_propagated_columns: PropagatedColumns;
  reagent_propagated_columns: PropagatedColumns;
  depth_first: boolean;
  num_rounds: number;
}

export interface EnumeratorConfig {
  keep_building_blocks_in_final_output: boolean;
  max_num_components: number;
  max_num_routes_per_compound: number;
  max_num_combinations_per_template: number;
  max_num_products_per_step: number;
  products_specs: ProductsSpecs;
  enumeration: EnumerationSpecs;
}

export const DEFAULT_CONFIG: EnumeratorConfig = {
  keep_building_blocks_in_final_output: false,
  max_num_components: -1,
  max_num_routes_per_compound: -1,
  max_num_combinations_per_template: -1,
  max_num_products_per_step: -1,
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
    max_num_aromatic_rings: -1,
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
    template_propagated_columns: {},
    bb_propagated_columns: {},
    reagent_propagated_columns: {},
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

/** Type-checks `partial` against DEFAULT_CONFIG's own shape, so there's no second schema to keep in
 * sync. Without it a wrong-typed value from a hand-edited YAML reaches a numeric comparison at
 * filter time — e.g. a string in a product-filter field makes its `>= 0` check permanently false,
 * silently disabling that filter. */
function validateShape(partial: any, defaults: any, path: string, errors: string[]): void {
  for (const k of Object.keys(defaults)) {
    if (!(k in partial)) continue;
    const expected = defaults[k];
    const actual = partial[k];
    if (Array.isArray(expected)) {
      if (!Array.isArray(actual) || !actual.every((x: unknown) => typeof x === 'string'))
        errors.push(`'${path}${k}' must be a list of strings.`);
    } else if (typeof expected === 'object') {
      // An array or primitive here would otherwise recurse as if it were the section's object,
      // validating nothing, or throw a raw TypeError instead of this function's own message. A null
      // top-level section reads as empty, but a null map inside a section would replace its default.
      const wrongType = actual != null && (typeof actual !== 'object' || Array.isArray(actual));
      if (wrongType || (actual === null && path !== ''))
        errors.push(`'${path}${k}' must be an object.`);
      else
        validateShape(actual ?? {}, expected, `${path}${k}.`, errors);
    } else if (typeof actual !== typeof expected || (typeof expected === 'number' && !Number.isFinite(actual)))
      errors.push(`'${path}${k}' must be a ${typeof expected}.`);
  }
}

/** validateShape checks these maps against their empty defaults, so it never reaches their entries. */
function validatePropagatedColumns(en: {[key: string]: unknown} | undefined, errors: string[]): void {
  for (const key of ['template_propagated_columns', 'bb_propagated_columns', 'reagent_propagated_columns']) {
    const map = en?.[key];
    if (!map || typeof map !== 'object' || Array.isArray(map)) continue;
    for (const [col, aggs] of Object.entries(map)) {
      if (!Array.isArray(aggs) || !aggs.every((a) => (AGGREGATIONS as readonly unknown[]).includes(a))) {
        errors.push(`'enumeration.${key}.${col}' must be a list drawn from ${AGGREGATIONS.join(', ')} ` +
          '(an empty list for none).');
      }
    }
  }
}

export function configFromYaml(text: string): EnumeratorConfig {
  const raw = yaml.load(text);
  // The Array check matters: a top-level YAML list is a truthy 'object', and every key lookup
  // below then misses on it, so the load would silently return pure defaults with no error.
  if (!raw || typeof raw !== 'object' || Array.isArray(raw))
    throw new Error('YAML did not parse to an object.');
  const errors: string[] = [];
  validateShape(raw, DEFAULT_CONFIG, '', errors);
  validatePropagatedColumns((raw as {enumeration?: {[key: string]: unknown}}).enumeration, errors);
  if (errors.length > 0) throw new Error(`Invalid config: ${errors.join('; ')}`);
  return mergeWithDefaults(raw as Partial<EnumeratorConfig>);
}

// Dropped from the schema, but a YAML saved before that still carries them; without this the
// Object.assign below lets them tag along and an old file re-saved leaks them forward.
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
