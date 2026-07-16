import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {EnumeratorConfig} from './config';

const DISABLED_HINT = 'Leave blank to disable this filter.';

// The platform's +/- steppers can't increment from blank and have no floor on decrement.
// Override both: blank "+" -> 1, "0 -" -> blank, otherwise a normal ±1.
export function fixNullableIntStepper(input: DG.InputBase<number | null>): void {
  const apply = (next: number | null): void => {
    input.value = next;
    const el = input.input as HTMLInputElement | undefined;
    if (el) el.value = next == null ? '' : String(next);
    try {input.fireChanged();} catch (e) {_package.logger.debug(`fireChanged unavailable: ${e}`);}
  };
  const override = (selector: string, step: 1 | -1): void => {
    input.root.querySelector(selector)?.addEventListener('click', (e) => {
      e.stopImmediatePropagation();
      e.preventDefault();
      const cur = input.value;
      if (step > 0) apply(cur == null ? 1 : cur + 1);
      else apply(cur != null && cur > 0 ? cur - 1 : null);
    }, true);
  };
  override('.ui-input-plus', 1);
  override('.ui-input-minus', -1);
}

// -1 is each field's "no cap" sentinel — shown as blank rather than a literal -1.
function intInput(
  label: string, value: number, tooltip?: string,
): {input: DG.InputBase<number | null>; get: () => number} {
  const input = ui.input.int(label, {value: value < 0 ? undefined : value, nullable: true, showPlusMinus: true});
  if (tooltip) input.setTooltip(tooltip);
  fixNullableIntStepper(input);
  return {input, get: () => (input.value == null ? -1 : input.value)};
}

function boolInput(label: string, value: boolean, tooltip?: string) {
  const input = ui.input.bool(label, {value});
  if (tooltip) input.setTooltip(tooltip);
  return {input, get: () => !!input.value};
}

function csvInput(label: string, items: string[], tooltip?: string) {
  const input = ui.input.string(label, {value: items.join(', ')});
  if (tooltip) input.setTooltip(tooltip);
  return {
    input,
    get: () => (input.value ?? '').split(',').map((s) => s.trim()).filter((s) => s.length > 0),
  };
}

// `syncToConfig` writes the current field values into the given config object.
export function buildCombinationLimitFields(initial: EnumeratorConfig): {
  root: HTMLElement;
  inputs: DG.InputBase<unknown>[];
  syncToConfig: (target: EnumeratorConfig) => void;
} {
  const keepBBs = boolInput('Keep building blocks in output', initial.keep_building_blocks_in_final_output,
    'Include the original building blocks (round 0) in the final product list.');
  const maxCombos = intInput('Max combinations per template', initial.max_num_combinations_per_template,
    'Per template per round: cap on the number of reactant combinations actually run. If the ' +
    'cartesian product exceeds this, the enumerator runs the first N and stops. Leave blank for no cap.');

  const root = ui.form([keepBBs.input, maxCombos.input]);

  const syncToConfig = (target: EnumeratorConfig): void => {
    target.keep_building_blocks_in_final_output = keepBBs.get();
    target.max_num_combinations_per_template = maxCombos.get();
  };

  return {root, inputs: [keepBBs.input, maxCombos.input], syncToConfig};
}

// Product-level filters: atom counts, charge/radical/isotope rejection.
export function buildProductFilterFields(initial: EnumeratorConfig): {
  root: HTMLElement;
  inputs: DG.InputBase<unknown>[];
  syncToConfig: (target: EnumeratorConfig) => void;
} {
  const ps = initial.products_specs;
  const products = {
    maxHeavy: intInput('Max heavy atoms', ps.max_num_heavy_atoms, DISABLED_HINT),
    minC: intInput('Min carbon atoms', ps.min_num_carbon_atoms, DISABLED_HINT),
    maxC: intInput('Max carbon atoms', ps.max_num_carbon_atoms, DISABLED_HINT),
    maxHetero: intInput('Max hetero atoms', ps.max_num_hetero_atoms, DISABLED_HINT),
    maxN: intInput('Max nitrogen atoms', ps.max_num_nitrogen, DISABLED_HINT),
    maxS: intInput('Max sulfur atoms', ps.max_num_sulfur, DISABLED_HINT),
    maxO: intInput('Max oxygen atoms', ps.max_num_oxygen, DISABLED_HINT),
    maxMetals: intInput('Max metals', ps.max_num_metals, DISABLED_HINT),
    maxHal: intInput('Max halogens', ps.max_num_halogens, DISABLED_HINT),
    maxArom: intInput('Max aromatic atoms', ps.max_num_aromatic_atoms, DISABLED_HINT),
    maxUnsat: intInput('Max unsaturated non-aromatic bonds', ps.max_num_unsaturated_nonaromatic_bonds, DISABLED_HINT),
    allowedAtoms: csvInput('Only these atoms allowed', ps.only_these_atoms_allowed,
      'Comma-separated element symbols. Leave empty to allow any atom.'),
    rmRadicals: boolInput('Reject radicals', ps.remove_radicals,
      'Reject products containing radical atoms.'),
    rmIsotopes: boolInput('Strip isotope information', ps.remove_isotope_information,
      'Strip isotope info from products and re-canonicalize.'),
    rmCharged: boolInput('Reject charged species', ps.remove_charged_species,
      'Reject products with non-zero formal charges.'),
  };

  const inputs = Object.values(products).map((p) => p.input);
  const root = ui.form(inputs);

  const syncToConfig = (target: EnumeratorConfig): void => {
    Object.assign(target.products_specs, {
      max_num_heavy_atoms: products.maxHeavy.get(),
      min_num_carbon_atoms: products.minC.get(),
      max_num_carbon_atoms: products.maxC.get(),
      max_num_hetero_atoms: products.maxHetero.get(),
      max_num_nitrogen: products.maxN.get(),
      max_num_sulfur: products.maxS.get(),
      max_num_oxygen: products.maxO.get(),
      max_num_metals: products.maxMetals.get(),
      max_num_halogens: products.maxHal.get(),
      max_num_aromatic_atoms: products.maxArom.get(),
      max_num_unsaturated_nonaromatic_bonds: products.maxUnsat.get(),
      only_these_atoms_allowed: products.allowedAtoms.get(),
      remove_radicals: products.rmRadicals.get(),
      remove_isotope_information: products.rmIsotopes.get(),
      remove_charged_species: products.rmCharged.get(),
    });
  };

  return {root, inputs, syncToConfig};
}
