import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {cloneConfig, EnumeratorConfig} from './config';

const DISABLED_HINT = 'Use -1 to disable this filter.';

function intInput(label: string, value: number, tooltip?: string): {input: DG.InputBase<number | null>; get: () => number} {
  const input = ui.input.int(label, {value, nullable: false});
  if (tooltip) input.setTooltip(tooltip);
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

function sectionHeader(text: string): HTMLElement {
  const h = ui.h2(text);
  h.style.marginTop = '16px';
  h.style.marginBottom = '6px';
  h.style.borderBottom = '1px solid var(--grey-2)';
  h.style.paddingBottom = '4px';
  return h;
}

export async function openConfigDialog(initial: EnumeratorConfig): Promise<EnumeratorConfig | null> {
  // Start from a clone so cancel/discard leaves `initial` untouched. Fields not exposed in this
  // dialog (file paths, column names, delimiter, output path, depth_first, num_rounds, etc. —
  // all configured in the main UI or unused at runtime) are preserved as-is from `initial`.
  const cfg = cloneConfig(initial);

  const general = {
    keepBBs: boolInput('Keep building blocks in output', cfg.keep_building_blocks_in_final_output,
      'Include the original building blocks (round 0) in the final product list.'),
    maxComponents: intInput('Max # components', cfg.max_num_components,
      'Max number of reactant components a template may have.'),
    maxRoutes: intInput('Max routes per compound', cfg.max_num_routes_per_compound,
      `Cap on number of routes saved per product. ${DISABLED_HINT}`),
    maxCombos: intInput('Max combinations per template', cfg.max_num_combinations_per_template,
      'Per template per round: cap on the number of reactant combinations actually run. If the cartesian product exceeds this, the enumerator runs the first N and stops.'),
  };

  const ps = cfg.products_specs;
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

  const buildSection = (title: string, inputs: {input: DG.InputBase<any>}[]) =>
    ui.div([sectionHeader(title), ui.form(inputs.map((i) => i.input))]);

  const body = ui.div([
    buildSection('General Settings', [general.keepBBs, general.maxComponents, general.maxRoutes, general.maxCombos]),
    buildSection('Product Filters', [
      products.maxHeavy, products.minC, products.maxC, products.maxHetero,
      products.maxN, products.maxS, products.maxO, products.maxMetals, products.maxHal,
      products.maxArom, products.maxUnsat,
      products.allowedAtoms,
      products.rmRadicals, products.rmIsotopes, products.rmCharged,
    ]),
  ]);
  body.style.minWidth = '420px';
  body.style.maxHeight = '70vh';
  body.style.overflowY = 'auto';
  body.style.paddingRight = '8px';

  return new Promise<EnumeratorConfig | null>((resolve) => {
    const dlg = ui.dialog({title: 'Enumerator — full config'});
    dlg.add(body);
    dlg.onOK(() => {
      // Mutate the clone in-place so any field NOT exposed in this dialog (column names, file
      // paths, delimiter, output path, depth_first, num_rounds, etc.) keeps its original value.
      cfg.keep_building_blocks_in_final_output = general.keepBBs.get();
      cfg.max_num_components = general.maxComponents.get();
      cfg.max_num_routes_per_compound = general.maxRoutes.get();
      cfg.max_num_combinations_per_template = general.maxCombos.get();
      Object.assign(cfg.products_specs, {
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
      resolve(cfg);
    });
    dlg.onCancel(() => resolve(null));
    dlg.show({resizable: true, width: 500});
  });
}
