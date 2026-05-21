/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {cloneConfig, configFromYaml, configToYaml, DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {openConfigDialog} from './config-form';
import {getRdKitModule} from '../chem-common-rdkit';
import {enumerate, EnumerationProgress, OutputRow, TemplateInput} from './enumerate';

const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
const BUNDLED_BBS = 'enumerations/bb.csv';
const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

// Drop-in replacement for the standalone package's rdkit-helper. By the time the app function
// runs, Chem's `init` handler has already initialized the RDKit module.
async function getRdKit() {
  return getRdKitModule();
}

// Sniff string columns and set semType so the grid renders reactions and molecules. We sample a
// handful of non-empty values per column: presence of `>>` wins as ChemicalReaction; otherwise if
// every sampled value looks like SMILES/SMARTS we mark the column as Molecule.
function detectChemSemTypes(df: DG.DataFrame): void {
  for (const col of df.columns.toList()) {
    if (col.type !== DG.COLUMN_TYPE.STRING) continue;
    if (col.semType) continue;
    const samples: string[] = [];
    const n = Math.min(col.length, 50);
    for (let i = 0; i < n && samples.length < 5; i++) {
      const v = col.get(i);
      if (v == null) continue;
      const s = String(v).trim();
      if (s.length > 0) samples.push(s);
    }
    if (samples.length === 0) continue;
    if (samples.some((s) => s.includes('>>')))
      col.semType = 'ChemicalReaction';
    df.meta.detectSemanticTypes();
  }
}

async function loadBundledCsv(name: string): Promise<DG.DataFrame | null> {
  try {
    const text = await _package.files.readAsText(name);
    const df = DG.DataFrame.fromCsv(text);
    df.name = name.replace(/\.csv$/i, '');
    detectChemSemTypes(df);
    await df.meta.detectSemanticTypes();
    return df;
  } catch (e) {
    console.warn(`Could not load bundled file ${name}: ${e}`);
    return null;
  }
}

function pickFile(accept: string): Promise<File | null> {
  return new Promise((resolve) => {
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = accept;
    input.style.display = 'none';
    document.body.appendChild(input);
    input.onchange = () => {
      const f = input.files?.[0] ?? null;
      document.body.removeChild(input);
      resolve(f);
    };
    input.click();
  });
}

function downloadText(text: string, filename: string, mime = 'text/plain'): void {
  const blob = new Blob([text], {type: mime});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function getStringColumn(df: DG.DataFrame, name: string): string[] {
  const col = df.col(name);
  if (!col) throw new Error(`Column "${name}" not found in "${df.name}". Available: ${df.columns.names().join(', ')}`);
  const out: string[] = new Array(col.length);
  for (let i = 0; i < col.length; i++) {
    const v = col.get(i);
    out[i] = v == null ? '' : String(v);
  }
  return out;
}

function buildResultDataFrame(rows: OutputRow[], name = 'Enumeration result'): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('product', rows.map((r) => r.product)),
    DG.Column.fromStrings('route', rows.map((r) => r.route)),
    DG.Column.fromStrings('template', rows.map((r) => r.template)),
    DG.Column.fromStrings('reaction_name', rows.map((r) => r.reaction_name)),
    DG.Column.fromInt32Array('round', new Int32Array(rows.map((r) => r.round))),
    DG.Column.fromInt32Array('n_routes', new Int32Array(rows.map((r) => r.n_routes))),
  ]);
  df.name = name;
  df.col('product')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('route')!.semType = 'ChemicalReaction';
  df.col('template')!.semType = 'ChemicalReaction';
  return df;
}

interface BuiltInputs {
  templates: TemplateInput[];
  buildingBlocks: string[];
  exclusionSmarts: string[];
  reagents: string[];
}

const isStringCol = (c: DG.Column) => c.type === DG.COLUMN_TYPE.STRING;

// Build a column input bound to a table. If the preferred column name exists and matches the
// filter, it's pre-selected; otherwise the first matching column is used (Datagrok default).
function makeColInput(
  label: string, table: DG.DataFrame | null, preferredName: string,
  filter: (c: DG.Column) => boolean, tooltip: string, nullable: boolean,
): DG.InputBase<DG.Column | null> {
  const opts: any = {filter};
  if (table) {
    opts.table = table;
    opts.nullable = nullable;
    const c = table.col(preferredName);
    if (c && filter(c)) opts.value = c;
  }
  const inp = ui.input.column(label, opts);
  inp.setTooltip(tooltip);
  return inp;
}

// When the parent table changes, swap the column input's table in-place (preserves DOM identity
// so form layout stays valid). Re-selects the preferred column name in the new table if it exists.
function bindColToTable(
  colInput: DG.InputBase<DG.Column | null>, tableInput: DG.InputBase<DG.DataFrame | null>,
  getPreferredName: () => string, filter: (c: DG.Column) => boolean,
): void {
  tableInput.onChanged.subscribe(() => {
    const t = tableInput.value;
    if (!t) return;
    ui.input.setColumnInputTable(colInput, t, filter);
    const wanted = getPreferredName();
    if (wanted) {
      const c = t.col(wanted);
      if (c && filter(c)) colInput.value = c;
    }
  });
}

function buildInputs(
  config: EnumeratorConfig, tDf: DG.DataFrame, bDf: DG.DataFrame,
  xDf: DG.DataFrame | null, rDf: DG.DataFrame | null,
): BuiltInputs {
  const smartsList = getStringColumn(tDf, config.enumeration.smarts_col);
  const blockingCol = config.enumeration.reactant_blocking_groups_per_template_column;
  const rxnNameCol = config.enumeration.reaction_name_col;
  const blockingListRaw = tDf.col(blockingCol) ? getStringColumn(tDf, blockingCol) : null;
  const rxnNameList = tDf.col(rxnNameCol) ? getStringColumn(tDf, rxnNameCol) : null;

  const templates: TemplateInput[] = [];
  for (let i = 0; i < smartsList.length; i++) {
    const smarts = (smartsList[i] ?? '').trim();
    if (!smarts) continue;
    const blockingRaw = blockingListRaw ? blockingListRaw[i] : '';
    const blockingSmartsList = blockingRaw ?
      blockingRaw.split(/[;|]/).map((s) => s.trim()).filter((s) => s.length > 0) : [];
    templates.push({smarts, blockingSmartsList, reactionName: rxnNameList?.[i] ?? ''});
  }

  const buildingBlocks = getStringColumn(bDf, config.enumeration.bb_smiles_column)
    .filter((s) => s.trim().length > 0);

  let exclusionSmarts: string[] = [];
  if (xDf) {
    const excCol = config.products_specs.exclusion_smarts_products_file_smarts_col;
    if (xDf.col(excCol))
      exclusionSmarts = getStringColumn(xDf, excCol).filter((s) => s.trim().length > 0);
  }

  let reagents: string[] = [];
  if (rDf) {
    const rCol = config.enumeration.reagent_smiles_column;
    if (rDf.col(rCol))
      reagents = getStringColumn(rDf, rCol).filter((s) => s.trim().length > 0);
  }

  return {templates, buildingBlocks, exclusionSmarts, reagents};
}

export async function buildEnumeratorView(): Promise<DG.ViewBase> {
  const view = DG.View.create();
  view.name = 'Chemical Enumerator';
  view.box = true;

  let config: EnumeratorConfig = cloneConfig(DEFAULT_CONFIG);

  const templatesDf = await loadBundledCsv(BUNDLED_TEMPLATES);
  const bbsDf = await loadBundledCsv(BUNDLED_BBS);
  const exclusionDf = await loadBundledCsv(BUNDLED_EXCLUSION);

  // ---- DATA inputs ----
  const templatesInput = ui.input.table('Reaction templates file', {
    value: templatesDf ?? undefined, nullable: false,
    tooltipText: 'Table with reaction SMARTS templates. Pick an open workspace table or upload a CSV via the file icon.',
  });
  const smartsColInput = makeColInput('Reaction SMARTS column', templatesDf, config.enumeration.smarts_col, isStringCol,
    'Column in the reaction templates file that contains the reaction SMARTS strings.', false);

  const blockingColInput = makeColInput('Blocking SMARTS column (Optional)', templatesDf,
    config.enumeration.reactant_blocking_groups_per_template_column, isStringCol,
    'Optional column whose values are SMARTS patterns (separated by ";" or "|"). A building block ' +
    'matching any blocking SMARTS for a given template is excluded from that template only.', true);

  const rxnNameColInput = makeColInput('Reaction name column (Optional)', templatesDf,
    config.enumeration.reaction_name_col, isStringCol,
    'Optional column with a friendly name for each reaction template. Surfaces in the output ' +
    '"reaction_name" column.', true);

  const bbsInput = ui.input.table('Building blocks file', {
    value: bbsDf ?? undefined, nullable: false,
    tooltipText: 'Table with the building-block library (SMILES). Pick a workspace table or upload a CSV.',
  });
  const bbColInput = makeColInput('Building blocks SMILES column', bbsDf, config.enumeration.bb_smiles_column, isStringCol,
    'Column in the building blocks file that contains SMILES.', false);

  const reagentsInput = ui.input.table('Reagents file (Optional)', {
    value: undefined, nullable: true,
    tooltipText: 'Optional table of reagent SMILES. When set, switches to reagents mode: every ' +
      'step uses exactly one building block (or product of an earlier round) and fills every ' +
      'remaining slot with reagents from this file — produces derivatives of each BB across rounds.',
  });
  const reagentsColInput = makeColInput('Reagent SMILES column', null,
    config.enumeration.reagent_smiles_column, isStringCol,
    'Column in the reagents file that contains the reagent SMILES.', true);

  const exclusionInput = ui.input.table('Exclusion substructures file (Optional)', {
    value: exclusionDf ?? undefined, nullable: true,
    tooltipText: 'Optional table of SMARTS patterns. Any product matching one of these is rejected.',
  });
  const exclusionColInput = makeColInput('Exclusion substructures column', exclusionDf,
    config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol,
    'Column in the exclusion substructures file that contains the SMARTS strings.', true);

  // Re-bind column inputs whenever the parent table changes. Preserves the column input identity
  // (so form layout stays valid) and re-selects the preferred column name in the new table.
  bindColToTable(smartsColInput, templatesInput, () => config.enumeration.smarts_col, isStringCol);
  bindColToTable(blockingColInput, templatesInput,
    () => config.enumeration.reactant_blocking_groups_per_template_column, isStringCol);
  bindColToTable(rxnNameColInput, templatesInput, () => config.enumeration.reaction_name_col, isStringCol);
  bindColToTable(bbColInput, bbsInput, () => config.enumeration.bb_smiles_column, isStringCol);
  bindColToTable(reagentsColInput, reagentsInput,
    () => config.enumeration.reagent_smiles_column, isStringCol);
  bindColToTable(exclusionColInput, exclusionInput,
    () => config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol);

  // ---- CONFIG inputs ----
  const numRoundsInput = ui.input.int('Number of rounds', {value: config.enumeration.num_rounds, min: 1});
  numRoundsInput.setTooltip(
    'Number of consecutive enumeration rounds. Round 1 reacts BBs only; round 2 takes round-1 ' +
    'products and (in depth-first mode) reacts each one with original BBs. Increase for deeper libraries.');

  const depthFirstInput = ui.input.bool('Depth first', {value: config.enumeration.depth_first});
  depthFirstInput.setTooltip(
    'When checked, each round-r > 1 step must combine EXACTLY ONE round-(r-1) product with original ' +
    'BBs (linear chain extension, no merging two complex products). Off (breadth-first) allows any ' +
    'combination from rounds 0..r-1 — typically explodes the search space and produces convergent routes.');

  const minCInput = ui.input.int('Min carbon atoms', {value: config.products_specs.min_num_carbon_atoms});
  minCInput.setTooltip('Minimum number of carbon atoms a product must have to pass. Set -1 to disable.');

  const maxCInput = ui.input.int('Max carbon atoms', {value: config.products_specs.max_num_carbon_atoms});
  maxCInput.setTooltip('Maximum number of carbon atoms a product may have. Set -1 to disable.');

  const maxCombosInput = ui.input.int('Max combinations / template',
    {value: config.max_num_combinations_per_template});
  maxCombosInput.setTooltip(
    'Per template, per round: cap on how many reactant combinations are actually run. If the ' +
    'cartesian product of matching BBs across slots exceeds this, the enumerator runs the first ' +
    'N combos and stops. Set -1 to disable.');

  const maxComponentsInput = ui.input.int('Max # components', {value: config.max_num_components, min: 1});
  maxComponentsInput.setTooltip(
    'Templates whose reaction SMARTS have more reactant slots than this are skipped. Default 4 ' +
    'covers all templates in the bundled set.');

  // True while we're pushing config → inputs (syncConfigToQuickInputs). Each setAndFire fires
  // onChanged, which triggers refreshTooltip → syncQuickInputsToConfig — and that read-back
  // happens MID-loop while only some inputs have been updated, overwriting config with the stale
  // values of the inputs we haven't reached yet. The flag short-circuits the read-back during the
  // sync; refreshTooltip just reads the live `config` instead.
  let pushingConfigToInputs = false;

  // ---- Info icons ----
  // Two distinct tooltips:
  //   - appInfoIcon (next to the view title) describes WHAT this app is and HOW it works.
  //   - configInfoIcon (next to the "Enumeration options" header) shows the FULL current
  //     config — every YAML field — as a structured, scannable card (not raw YAML).
  // Both are bound to live factories so the rendered HTML reflects the current state every time
  // the user hovers over the icon; no proactive refresh is needed.
  const mkIcon = (): HTMLElement => {
    const i = ui.iconFA('info-circle', () => {});
    i.style.marginLeft = '8px';
    i.style.color = 'var(--blue-2)';
    i.style.cursor = 'help';
    return i;
  };
  const appInfoIcon = mkIcon();
  const configInfoIcon = mkIcon();

  function buildAppHelp(): HTMLElement {
    const card = ui.div();
    card.style.fontSize = '12px';
    card.style.maxWidth = '520px';
    card.style.lineHeight = '1.5';
    card.style.padding = '4px 2px';
    card.innerHTML = `
      <div style="font-weight: bold; font-size: 13px; margin-bottom: 4px;">Chemical library enumeration</div>
      <p style="margin: 0 0 6px 0;">Generate a product library from reaction SMARTS templates and a set of starting materials. Each round can take products from the previous round and grow them further.</p>
      <div style="font-weight: bold; margin-top: 8px; margin-bottom: 2px;">Inputs</div>
      <ul style="margin: 0 0 6px 16px; padding: 0;">
        <li><b>Reaction templates file</b> — reaction SMARTS (LHS&gt;&gt;RHS). Optional blocking SMARTS exclude unwanted matches per template; optional reaction names surface in the output.</li>
        <li><b>Building blocks file</b> — SMILES of starting materials.</li>
        <li><b>Reagents file</b> (optional) — switches to <i>reagents mode</i>: every step uses exactly one BB or earlier-round product, with reagents in the remaining slots. Yields derivatives of each BB across rounds.</li>
        <li><b>Exclusion substructures</b> (optional) — SMARTS patterns; any product matching one is rejected.</li>
      </ul>
      <div style="font-weight: bold; margin-top: 8px; margin-bottom: 2px;">Enumeration modes</div>
      <ul style="margin: 0 0 6px 16px; padding: 0;">
        <li><b>Depth-first</b> — round 1 combines BBs; round R extends each round-(R-1) product with original BBs (linear chains, no merging of two complex products).</li>
        <li><b>Breadth-first</b> — each round may combine any products from earlier rounds with BBs (convergent routes possible).</li>
        <li><b>Reagents</b> — active whenever a reagents file is selected; overrides depth/breadth-first.</li>
      </ul>
      <div style="font-weight: bold; margin-top: 8px; margin-bottom: 2px;">Tips</div>
      <ul style="margin: 0 0 0 16px; padding: 0;">
        <li>Select rows on the right-pane grids and click <i>Subset by selection</i> to enumerate only a subset.</li>
        <li>Open the <i>Preview</i> tab to sample products at a reduced budget before kicking off a full run.</li>
        <li>The <i>i</i> icon next to <i>Enumeration options</i> shows the full current config.</li>
      </ul>`;
    return card;
  }

  function buildConfigCard(): HTMLElement {
    if (!pushingConfigToInputs) syncQuickInputsToConfig();
    const en = config.enumeration;
    const ps = config.products_specs;
    const reagentsActive = reagentsInput.value != null;
    const modeLabel = reagentsActive ? 'Reagents (linear chain extension with reagents in other slots)' :
      en.depth_first ? 'Depth-first (linear chain extension)' : 'Breadth-first (convergent allowed)';
    const fmtNum = (n: number, hint = 'unlimited') => n < 0 ? hint : String(n);
    const yn = (b: boolean) => b ? 'Yes' : 'No';

    const card = ui.div();
    card.style.fontSize = '12px';
    card.style.maxHeight = '500px';
    card.style.maxWidth = '440px';
    card.style.overflow = 'auto';
    card.style.padding = '4px 2px';
    card.style.lineHeight = '1.5';

    const sectionTitle = (text: string) => {
      const h = ui.div();
      h.textContent = text;
      h.style.fontWeight = 'bold';
      h.style.marginTop = '8px';
      h.style.marginBottom = '2px';
      h.style.paddingBottom = '2px';
      h.style.borderBottom = '1px solid var(--grey-3)';
      return h;
    };
    const row = (label: string, value: string) => {
      const r = ui.divH([], {style: {justifyContent: 'space-between', padding: '1px 0', gap: '12px'}});
      const l = ui.div(); l.textContent = label; l.style.color = 'var(--grey-6)';
      const v = ui.div(); v.textContent = value; v.style.color = 'var(--text-color)';
      v.style.textAlign = 'right';
      r.appendChild(l); r.appendChild(v);
      return r;
    };

    card.appendChild(sectionTitle('Enumeration'));
    card.appendChild(row('Rounds', String(en.num_rounds)));
    card.appendChild(row('Mode', modeLabel));
    card.appendChild(row('Max components', String(config.max_num_components)));
    card.appendChild(row('Max combinations / template', fmtNum(config.max_num_combinations_per_template)));
    card.appendChild(row('Max routes / compound', fmtNum(config.max_num_routes_per_compound)));
    card.appendChild(row('Keep BBs in final output', yn(config.keep_building_blocks_in_final_output)));

    card.appendChild(sectionTitle('Columns'));
    card.appendChild(row('Reaction SMARTS', en.smarts_col));
    card.appendChild(row('Blocking SMARTS', en.reactant_blocking_groups_per_template_column));
    card.appendChild(row('Reaction name', en.reaction_name_col));
    card.appendChild(row('Building block SMILES', en.bb_smiles_column));
    card.appendChild(row('Reagent SMILES', en.reagent_smiles_column));
    card.appendChild(row('Exclusion SMARTS', ps.exclusion_smarts_products_file_smarts_col));

    card.appendChild(sectionTitle('Product filters'));
    card.appendChild(row('Carbons (min – max)',
      `${fmtNum(ps.min_num_carbon_atoms, 'any')} – ${fmtNum(ps.max_num_carbon_atoms, 'any')}`));
    card.appendChild(row('Heavy atoms (max)', fmtNum(ps.max_num_heavy_atoms, 'any')));
    card.appendChild(row('Hetero atoms (max)', fmtNum(ps.max_num_hetero_atoms, 'any')));
    card.appendChild(row('Nitrogen (max)', fmtNum(ps.max_num_nitrogen, 'any')));
    card.appendChild(row('Sulfur (max)', fmtNum(ps.max_num_sulfur, 'any')));
    card.appendChild(row('Oxygen (max)', fmtNum(ps.max_num_oxygen, 'any')));
    card.appendChild(row('Metals (max)', fmtNum(ps.max_num_metals, 'any')));
    card.appendChild(row('Halogens (max)', fmtNum(ps.max_num_halogens, 'any')));
    card.appendChild(row('Aromatic atoms (max)', fmtNum(ps.max_num_aromatic_atoms, 'any')));
    card.appendChild(row('Unsat. non-arom. bonds (max)', fmtNum(ps.max_num_unsaturated_nonaromatic_bonds, 'any')));
    card.appendChild(row('Allowed atoms',
      ps.only_these_atoms_allowed.length ? ps.only_these_atoms_allowed.join(', ') : 'any'));
    card.appendChild(row('Reject radicals', yn(ps.remove_radicals)));
    card.appendChild(row('Strip isotopes', yn(ps.remove_isotope_information)));
    card.appendChild(row('Reject charged species', yn(ps.remove_charged_species)));

    return card;
  }

  ui.tooltip.bind(appInfoIcon, () => buildAppHelp());
  ui.tooltip.bind(configInfoIcon, () => buildConfigCard());

  // Kept for back-compat with existing call sites — the tooltip is bound to live factories so
  // hovering always shows up-to-date state; no rebind is needed. We still sync inputs → config
  // so callers that rely on that side effect keep working.
  const refreshTooltip = (): void => {
    if (!pushingConfigToInputs) syncQuickInputsToConfig();
  };

  const syncQuickInputsToConfig = () => {
    config.enumeration.num_rounds = numRoundsInput.value ?? config.enumeration.num_rounds;
    config.enumeration.depth_first = !!depthFirstInput.value;
    config.products_specs.min_num_carbon_atoms = minCInput.value ?? -1;
    config.products_specs.max_num_carbon_atoms = maxCInput.value ?? -1;
    config.max_num_combinations_per_template = maxCombosInput.value ?? config.max_num_combinations_per_template;
    config.max_num_components = maxComponentsInput.value ?? config.max_num_components;
    // Column inputs hold a Column object; persist its name in the YAML config. If the input has no
    // selection, keep the previous config value so YAML round-trip stays stable.
    config.enumeration.smarts_col = smartsColInput.value?.name ?? config.enumeration.smarts_col;
    config.enumeration.reactant_blocking_groups_per_template_column =
      blockingColInput.value?.name ?? config.enumeration.reactant_blocking_groups_per_template_column;
    config.enumeration.reaction_name_col =
      rxnNameColInput.value?.name ?? config.enumeration.reaction_name_col;
    config.enumeration.bb_smiles_column = bbColInput.value?.name ?? config.enumeration.bb_smiles_column;
    config.enumeration.reagent_smiles_column =
      reagentsColInput.value?.name ?? config.enumeration.reagent_smiles_column;
    config.products_specs.exclusion_smarts_products_file_smarts_col =
      exclusionColInput.value?.name ?? config.products_specs.exclusion_smarts_products_file_smarts_col;
  };

  // Setting `input.value = X` updates the model but, on int/float/string inputs, does NOT always
  // refresh the rendered <input type="number"> text — the Dart-side widget skips the re-render
  // when the value comes via the API rather than user typing. The fix is to also push the new
  // value into the underlying HTMLInputElement so its visible text matches the model.
  const setAndFire = <T>(input: DG.InputBase<T>, v: T) => {
    input.value = v;
    try {
      const el = input.input as HTMLInputElement | undefined;
      if (el?.tagName === 'INPUT' && el.type !== 'checkbox') {
        const desired = v == null ? '' : String(v);
        if (el.value !== desired) el.value = desired;
      }
    } catch {/* ignore — non-textual inputs (column/table/bool/etc.) */}
    try {input.fireChanged();} catch {/* ignore — older API versions */}
  };

  const syncConfigToQuickInputs = () => {
    pushingConfigToInputs = true;
    try {
      setAndFire(numRoundsInput, config.enumeration.num_rounds);
      setAndFire(depthFirstInput, config.enumeration.depth_first);
      setAndFire(minCInput, config.products_specs.min_num_carbon_atoms);
      setAndFire(maxCInput, config.products_specs.max_num_carbon_atoms);
      setAndFire(maxCombosInput, config.max_num_combinations_per_template);
      setAndFire(maxComponentsInput, config.max_num_components);
      // For column inputs, look up the column object by name on the currently-selected table.
      const tDf = templatesInput.value;
      if (tDf) {
        const sc = tDf.col(config.enumeration.smarts_col);
        if (sc) setAndFire(smartsColInput, sc);
        const bc = tDf.col(config.enumeration.reactant_blocking_groups_per_template_column);
        if (bc) setAndFire(blockingColInput, bc);
        const rc = tDf.col(config.enumeration.reaction_name_col);
        if (rc) setAndFire(rxnNameColInput, rc);
      }
      const bDf = bbsInput.value;
      if (bDf) {
        const c = bDf.col(config.enumeration.bb_smiles_column);
        if (c) setAndFire(bbColInput, c);
      }
      const rDf = reagentsInput.value;
      if (rDf) {
        const c = rDf.col(config.enumeration.reagent_smiles_column);
        if (c) setAndFire(reagentsColInput, c);
      }
      const xDf = exclusionInput.value;
      if (xDf) {
        const c = xDf.col(config.products_specs.exclusion_smarts_products_file_smarts_col);
        if (c) setAndFire(exclusionColInput, c);
      }
    } finally {
      pushingConfigToInputs = false;
    }
  };

  // ---- Validation ----
  const validationDiv = ui.divText('', {style: {color: 'var(--red-3)', fontSize: '12px', minHeight: '16px', flex: '0 0 auto'}});

  function validate(): string | null {
    syncQuickInputsToConfig();
    const tDf = templatesInput.value;
    if (!tDf) return 'Select a reaction templates file.';
    if (!smartsColInput.value) return 'Select a reaction SMARTS column.';

    const bDf = bbsInput.value;
    if (!bDf) return 'Select a building blocks file.';
    if (!bbColInput.value) return 'Select a SMILES column.';

    const rDf = reagentsInput.value;
    if (rDf && !reagentsColInput.value)
      return 'Select a reagent SMILES column or clear the reagents file.';

    const xDf = exclusionInput.value;
    if (xDf && !exclusionColInput.value)
      return 'Select an exclusion substructures column or clear the exclusion file.';

    const rounds = numRoundsInput.value ?? 0;
    if (rounds < 1) return 'Number of rounds must be at least 1.';

    const maxComp = maxComponentsInput.value ?? 0;
    if (maxComp < 1) return 'Max # components must be at least 1.';

    return null;
  }

  function refreshValidation(): void {
    const err = validate();
    validationDiv.textContent = err ?? '';
    runBtn.disabled = err != null;
  }
  refreshTooltip();

  // ---- Side grids with explicit selection-driven subsetting ----
  // Each side grid renders the DataFrame currently held in the corresponding table input. The
  // user picks rows in the grid (Ctrl/Shift+click), then presses "Subset by selection": we clone
  // the current DataFrame with the selection mask, register the clone with the workspace (the
  // table input is a choice widget backed by the workspace tables list — values not registered
  // are silently rejected and the input snaps back), set it as the input value, and remount the
  // grid against the subset. Previous subsets we created are closed to keep the workspace tidy.
  let suppressTemplatesChange = false;
  let suppressBbsChange = false;
  let suppressReagentsChange = false;
  let prevTemplatesSubset: DG.DataFrame | null = null;
  let prevBbsSubset: DG.DataFrame | null = null;
  let prevReagentsSubset: DG.DataFrame | null = null;

  const templatesGridHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', overflow: 'hidden'}});
  const bbsGridHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', overflow: 'hidden'}});
  const reagentsGridHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', overflow: 'hidden'}});

  const templatesSubsetStatus = ui.divText('',
    {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
  const bbsSubsetStatus = ui.divText('',
    {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
  const reagentsSubsetStatus = ui.divText('',
    {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});

  function updateTemplatesStatus(): void {
    const cur = templatesInput.value;
    templatesSubsetStatus.textContent = cur ? `${cur.rowCount} rows` : '';
  }
  function updateBbsStatus(): void {
    const cur = bbsInput.value;
    bbsSubsetStatus.textContent = cur ? `${cur.rowCount} rows` : '';
  }
  function updateReagentsStatus(): void {
    const cur = reagentsInput.value;
    reagentsSubsetStatus.textContent = cur ? `${cur.rowCount} rows` : '';
  }

  function mountTemplatesGrid(): void {
    templatesGridHost.innerHTML = '';
    const df = templatesInput.value;
    if (!df) {updateTemplatesStatus(); return;}
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    templatesGridHost.appendChild(grid.root);
    updateTemplatesStatus();
  }

  function mountBbsGrid(): void {
    bbsGridHost.innerHTML = '';
    const df = bbsInput.value;
    if (!df) {updateBbsStatus(); return;}
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    bbsGridHost.appendChild(grid.root);
    updateBbsStatus();
  }

  function mountReagentsGrid(): void {
    reagentsGridHost.innerHTML = '';
    const df = reagentsInput.value;
    if (!df) {
      reagentsGridHost.appendChild(ui.divText('No reagents file selected. Pick one in the Data ' +
        'section to enable reagents-mode enumeration (every step uses exactly one BB / earlier ' +
        'product plus reagents in the remaining slots).',
      {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}}));
      updateReagentsStatus();
      return;
    }
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    reagentsGridHost.appendChild(grid.root);
    updateReagentsStatus();
  }

  // ui.input.table is a Dart-side ChoiceInput whose items list is the workspace tables. Setting
  // .value to a DataFrame that's not in that list is silently rejected (the input snaps back to
  // the previous value). We side-step this by registering the clone with the workspace via
  // grok.shell.addTable before assigning, and we briefly dip the value through null to force the
  // widget to drop any cached pointer-equality with the previous selection.
  function assignTableInput(input: DG.InputBase<DG.DataFrame | null>, df: DG.DataFrame,
    setSuppress: (v: boolean) => void): void {
    grok.shell.addTable(df);
    setSuppress(true);
    try {
      try {(input as any).value = null;} catch {/* nullable: false rejects */}
      input.value = df;
    } finally {setSuppress(false);}
  }

  function subsetTemplatesBySelection(): void {
    const df = templatesInput.value;
    if (!df) return;
    const sel = df.selection;
    if (sel.trueCount === 0) {
      grok.shell.info('Select rows in the reaction templates grid first (Ctrl/Shift+click).');
      return;
    }
    if (sel.trueCount === df.rowCount) {
      grok.shell.info('All rows are selected — nothing to subset.');
      return;
    }
    const subset = df.clone(sel);
    subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
    detectChemSemTypes(subset);

    const prev = prevTemplatesSubset;
    prevTemplatesSubset = subset;
    assignTableInput(templatesInput, subset, (v) => {suppressTemplatesChange = v;});
    mountTemplatesGrid();
    refreshValidation();
    // Close the previous subset only after the input has switched away from it.
    if (prev && prev !== subset && prev !== df)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  function subsetBbsBySelection(): void {
    const df = bbsInput.value;
    if (!df) return;
    const sel = df.selection;
    if (sel.trueCount === 0) {
      grok.shell.info('Select rows in the building blocks grid first (Ctrl/Shift+click).');
      return;
    }
    if (sel.trueCount === df.rowCount) {
      grok.shell.info('All rows are selected — nothing to subset.');
      return;
    }
    const subset = df.clone(sel);
    subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
    detectChemSemTypes(subset);

    const prev = prevBbsSubset;
    prevBbsSubset = subset;
    assignTableInput(bbsInput, subset, (v) => {suppressBbsChange = v;});
    mountBbsGrid();
    refreshValidation();
    if (prev && prev !== subset && prev !== df)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  function subsetReagentsBySelection(): void {
    const df = reagentsInput.value;
    if (!df) {
      grok.shell.info('No reagents file selected.');
      return;
    }
    const sel = df.selection;
    if (sel.trueCount === 0) {
      grok.shell.info('Select rows in the reagents grid first (Ctrl/Shift+click).');
      return;
    }
    if (sel.trueCount === df.rowCount) {
      grok.shell.info('All rows are selected — nothing to subset.');
      return;
    }
    const subset = df.clone(sel);
    subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
    detectChemSemTypes(subset);

    const prev = prevReagentsSubset;
    prevReagentsSubset = subset;
    assignTableInput(reagentsInput, subset, (v) => {suppressReagentsChange = v;});
    mountReagentsGrid();
    refreshValidation();
    if (prev && prev !== subset && prev !== df)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  // When the user picks a different table or uploads a new CSV through the table input control,
  // detect chem semtypes (if needed) and re-mount the side grid. Programmatic value sets (our own
  // subset path) are guarded by `suppressTemplatesChange` / `suppressBbsChange`.
  templatesInput.onChanged.subscribe(() => {
    if (suppressTemplatesChange) return;
    const df = templatesInput.value;
    if (df) detectChemSemTypes(df);
    mountTemplatesGrid();
    refreshValidation();
  });
  bbsInput.onChanged.subscribe(() => {
    if (suppressBbsChange) return;
    const df = bbsInput.value;
    if (df) detectChemSemTypes(df);
    mountBbsGrid();
    refreshValidation();
  });
  reagentsInput.onChanged.subscribe(() => {
    if (suppressReagentsChange) return;
    const df = reagentsInput.value;
    if (df) detectChemSemTypes(df);
    mountReagentsGrid();
    refreshValidation();
    updateModeSummary();
  });

  // Re-validate on every input change so the Run button stays accurate.
  const wireValidation = <T>(input: DG.InputBase<T>) => input.onChanged.subscribe(() => {
    refreshTooltip();
    refreshValidation();
  });
  wireValidation(smartsColInput); wireValidation(blockingColInput); wireValidation(rxnNameColInput);
  wireValidation(bbColInput); wireValidation(reagentsColInput);
  wireValidation(exclusionInput); wireValidation(exclusionColInput);
  wireValidation(numRoundsInput); wireValidation(depthFirstInput);
  wireValidation(minCInput); wireValidation(maxCInput);
  wireValidation(maxCombosInput); wireValidation(maxComponentsInput);
  depthFirstInput.onChanged.subscribe(() => updateModeSummary());

  // ---- Buttons ----
  const editConfigBtn = ui.button('Edit all parameters', async () => {
    syncQuickInputsToConfig();
    const updated = await openConfigDialog(config);
    if (updated) {config = updated; syncConfigToQuickInputs(); refreshTooltip(); refreshValidation();}
  });
  ui.tooltip.bind(editConfigBtn, 'Open the full config form (every YAML field).');

  const loadYamlBtn = ui.button('Load YAML…', async () => {
    const f = await pickFile('.yaml,.yml');
    if (!f) return;
    try {
      const text = await f.text();
      config = configFromYaml(text);
      syncConfigToQuickInputs();
      refreshTooltip();
      refreshValidation();
      grok.shell.info(`Loaded config from ${f.name}.`);
    } catch (e) {
      grok.shell.error(`Could not load YAML: ${e instanceof Error ? e.message : String(e)}`);
    }
  });
  ui.tooltip.bind(loadYamlBtn, 'Load a YAML config file from disk and apply it to the form.');

  const saveYamlBtn = ui.button('Save YAML', () => {
    syncQuickInputsToConfig();
    downloadText(configToYaml(config), 'enumerator-config.yaml', 'text/yaml');
  });
  ui.tooltip.bind(saveYamlBtn, 'Download the current config as a YAML file.');

  // ---- Run / Cancel ----
  const progressLabel = ui.divText('', {style: {fontSize: '12px', color: 'var(--grey-5)'}});
  let cancelled = false;
  const runBtn = ui.bigButton('Run enumeration', async () => {
    if (validate() != null) return;
    syncQuickInputsToConfig();
    cancelled = false;
    runBtn.disabled = true;
    cancelBtn.style.display = '';
    progressLabel.textContent = 'Initializing…';
    try {
      await runFullEnumeration();
    } catch (e) {
      grok.shell.error(`Enumeration failed: ${e instanceof Error ? e.message : String(e)}`);
      console.error(e);
    } finally {
      runBtn.disabled = false;
      cancelBtn.style.display = 'none';
      progressLabel.textContent = '';
    }
  });
  ui.tooltip.bind(runBtn, 'Run the full enumeration with the current config and add the result to the workspace.');

  const cancelBtn = ui.button('Cancel', () => {cancelled = true;});
  cancelBtn.style.display = 'none';

  async function runFullEnumeration(): Promise<void> {
    progressLabel.textContent = 'Loading RDKit…';
    const rdkit = await getRdKit();
    const tDf = templatesInput.value!;
    const bDf = bbsInput.value!;
    const xDf = exclusionInput.value;
    const rDf = reagentsInput.value;
    const inputs = buildInputs(config, tDf, bDf, xDf, rDf);

    const reagentsPart = inputs.reagents.length > 0 ? ` × ${inputs.reagents.length} reagents` : '';
    progressLabel.textContent =
      `Running: ${inputs.templates.length} templates × ${inputs.buildingBlocks.length} BBs${reagentsPart} × ${config.enumeration.num_rounds} rounds`;
    const onProgress = (p: EnumerationProgress) => {
      const combo = p.combosTotal && p.combosTotal > 0 ?
        `, combos ${p.combosDone}/${p.combosTotal}` : '';
      progressLabel.textContent =
        `Round ${p.round}/${p.numRounds}, template ${p.templateIndex + 1}/${p.numTemplates}${combo}, products: ${p.productsSoFar}`;
    };

    const start = performance.now();
    const {rows, warnings} = await enumerate({
      rdkit, config, ...inputs, onProgress, isCancelled: () => cancelled,
    });
    const elapsed = ((performance.now() - start) / 1000).toFixed(1);

    if (cancelled)
      grok.shell.warning(`Enumeration cancelled. Partial results: ${rows.length} rows.`);
    else
      grok.shell.info(`Enumeration done in ${elapsed}s — ${rows.length} rows.`);
    if (warnings.length > 0) {
      console.warn('Enumeration warnings:', warnings);
      grok.shell.warning(`${warnings.length} warning(s); see console for details.`);
    }
    if (rows.length > 0) grok.shell.addTableView(buildResultDataFrame(rows));
  }

  // ---- Layout ----
  // Top-level: header (auto), main content (fills), validation (auto), action bar (auto).
  // The main content is a horizontal split: inputs on the left, side grids on the right.
  const headerRow = ui.divH([
    ui.divText('Chemical library enumeration', {style: {fontSize: '16px', fontWeight: 'bold'}}),
    appInfoIcon,
  ], {style: {alignItems: 'center', flex: '0 0 auto', marginBottom: '8px'}});

  // Three labelled sub-sections inside the Data column, each grouping the inputs that drive
  // one logical part of the run: reactions, starting materials, and the optional extras.
  // An optional trailing element (typically an info icon) renders inline with the title.
  const sectionHeader = (text: string, trailing?: HTMLElement): HTMLElement => {
    const t = ui.divText(text, {style: {fontWeight: 'bold', fontSize: '13px',
      color: 'var(--text-color)'}});
    if (!trailing) {
      t.style.marginTop = '10px';
      t.style.marginBottom = '4px';
      t.style.flex = '0 0 auto';
      return t;
    }
    return ui.divH([t, trailing],
      {style: {alignItems: 'center', marginTop: '10px', marginBottom: '4px', flex: '0 0 auto'}});
  };

  const dataSection = ui.divV([
    sectionHeader('Select Reactions'),
    ui.form([templatesInput, smartsColInput, blockingColInput, rxnNameColInput]),
    sectionHeader('Select Starting Materials'),
    ui.form([bbsInput, bbColInput]),
    sectionHeader('Additional settings'),
    ui.form([exclusionInput, exclusionColInput, reagentsInput, reagentsColInput]),
  ], {style: {flex: '0 0 auto'}});

  const configButtonsRow = ui.divH([editConfigBtn, loadYamlBtn, saveYamlBtn],
    {style: {gap: '8px', alignItems: 'center', marginTop: '6px', flex: '0 0 auto', flexWrap: 'wrap'}});

  // Plain-English summary of what the enumerator will do given the current selections — refreshed
  // whenever the reagents file or the depth-first toggle changes.
  const modeSummary = ui.divText('', {style: {fontSize: '11px', color: 'var(--grey-6)',
    marginTop: '6px', padding: '6px 8px', background: 'var(--grey-1)', borderRadius: '4px',
    border: '1px solid var(--grey-2)', flex: '0 0 auto', lineHeight: '1.4'}});

  function updateModeSummary(): void {
    const hasReagents = reagentsInput.value != null;
    if (hasReagents) {
      modeSummary.textContent = 'Reagents mode: each step combines EXACTLY ONE building block (or ' +
        'a product from an earlier round) with reagents in every remaining slot. Builds derivatives ' +
        'of each BB across rounds. Depth-first / breadth-first toggle is ignored.';
    } else if (depthFirstInput.value) {
      modeSummary.textContent = 'Depth-first mode: round 1 combines original BBs; round R > 1 ' +
        'extends each round-(R-1) product with original BBs (one prev-round product per step — ' +
        'linear chains, no convergent merging).';
    } else {
      modeSummary.textContent = 'Breadth-first mode: each round may combine any products from ' +
        'earlier rounds with BBs (convergent routes possible — the search space grows quickly).';
    }
  }

  const configSection = ui.divV([
    sectionHeader('Enumeration options', configInfoIcon),
    ui.form([numRoundsInput, depthFirstInput, maxComponentsInput, minCInput, maxCInput, maxCombosInput]),
    configButtonsRow,
    modeSummary,
  ], {style: {flex: '0 0 auto', marginTop: '8px'}});

  // Left pane scrolls vertically if it overflows so the action bar stays visible.
  const leftPane = ui.divV([dataSection, configSection],
    {style: {minWidth: '320px', overflowY: 'auto', overflowX: 'hidden', paddingRight: '8px'}});

  const selectionHint = (text: string) => ui.divText(text,
    {style: {fontSize: '11px', color: 'var(--grey-5)', marginBottom: '4px', flex: '0 0 auto'}});

  const panelHeader = (title: string, subsetBtn: HTMLElement, status: HTMLElement) =>
    ui.divH([
      ui.divText(title, {style: {fontWeight: 'bold', fontSize: '13px'}}),
      ui.div([], {style: {flex: '1 1 auto'}}),
      status, subsetBtn,
    ], {style: {alignItems: 'center', gap: '8px', flex: '0 0 auto', marginBottom: '4px'}});

  const templatesSubsetBtn = ui.button('Subset by selection', () => subsetTemplatesBySelection());
  ui.tooltip.bind(templatesSubsetBtn, 'Replace the reaction templates with only the rows currently ' +
    'selected in the grid. To restore the full set, reload the table from the input on the left.');

  const bbsSubsetBtn = ui.button('Subset by selection', () => subsetBbsBySelection());
  ui.tooltip.bind(bbsSubsetBtn, 'Replace the building blocks with only the rows currently selected ' +
    'in the grid. To restore the full set, reload the table from the input on the left.');

  const reagentsSubsetBtn = ui.button('Subset by selection', () => subsetReagentsBySelection());
  ui.tooltip.bind(reagentsSubsetBtn, 'Replace the reagents with only the rows currently selected ' +
    'in the grid. To restore the full set, reload the table from the input on the left.');

  const tabPanel = (header: HTMLElement, hint: string, gridHost: HTMLElement): HTMLElement =>
    ui.divV([
      header,
      selectionHint(hint),
      gridHost,
    ], {style: {height: '100%', display: 'flex', flexDirection: 'column', padding: '8px',
      background: 'var(--white)', boxSizing: 'border-box'}});

  const templatesPanel = tabPanel(
    panelHeader('Reaction templates', templatesSubsetBtn, templatesSubsetStatus),
    'Select rows (Ctrl/Shift+click) and click "Subset by selection" to enumerate only the ' +
      'chosen templates.',
    templatesGridHost);

  const bbsPanel = tabPanel(
    panelHeader('Building blocks', bbsSubsetBtn, bbsSubsetStatus),
    'Select rows (Ctrl/Shift+click) and click "Subset by selection" to enumerate only the ' +
      'chosen building blocks.',
    bbsGridHost);

  const reagentsPanel = tabPanel(
    panelHeader('Reagents', reagentsSubsetBtn, reagentsSubsetStatus),
    'Select rows (Ctrl/Shift+click) and click "Subset by selection" to enumerate only the ' +
      'chosen reagents.',
    reagentsGridHost);

  // ---- Preview tab (lazy) ----
  // Same idea as the old persistent preview: build inputs, run a budgeted enumerate, show the
  // result in a grid. Re-renders whenever the user switches to this tab. Subsequent activations
  // bump previewRunId so any in-flight preview is short-circuited via isCancelled.
  const PREVIEW_TARGET_ROWS = 20;
  const PREVIEW_MAX_COMBOS_PER_TEMPLATE = 3;
  const PREVIEW_MAX_ROUNDS = 2;
  const previewHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', overflow: 'hidden'}});
  const previewStatus = ui.divText('',
    {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
  const previewHeader = ui.divH([
    ui.divText('Sample products', {style: {fontWeight: 'bold', fontSize: '13px'}}),
    ui.div([], {style: {flex: '1 1 auto'}}),
    previewStatus,
  ], {style: {alignItems: 'center', gap: '8px', flex: '0 0 auto', marginBottom: '4px'}});
  const previewPanel = tabPanel(
    previewHeader,
    'Quick preview — runs a small subset of the enumeration (≤ 2 rounds, ≤ 3 combos / template) ' +
      'to give a flavour of the products before kicking off the full run.',
    previewHost);

  let previewRunId = 0;
  function showInPreview(content: HTMLElement | null): void {
    previewHost.innerHTML = '';
    if (content) previewHost.appendChild(content);
  }
  function shuffleInPlace<T>(arr: T[]): T[] {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [arr[i], arr[j]] = [arr[j], arr[i]];
    }
    return arr;
  }
  function pickPreviewSamples(rows: OutputRow[], n: number): OutputRow[] {
    if (rows.length <= n) return shuffleInPlace(rows.slice());
    const stepCount = (r: OutputRow) => r.route ? Math.max(0, r.route.split('>>').length - 1) : 0;
    const multi = shuffleInPlace(rows.filter((r) => stepCount(r) > 1));
    const single = shuffleInPlace(rows.filter((r) => stepCount(r) <= 1));
    const targetMulti = Math.min(multi.length, Math.ceil(n * 0.7));
    const combined = [...multi.slice(0, targetMulti), ...single.slice(0, n - targetMulti)];
    if (combined.length < n) combined.push(...multi.slice(targetMulti), ...single.slice(n - targetMulti));
    return shuffleInPlace(combined.slice(0, n));
  }

  async function refreshPreview(): Promise<void> {
    const myRunId = ++previewRunId;
    const err = validate();
    if (err) {
      previewStatus.textContent = '';
      showInPreview(ui.divText(`Fix the validation error first: ${err}`,
        {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}}));
      return;
    }
    previewStatus.textContent = 'running preview…';
    showInPreview(ui.divText('Computing preview…',
      {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}}));

    const tDf = templatesInput.value!;
    const bDf = bbsInput.value!;
    const xDf = exclusionInput.value;
    const rDf = reagentsInput.value;

    let inputs: BuiltInputs;
    try {
      inputs = buildInputs(config, tDf, bDf, xDf, rDf);
    } catch (e) {
      previewStatus.textContent = '';
      showInPreview(ui.divText(e instanceof Error ? e.message : String(e),
        {style: {color: 'var(--red-3)', padding: '20px', textAlign: 'center'}}));
      return;
    }

    let rdkit: any;
    try {rdkit = await getRdKit();} catch (e) {
      previewStatus.textContent = '';
      showInPreview(ui.divText(`Could not load RDKit: ${e instanceof Error ? e.message : String(e)}`,
        {style: {color: 'var(--red-3)', padding: '20px', textAlign: 'center'}}));
      return;
    }
    if (myRunId !== previewRunId) return;

    const previewConfig = cloneConfig(config);
    previewConfig.enumeration.num_rounds = Math.min(previewConfig.enumeration.num_rounds, PREVIEW_MAX_ROUNDS);
    previewConfig.max_num_combinations_per_template =
      previewConfig.max_num_combinations_per_template < 0 ?
        PREVIEW_MAX_COMBOS_PER_TEMPLATE :
        Math.min(previewConfig.max_num_combinations_per_template, PREVIEW_MAX_COMBOS_PER_TEMPLATE);
    previewConfig.max_num_routes_per_compound = 1;
    previewConfig.keep_building_blocks_in_final_output = false;

    let rows: OutputRow[] = [];
    try {
      const result = await enumerate({
        rdkit, config: previewConfig, ...inputs,
        isCancelled: () => myRunId !== previewRunId,
      });
      rows = result.rows;
    } catch (e) {
      if (myRunId !== previewRunId) return;
      previewStatus.textContent = '';
      showInPreview(ui.divText(`Preview failed: ${e instanceof Error ? e.message : String(e)}`,
        {style: {color: 'var(--red-3)', padding: '20px', textAlign: 'center'}}));
      return;
    }
    if (myRunId !== previewRunId) return;

    if (rows.length === 0) {
      previewStatus.textContent = '';
      showInPreview(ui.divText('No products produced within the preview budget. Try relaxing ' +
        'filters or verifying that templates and building blocks are compatible.',
      {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}}));
      return;
    }

    const samples = pickPreviewSamples(rows, PREVIEW_TARGET_ROWS);
    const df = buildResultDataFrame(samples, 'Preview');
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    try {
      grid.props.rowHeight = 110;
      const setW = (n: string, w: number) => {const c = grid.col(n); if (c) c.width = w;};
      setW('product', 180); setW('route', 460); setW('reaction_name', 140);
      setW('round', 56); setW('n_routes', 70); setW('template', 280);
    } catch (e) {
      console.warn('Preview grid styling failed:', e);
    }
    showInPreview(grid.root);
    previewStatus.textContent =
      `${samples.length} samples of ${rows.length} preview rows (≤ ${previewConfig.enumeration.num_rounds} rounds, ≤ ${PREVIEW_MAX_COMBOS_PER_TEMPLATE} combos / template)`;
  }

  // ---- Right pane: TabControl with lazy panes ----
  // addPane's factory is called once when the tab is first activated; we then keep the panel DOM
  // alive across switches. The preview tab additionally listens on onTabChanged to re-run whenever
  // the user comes back to it (so it picks up edits made on the left while another tab was shown).
  const tabs = ui.tabControl(null, false);
  tabs.addPane('Reaction templates', () => templatesPanel);
  tabs.addPane('Building blocks', () => bbsPanel);
  tabs.addPane('Reagents', () => reagentsPanel);
  tabs.addPane('Preview', () => previewPanel);

  tabs.onTabChanged.subscribe(() => {
    if (tabs.currentPane?.name === 'Preview') refreshPreview();
  });

  const rightPane = tabs.root;
  // rightPane.style.flex = '1 1 0';
  // rightPane.style.minWidth = '0';

  // Resizable horizontal split — drag the divider to rebalance inputs vs side grids.
  const mainRow = ui.splitH([leftPane, rightPane],
    {style: {flex: '1 1 0', minHeight: '0', width: '100%'}}, true);

  // Bottom action bar — Run/Cancel pinned at the bottom; config buttons live with the config form.
  const actionBar = ui.divH([runBtn, cancelBtn, progressLabel],
    {style: {gap: '8px', alignItems: 'center', marginTop: '10px', flex: '0 0 auto'}});

  const root = ui.divV([
    headerRow,
    mainRow,
    validationDiv,
    actionBar,
  ], {style: {padding: '16px', height: '100%', boxSizing: 'border-box', overflow: 'hidden'}});

  view.append(root);

  // Mount the initial grids, populate the mode summary, and run validation once everything is
  // wired up.
  mountTemplatesGrid();
  mountBbsGrid();
  mountReagentsGrid();
  updateModeSummary();
  refreshValidation();
  return view;
}
