/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {cloneConfig, configFromYaml, configToYaml, DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {openConfigDialog} from './config-form';
import {getRdKitModule} from '../chem-common-rdkit';
import {enumerate, EnumerationProgress, OutputRow, PerRoundOverride, TemplateInput} from './enumerate';

const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
const BUNDLED_BBS = 'enumerations/bb.csv';
const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

// Sniff string columns and set semType so the grid renders reactions and molecules: presence of
// `>>` in sampled values wins as ChemicalReaction, else auto-detection handles Molecule etc.
function detectChemSemTypes(df: DG.DataFrame): void {
  // detectSemanticTypes() scans the WHOLE dataframe; calling it per-column made this O(columns²)
  // and it ran on every step-clone. Tag ChemicalReaction columns first, then auto-detect once.
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
  }
  df.meta.detectSemanticTypes();
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

// Swap the column input's table in-place on parent-table change (preserves DOM identity so
// form layout stays valid), re-selecting the preferred column name if it exists in the new table.
function bindColToTable(
  colInput: DG.InputBase<DG.Column | null>, tableInput: DG.InputBase<DG.DataFrame | null>,
  getPreferredName: () => string, filter: (c: DG.Column) => boolean,
) {
  return tableInput.onChanged.subscribe(() => {
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

// Factored out of buildInputs so the same extraction also runs on a per-step subset table.
function extractTemplates(config: EnumeratorConfig, tDf: DG.DataFrame): TemplateInput[] {
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
  return templates;
}

function extractBuildingBlocks(config: EnumeratorConfig, bDf: DG.DataFrame): string[] {
  return getStringColumn(bDf, config.enumeration.bb_smiles_column).filter((s) => s.trim().length > 0);
}

function extractReagents(config: EnumeratorConfig, rDf: DG.DataFrame): string[] {
  const rCol = config.enumeration.reagent_smiles_column;
  return rDf.col(rCol) ? getStringColumn(rDf, rCol).filter((s) => s.trim().length > 0) : [];
}

function buildInputs(
  config: EnumeratorConfig, tDf: DG.DataFrame, bDf: DG.DataFrame,
  xDf: DG.DataFrame | null, rDf: DG.DataFrame | null,
): BuiltInputs {
  const templates = extractTemplates(config, tDf);
  const buildingBlocks = extractBuildingBlocks(config, bDf);

  let exclusionSmarts: string[] = [];
  if (xDf) {
    const excCol = config.products_specs.exclusion_smarts_products_file_smarts_col;
    if (xDf.col(excCol))
      exclusionSmarts = getStringColumn(xDf, excCol).filter((s) => s.trim().length > 0);
  }

  const reagents = rDf ? extractReagents(config, rDf) : [];

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

  // Re-bind column inputs whenever the parent table changes.
  view.subs.push(
    bindColToTable(smartsColInput, templatesInput, () => config.enumeration.smarts_col, isStringCol),
    bindColToTable(blockingColInput, templatesInput,
      () => config.enumeration.reactant_blocking_groups_per_template_column, isStringCol),
    bindColToTable(rxnNameColInput, templatesInput, () => config.enumeration.reaction_name_col, isStringCol),
    bindColToTable(bbColInput, bbsInput, () => config.enumeration.bb_smiles_column, isStringCol),
    bindColToTable(reagentsColInput, reagentsInput,
      () => config.enumeration.reagent_smiles_column, isStringCol),
    bindColToTable(exclusionColInput, exclusionInput,
      () => config.products_specs.exclusion_smarts_products_file_smarts_col, isStringCol),
  );

  // ---- CONFIG inputs ----
  const MAX_ROUNDS = 10;
  const numRoundsInput = ui.input.int('Number of rounds', {value: config.enumeration.num_rounds, min: 1, max: MAX_ROUNDS});
  numRoundsInput.setTooltip(
    'Number of consecutive enumeration rounds. Round 1 reacts BBs only; round 2 takes round-1 ' +
    `products and (in depth-first mode) reacts each one with original BBs. Increase for deeper ` +
    `libraries (capped at ${MAX_ROUNDS} — a step tab is built for every round, and product counts ` +
    `grow combinatorially with each one).`);

  const depthFirstInput = ui.input.bool('Depth first', {value: config.enumeration.depth_first});
  depthFirstInput.setTooltip(
    'When checked, each round-r > 1 step must combine EXACTLY ONE round-(r-1) product with original ' +
    'BBs (linear chain extension, no merging two complex products). Off (breadth-first) allows any ' +
    'combination from rounds 0..r-1 — typically explodes the search space and produces convergent routes.');

  // Promoted out of "Advanced limits & product filters" — used often enough to live at top level.
  const maxComponentsInput = ui.input.int('Max # components', {value: config.max_num_components, min: 1});
  maxComponentsInput.setTooltip('Max number of reactant components a template may have.');
  const maxRoutesInput = ui.input.int('Max routes per compound', {value: config.max_num_routes_per_compound});
  maxRoutesInput.setTooltip('Cap on the number of routes saved per product. Default -1 disables the cap.');

  // True while pushing config → inputs (syncConfigToQuickInputs). Each setAndFire fires onChanged,
  // which triggers a read-back (syncQuickInputsToConfig) mid-loop — before all inputs are updated —
  // overwriting config with stale values. The flag short-circuits that read-back during the sync.
  let pushingConfigToInputs = false;

  // ---- Info icons ----
  // appInfoIcon: what this app is / how it works. configInfoIcon: full current config as a card.
  // Both bind to live factories so hovering always reflects current state.
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
    const card = ui.div([], {style: {fontSize: '12px', maxWidth: '520px', lineHeight: '1.5', padding: '4px 2px'}});
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

  function currentMode(): 'depth' | 'breadth' | 'reagents' {
    return reagentsInput.value != null ? 'reagents' : (depthFirstInput.value ? 'depth' : 'breadth');
  }
  const MODE_LABEL = {depth: 'Depth-first', breadth: 'Breadth-first', reagents: 'Reagents'} as const;

  function buildConfigCard(): HTMLElement {
    const en = config.enumeration;
    const ps = config.products_specs;
    const mode = currentMode();
    const MODE_DESC = {
      depth: '(linear chain extension)',
      breadth: '(convergent allowed)',
      reagents: '(with reagents in other slots)',
    } as const;
    const modeLabel = `${MODE_LABEL[mode]} ${MODE_DESC[mode]}`;
    const fmtNum = (n: number, hint = 'unlimited') => n < 0 ? hint : String(n);
    const yn = (b: boolean) => b ? 'Yes' : 'No';

    const card = ui.div([], {style: {fontSize: '12px', maxHeight: '500px', maxWidth: '440px', overflow: 'auto', padding: '4px 2px', lineHeight: '1.5'}});

    const sectionTitle = (text: string) =>
      ui.divText(text, {style: {fontWeight: 'bold', marginTop: '8px', marginBottom: '2px', paddingBottom: '2px', borderBottom: '1px solid var(--grey-3)'}});
    const row = (label: string, value: string) => ui.divH([
      ui.divText(label, {style: {color: 'var(--grey-6)'}}),
      ui.divText(value, {style: {color: 'var(--text-color)', textAlign: 'right'}}),
    ], {style: {justifyContent: 'space-between', padding: '1px 0', gap: '12px'}});

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

  const syncQuickInputsToConfig = () => {
    config.enumeration.num_rounds = numRoundsInput.value ?? config.enumeration.num_rounds;
    config.enumeration.depth_first = !!depthFirstInput.value;
    config.max_num_components = maxComponentsInput.value ?? config.max_num_components;
    config.max_num_routes_per_compound = maxRoutesInput.value ?? config.max_num_routes_per_compound;
    // Column inputs hold a Column object; persist its name. Keep the previous value if unselected.
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

  // `input.value = X` updates the model but the Dart widget doesn't always re-render the visible
  // <input> text when set via API rather than typing — also push the value into the DOM element.
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
      // Clamp on load too — a hand-edited/older YAML could carry num_rounds above the UI's max.
      if (config.enumeration.num_rounds > MAX_ROUNDS) config.enumeration.num_rounds = MAX_ROUNDS;
      setAndFire(numRoundsInput, config.enumeration.num_rounds);
      setAndFire(depthFirstInput, config.enumeration.depth_first);
      setAndFire(maxComponentsInput, config.max_num_components);
      setAndFire(maxRoutesInput, config.max_num_routes_per_compound);
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
  const validationDiv = ui.divText('', {style: {color: 'var(--red-3)', fontSize: '12px', flex: '0 0 auto'}});

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
    if (rounds > MAX_ROUNDS) return `Number of rounds must be at most ${MAX_ROUNDS}.`;

    if (config.max_num_components < 1) return 'Max # components must be at least 1.';

    return null;
  }

  function refreshValidation(): void {
    refreshCfgRibbon();
    refreshStrategyCards();
    const err = validate();
    validationDiv.textContent = err ?? '';
    runBtn.disabled = err != null;
  }
  if (!pushingConfigToInputs) syncQuickInputsToConfig();

  // ---- Side grids with explicit selection-driven subsetting ----
  // "Subset by selection" clones the DataFrame with the selection mask, registers the clone with
  // the workspace (table input is a choice widget backed by the workspace list — unregistered
  // values are silently rejected), and remounts the grid. `original` remembers the user-loaded
  // full table so "Use all" can restore it without a file reload; resets on file swap.
  type SubsetState = {prev: DG.DataFrame | null; original: DG.DataFrame | null; suppress: boolean};
  const templatesState: SubsetState = {prev: null, original: null, suppress: false};
  const bbsState: SubsetState = {prev: null, original: null, suppress: false};
  const reagentsState: SubsetState = {prev: null, original: null, suppress: false};

  // ---- Per-step (per-round) subsetting state ----
  // Each component (0 = reactions, 1 = BBs, 2 = reagents) can be narrowed per round. We keep a
  // display-only clone (NOT registered in the workspace) per (component, round); its row SELECTION
  // is the per-step subset — empty/full selection means "use the global table". null = not built yet.
  // `stepState` lives INSIDE makeDataPanel (one array per panel, not one shared array indexed by
  // component) — see DataPanel.applyOverrideForRound for how buildPerRoundOverrides reaches in.

  // Tab row-count badges. Created upfront (they're standalone DOM + closure, no pane dependency) so
  // makeDataPanel can take them directly rather than through a nullable thunk; `.el` gets attached to
  // its actual pane header later, once the panes exist.
  type TabBadge = {el: HTMLSpanElement; refresh: (n: number | null) => void};
  const makeTabBadge = (): TabBadge => {
    const el = document.createElement('span');
    el.className = 'chem-enum-tab-badge';
    return {el, refresh: (n: number | null) => {
      el.textContent = n != null ? String(n) : '';
      el.style.display = n != null ? '' : 'none';
    }};
  };
  const templatesBadge = makeTabBadge();
  const bbsBadge = makeTabBadge();
  const reagentsBadge = makeTabBadge();

  const GRID_ROW_HEIGHT = 75;

  function applyGridColumnSizing(grid: DG.Grid, extendLast = true): void {
    try {
      grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
      if (extendLast) grid.props.extendLastColumn = true;
    } catch { /* setColumnsWidthType not available on older Dart builds */ }
  }

  // A single RAF after appending a resizable ui.splitH isn't enough for children to report a real
  // clientWidth (view may not be attached yet). ui.onSizeChanged is ResizeObserver-backed, so it
  // fires once real layout completes (and again on later resizes) instead of polling — subscribe,
  // size on the first nonzero total, then unsubscribe (this only needs to run once, on initial
  // layout). Tracked in view.subs so a view closed before that first size ever arrives doesn't leak.
  function sizeSplitOnceLaidOut(a: HTMLElement, b: HTMLElement, computeAWidth: (total: number) => number): void {
    const sub = ui.onSizeChanged(a).subscribe(() => {
      const total = a.clientWidth + b.clientWidth;
      if (total === 0) return;
      const aWidth = computeAWidth(total);
      a.style.width = aWidth + 'px';
      a.style.flexGrow = String(aWidth / total);
      b.style.width = (total - aWidth) + 'px';
      b.style.flexGrow = String((total - aWidth) / total);
      sub.unsubscribe();
    });
    view.subs.push(sub);
  }

  // REVERTED (live-verified): DG.Viewer.filters(df) was tried here — a standalone Viewer needs no
  // TableView at all — but with no explicit `filters`/`columnNames` config it renders a genuinely
  // EMPTY panel (confirmed live: toggling "Show filters" showed no per-column filters whatsoever).
  // `createDefaultFilters` (the auto-detect-and-populate behavior we actually need, including
  // substructure search for molecule/reaction columns) is exclusively a `TableView.getFiltersGroup()`
  // option (view.ts:404, bound to `api.grok_TableView_GetFilters`) — a different Dart call than the
  // standalone `Viewer.filters()` binding (`api.grok_Viewer_Filters`), despite view.ts's deprecated
  // `filters()` method looking like the same call at the JS layer. A TableView is genuinely required
  // for auto-populated defaults; back to that.
  //
  // mountDf creates a detached DG.TableView per call when withFilters is on; `host.innerHTML = ''`
  // only drops the DOM node, so the TableView leaks unless explicitly closed. Track the live view
  // per host and close the previous one before mounting a replacement.
  const liveFilterViews = new Map<HTMLElement, DG.TableView>();
  function closeLiveFilterView(host: HTMLElement): void {
    const prevTv = liveFilterViews.get(host);
    if (prevTv) {
      liveFilterViews.delete(host);
      try {prevTv.close();} catch (e) {console.warn('Could not close previous filters view:', e);}
    }
  }
  // Neither of the above runs on VIEW close — only on remount/no-table. Close every still-live
  // filters view when the view itself tears down, or their Dart-side resources leak for the session.
  view.subs.push(new Subscription(() => {
    for (const tv of liveFilterViews.values()) {
      try {tv.close();} catch (e) {console.warn('Could not close filters view on teardown:', e);}
    }
    liveFilterViews.clear();
  }));

  // `withFilters` adds the platform's real auto-detected default filters panel. A detached TableView
  // (addToWorkspace=false) never gets the shell's "added to dock tree" call, so grid/filters silently
  // fail to init (verified live: blank grid AND an empty filter panel — 0 filters, not just missing
  // substructure search). Firing `_onAdded()` manually (same trick as `mpo-create-profile.ts`) fixes
  // the grid, but ONLY fixes the filters too if it's DEFERRED a tick, exactly like
  // `mpo-create-profile.ts`'s own `setTimeout(() => this.tableView._onAdded(), 0)` — calling it
  // synchronously in the same tick as `DG.TableView.create()` (what this code did before) leaves
  // `getFiltersGroup({createDefaultFilters: true})` producing a real FilterGroup object with zero
  // populated filters inside it (confirmed live via DOM inspection: innerHTML empty vs. populated
  // after a deferred call). We keep only the grid + default FilterGroup, not tv.root.
  function mountDf(host: HTMLElement, df: DG.DataFrame, withFilters: boolean): void {
    closeLiveFilterView(host);
    host.innerHTML = '';
    const grid = DG.Viewer.grid(df);
    grid.props.rowHeight = GRID_ROW_HEIGHT;
    grid.root.style.cssText += ';width:100%;height:100%';
    if (!withFilters) {
      host.appendChild(grid.root);
      applyGridColumnSizing(grid);
      return;
    }
    const tv = DG.TableView.create(df, false);
    liveFilterViews.set(host, tv);
    const gridBox = ui.div([grid.root], {style: {flex: '1 1 0', minWidth: '0', height: '100%', overflow: 'hidden'}});
    const filtersBox = ui.div([], {style: {flex: '0 0 auto', width: '220px', height: '100%',
      overflow: 'auto', borderRight: '1px solid var(--grey-2)'}});
    const split = ui.splitH([filtersBox, gridBox], {style: {width: '100%', height: '100%', minHeight: '0'}}, true);
    host.appendChild(split);
    applyGridColumnSizing(grid);
    setTimeout(() => {
      // Bail if this host has moved on to a different mount since this was scheduled (fast tab
      // switching, or the panel closed) — liveFilterViews no longer pointing at this exact tv means
      // `host` isn't showing it anymore.
      if (liveFilterViews.get(host) !== tv) return;
      tv._onAdded();
      const filterGroup = tv.getFiltersGroup({createDefaultFilters: true});
      filterGroup.root.style.cssText += ';width:100%;height:100%;overflow:auto';
      filtersBox.appendChild(filterGroup.root);
      // ui.splitH ignores child width/flex style on first layout; its resize handler reads flexGrow
      // off the wrapper boxes it creates (children[0]/[2], skipping the divider at [1]) — set size there.
      const filtersWrap = split.children[0] as HTMLElement;
      const gridWrap = split.children[2] as HTMLElement;
      if (gridWrap && filtersWrap)
        sizeSplitOnceLaidOut(filtersWrap, gridWrap, (total) => Math.min(260, Math.round(total * 0.25)));
    }, 0);
  }

  // ui.input.table is a ChoiceInput over the workspace tables list; setting .value to a DataFrame
  // not in that list is silently rejected. Register via grok.shell.addTable first, and dip through
  // null to drop any cached pointer-equality with the previous selection.
  function assignTableInput(input: DG.InputBase<DG.DataFrame | null>, df: DG.DataFrame,
    setSuppress: (v: boolean) => void): void {
    grok.shell.addTable(df);
    setSuppress(true);
    try {
      try {input.value = null;} catch {/* nullable: false rejects */}
      input.value = df;
    } finally {setSuppress(false);}
  }

  // Shared clone/rename/detect core for both the global ("All steps") and per-step "Subset by
  // selection" actions — returns null (after an info toast) when there's nothing to subset.
  function cloneSubsetByRows(df: DG.DataFrame, emptyMsg: string): DG.DataFrame | null {
    const sel = df.selection;
    if (sel.trueCount === 0) {
      grok.shell.info(emptyMsg);
      return null;
    }
    if (sel.trueCount === df.rowCount) {
      grok.shell.info('All rows are selected — nothing to subset.');
      return null;
    }
    const subset = df.clone(sel);
    subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
    detectChemSemTypes(subset);
    return subset;
  }

  function subsetBySelection(
    input: DG.InputBase<DG.DataFrame | null>, state: SubsetState,
    mountFn: () => void, gridLabel: string, noTableMsg?: string,
  ): void {
    const df = input.value;
    if (!df) {
      if (noTableMsg) grok.shell.info(noTableMsg);
      return;
    }
    const subset = cloneSubsetByRows(df, `Select rows in the ${gridLabel} first (Ctrl/Shift+click).`);
    if (!subset) return;
    if (df !== state.prev) state.original = df; // remember the user's own table for "Use all"
    const prev = state.prev;
    state.prev = subset;
    assignTableInput(input, subset, (v) => { state.suppress = v; });
    mountFn();
    refreshValidation();
    // Close the previous subset only after the input has switched away from it.
    if (prev && prev !== subset && prev !== df)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  // Undo "Subset by selection": put the remembered full table back into the input and close the
  // subset clone. No-op (with a hint) when the full library is already in use.
  function restoreFullTable(
    input: DG.InputBase<DG.DataFrame | null>, state: SubsetState, mountFn: () => void, noun: string,
  ): void {
    const orig = state.original;
    if (!orig || input.value === orig) {
      grok.shell.info(`The full ${noun} library is already in use.`);
      return;
    }
    const prev = state.prev;
    state.prev = null;
    assignTableInput(input, orig, (v) => { state.suppress = v; });
    mountFn();
    refreshValidation();
    if (prev && prev !== orig)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  // Re-validate on every input change so the Run button stays accurate.
  const wireValidation = (input: DG.InputBase<unknown>): void => {
    view.subs.push(input.onChanged.subscribe(() => {
      refreshValidation();
    }));
  };
  [smartsColInput, blockingColInput, rxnNameColInput, bbColInput, reagentsColInput,
    exclusionInput, exclusionColInput, numRoundsInput, depthFirstInput,
    maxComponentsInput, maxRoutesInput,
  ].forEach((inp) => wireValidation(inp));

  // ---- Buttons ----
  const editConfigBtn = ui.button('Advanced limits & product filters', async () => {
    syncQuickInputsToConfig();
    const updated = await openConfigDialog(config);
    if (updated) {config = updated; syncConfigToQuickInputs(); refreshValidation();}
  });
  ui.tooltip.bind(editConfigBtn,
    'Caps (components, combinations, routes per compound) and product filters — every config field.');

  const loadYamlBtn = ui.button('Load YAML…', async () => {
    const f = await pickFile('.yaml,.yml');
    if (!f) return;
    try {
      const text = await f.text();
      config = configFromYaml(text);
      syncConfigToQuickInputs();
      refreshValidation();
      grok.shell.info(`Loaded config from ${f.name}.`);
    } catch (e) {
      grok.shell.error(`Could not load YAML: ${e instanceof Error ? e.message : String(e)}`);
    }
  });
  ui.tooltip.bind(loadYamlBtn, 'Load a YAML config file from disk and apply it to the form.');

  const saveYamlBtn = ui.button('Save YAML', () => {
    syncQuickInputsToConfig();
    DG.Utils.download('enumerator-config.yaml', configToYaml(config), 'text/yaml');
  });
  ui.tooltip.bind(saveYamlBtn, 'Download the current config as a YAML file.');

  // ---- Run / Cancel ----
  const progressLabel = ui.divText('', {style: {fontSize: '12px', color: 'var(--grey-5)'}});
  let cancelled = false;
  const runBtn = ui.bigButton('Enumerate', () => runWithUi(runEnumeration));
  ui.tooltip.bind(runBtn, 'Run enumeration with the current config and add the result to the workspace.');

  const cancelBtn = ui.button('Cancel', () => {cancelled = true;});
  cancelBtn.style.display = 'none';

  // Shared run chrome: validate, disable the Run button, show Cancel + progress, restore on finish.
  async function runWithUi(fn: () => Promise<void>): Promise<void> {
    if (validate() != null) return;
    syncQuickInputsToConfig();
    cancelled = false;
    runBtn.disabled = true;
    cancelBtn.style.display = '';
    progressLabel.textContent = 'Initializing…';
    try {
      await fn();
    } catch (e) {
      grok.shell.error(`Enumeration failed: ${e instanceof Error ? e.message : String(e)}`);
      console.error(e);
    } finally {
      runBtn.disabled = false;
      cancelBtn.style.display = 'none';
      progressLabel.textContent = '';
    }
  }

  async function runEnumeration(): Promise<void> {
    progressLabel.textContent = 'Loading RDKit…';
    const rdkit = await getRdKitModule();
    const tDf = templatesInput.value!;
    const bDf = bbsInput.value!;
    const xDf = exclusionInput.value;
    const rDf = reagentsInput.value;
    const inputs = buildInputs(config, tDf, bDf, xDf, rDf);

    const reagentsPart = inputs.reagents.length > 0 ? ` × ${inputs.reagents.length} reagents` : '';
    progressLabel.textContent =
      `Running: ${inputs.templates.length} templates × ${inputs.buildingBlocks.length} BBs${reagentsPart} × ${config.enumeration.num_rounds} round(s)`;
    const onProgress = (p: EnumerationProgress) => {
      const combo = p.combosTotal && p.combosTotal > 0 ?
        `, combos ${p.combosDone}/${p.combosTotal}` : '';
      progressLabel.textContent =
        `Round ${p.round}/${p.numRounds}, template ${p.templateIndex + 1}/${p.numTemplates}${combo}, products: ${p.productsSoFar}`;
    };

    const perRoundOverrides = buildPerRoundOverrides(config);
    if (perRoundOverrides)
      progressLabel.textContent += ' · per-step subsets active';

    const start = performance.now();
    const {rows, warnings} = await enumerate({
      rdkit, config, ...inputs, perRoundOverrides, onProgress, isCancelled: () => cancelled,
    });
    const elapsed = ((performance.now() - start) / 1000).toFixed(1);

    if (cancelled) {
      grok.shell.warning(`Enumeration cancelled. Partial results: ${rows.length} rows.`);
    } else {
      grok.shell.info(`Enumeration done in ${elapsed}s — ${rows.length} rows.`);
    }
    if (warnings.length > 0) {
      console.warn('Enumeration warnings:', warnings);
      // Surface the actual warning TEXT, not just a count — e.g. a per-step override silently not
      // applying is only visible this way, not as a bare count.
      const preview = warnings.slice(0, 3).join(' | ');
      const more = warnings.length > 3 ? ` (+${warnings.length - 3} more; see console)` : '';
      grok.shell.warning(`${preview}${more}`);
    }
    if (rows.length > 0) grok.shell.addTableView(buildResultDataFrame(rows));
  }

  // ---- Layout ----
  // Top-level: cfgRibbon (chips + run, auto), main content (fills), validation (auto).
  // The main content is a horizontal split: inputs on the left, side grids on the right.
  // The view title + info icon live in the view ribbon (setRibbonPanels).

  type StratCard = {root: HTMLElement; icon: HTMLElement};

  const applyStratCardStyle = (card: StratCard, mode: string, cur: string, enabled: boolean): void => {
    const sel = cur === mode;
    card.root.style.border = sel ? '2px solid var(--blue-2)' : '1px solid var(--grey-3)';
    card.root.style.opacity = enabled ? '1' : '0.5';
    card.root.style.cursor = enabled ? 'pointer' : 'not-allowed';
    card.icon.style.color = sel ? 'var(--blue-2)' : 'var(--grey-5)';
  };

  // Strategy cards replace the depth-first checkbox: depth/breadth drive the hidden `depthFirstInput`
  // (existing sync/validation keeps working); reagents card is active only when a reagents file is set.
  const buildStratCard = (icon: string, title: string, desc: string): StratCard => {
    const ic = ui.iconFA(icon);
    ic.style.marginTop = '2px';
    const root = ui.divH([ic, ui.divV([
      ui.divText(title, {style: {fontWeight: 'bold', fontSize: '12px'}}),
      ui.divText(desc, {style: {fontSize: '11px', color: 'var(--grey-6)', lineHeight: '1.3'}}),
    ], {style: {gap: '1px'}})], {style: {gap: '8px', alignItems: 'flex-start', padding: '8px 10px',
      border: '1px solid var(--grey-3)', borderRadius: '4px', cursor: 'pointer'}});
    return {root, icon: ic};
  };
  const stratDepthCard = buildStratCard('arrow-right', 'Depth-first',
    'Grow each product with original blocks — linear chains.');
  const stratBreadthCard = buildStratCard('sitemap', 'Breadth-first',
    'Combine any earlier products — convergent routes.');
  const stratReagentsCard = buildStratCard('flask', 'Reagents mode',
    'One block per step + reagents in remaining slots. Needs a reagents file.');

  function refreshStrategyCards(): void {
    const cur = currentMode();
    const hasReagents = cur === 'reagents';
    applyStratCardStyle(stratDepthCard, 'depth', cur, !hasReagents);
    applyStratCardStyle(stratBreadthCard, 'breadth', cur, !hasReagents);
    applyStratCardStyle(stratReagentsCard, 'reagents', cur, hasReagents);
  }
  stratDepthCard.root.onclick = (): void => {
    if (reagentsInput.value != null) return;
    setAndFire(depthFirstInput, true);
  };
  stratBreadthCard.root.onclick = (): void => {
    if (reagentsInput.value != null) return;
    setAndFire(depthFirstInput, false);
  };
  stratReagentsCard.root.onclick = (): void => {
    if (reagentsInput.value == null) {
      grok.shell.warning('Add a reagents file in “Extras” to switch to reagents mode.');
    }
  };

  // "Map columns" hides the optional template columns (blocking + reaction name) until needed.
  const mapColsBody = ui.form([blockingColInput, rxnNameColInput]);
  mapColsBody.style.display = 'none';
  const mapColsChevron = ui.iconFA('chevron-right');
  mapColsChevron.style.transition = 'transform 0.15s';
  const mapColsLink = ui.divH([mapColsChevron, ui.span([' Map columns (optional)'])],
    {style: {fontSize: '12px', color: 'var(--blue-2)', cursor: 'pointer', marginTop: '2px',
      gap: '2px', alignItems: 'center'}});
  let mapColsOpen = false;
  mapColsLink.onclick = (): void => {
    mapColsOpen = !mapColsOpen;
    mapColsBody.style.display = mapColsOpen ? '' : 'none';
    mapColsChevron.style.transform = mapColsOpen ? 'rotate(90deg)' : '';
  };

  editConfigBtn.style.width = '100%';
  const yamlRow = ui.divH([loadYamlBtn, saveYamlBtn],
    {style: {gap: '8px', alignItems: 'center', marginTop: '2px', flexWrap: 'wrap'}});

  // Right-pane tab references — assigned when tabs are built; used by section-open handlers for
  // context-sensitive tab switching. Declared here so openAccPaneExclusive can close over them.
  let templatesPane: DG.TabPane | undefined;
  let bbsPane: DG.TabPane | undefined;
  let reagentsPane: DG.TabPane | undefined;

  // Target pane resolved lazily via a thunk: Reactions pane is added expanded, so its factory runs
  // synchronously inside addPane, before later panes exist — capturing directly would hit the TDZ.
  const mkNextBtn = (getTarget: () => DG.AccordionPane): HTMLElement => {
    const btn = ui.button('Next →', () => openAccPaneExclusive(getTarget()));
    btn.classList.add('chem-enum-next-btn');
    return btn;
  };

  const accordion = ui.accordion();
  accordion.root.classList.add('chem-enum-accordion');
  const accReactionsPane = accordion.addPane('Reactions', () =>
    ui.divV([ui.form([templatesInput, smartsColInput]), mapColsLink, mapColsBody, mkNextBtn(() => accBbsPane)]), true);
  const accBbsPane = accordion.addPane('Building blocks', () =>
    ui.divV([ui.form([bbsInput, bbColInput]), mkNextBtn(() => accCombinePane)]), false);
  const accCombinePane = accordion.addPane('How to combine', () => ui.divV([
    ui.divH([
      ui.divText('Strategy', {style: {fontSize: '11px', color: 'var(--grey-6)', marginBottom: '2px'}}),
      configInfoIcon,
    ], {style: {alignItems: 'center', gap: '4px'}}),
    ui.divV([stratDepthCard.root, stratBreadthCard.root, stratReagentsCard.root], {style: {gap: '6px'}}),
    ui.form([numRoundsInput, maxComponentsInput, maxRoutesInput]),
    editConfigBtn,
    yamlRow,
    mkNextBtn(() => accExtrasPane),
  ], {style: {gap: '8px'}}), false);
  const accExtrasPane = accordion.addPane('Extras (optional)',
    () => ui.form([reagentsInput, reagentsColInput, exclusionInput, exclusionColInput]), false);
  const accPanes = [accReactionsPane, accBbsPane, accCombinePane, accExtrasPane];

  // Subtitle spans injected into each pane header — updated by refreshCfgRibbon().
  const injectPaneSub = (pane: DG.AccordionPane): HTMLElement => {
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    const sub = document.createElement('span');
    sub.className = 'chem-enum-pane-subtitle';
    header?.appendChild(sub);
    return sub;
  };
  const subReactions = injectPaneSub(accReactionsPane);
  const subBbs = injectPaneSub(accBbsPane);
  const subCombine = injectPaneSub(accCombinePane);
  const subExtras = injectPaneSub(accExtrasPane);

  // Left pane scrolls vertically if it overflows so the action bar stays visible.
  const leftPane = ui.divV([accordion.root],
    {style: {minWidth: '320px', overflowY: 'auto', overflowX: 'hidden', paddingRight: '8px'}});

  // === Config-summary ribbon (shown above the right-pane tabs) ===
  const cfgChipEl = (text: string): HTMLElement => {
    const chip = ui.divText(text);
    chip.className = 'chem-enum-chip';
    return chip;
  };
  const chipReactions = cfgChipEl('');
  const chipBbs = cfgChipEl('');
  const chipCombine = cfgChipEl('');
  const cfgEstEl = ui.divText('');
  cfgEstEl.className = 'chem-enum-cfg-est';
  const ribbonArrow = ui.iconFA('arrow-right');
  ribbonArrow.classList.add('chem-enum-ribbon-arrow');
  const runGroup = ui.divH([runBtn, cancelBtn, progressLabel],
    {style: {alignItems: 'center', gap: '6px'}});

  function refreshCfgRibbon(): void {
    const tDf = templatesInput.value; const bDf = bbsInput.value;
    const combineText = `${MODE_LABEL[currentMode()]} · ${numRoundsInput.value ?? 0} rounds`;
    const setChip = (chip: HTMLElement, text: string, err: boolean): void => {
      chip.textContent = text;
      chip.classList.toggle('chem-enum-chip--err', err);
    };
    setChip(chipReactions, tDf ? `${tDf.rowCount} reactions` : 'No reaction table', !tDf || !smartsColInput.value);
    setChip(chipBbs, bDf ? `${bDf.rowCount} Building Blocks` : 'No Building Blocks table', !bDf || !bbColInput.value);
    chipCombine.textContent = combineText;
    const n = (tDf && bDf) ? tDf.rowCount * bDf.rowCount : 0;
    cfgEstEl.textContent = n > 0 ? `≈ ${n.toLocaleString()} products` : '';
    subReactions.textContent = tDf ? `${tDf.rowCount} reactions` : 'No table selected';
    subBbs.textContent = bDf ? `${bDf.rowCount} building blocks` : 'No table selected';
    subCombine.textContent = combineText;
    subExtras.textContent = reagentsInput.value != null ? 'Reagents added' : '';
  }
  // === Exclusive accordion ===

  function switchTabForAccPane(pane: DG.AccordionPane): void {
    if (pane === accReactionsPane && templatesPane) {
      tabs.currentPane = templatesPane;
    } else if (pane === accBbsPane && bbsPane) {
      tabs.currentPane = bbsPane;
    } else if (pane === accExtrasPane && reagentsPane && reagentsInput.value != null) {
      tabs.currentPane = reagentsPane;
    }
  }

  function openAccPaneExclusive(pane: DG.AccordionPane): void {
    accPanes.forEach((p) => { p.expanded = p === pane; });
    switchTabForAccPane(pane);
  }

  // Wire native pane-header clicks for exclusive-open (Dart toggles pane.expanded first,
  // then our listener fires and collapses the others).
  accPanes.forEach((pane) => {
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    header?.addEventListener('click', () => { if (pane.expanded) openAccPaneExclusive(pane); });
  });

  chipReactions.onclick = () => openAccPaneExclusive(accReactionsPane);
  chipBbs.onclick = () => openAccPaneExclusive(accBbsPane);
  chipCombine.onclick = () => openAccPaneExclusive(accCombinePane);

  openAccPaneExclusive(accReactionsPane);

  const panelHeader = (hint: string, subsetBtn?: HTMLElement, status?: HTMLElement): HTMLElement => {
    const hintEl = ui.divText(hint, {style: {
      fontSize: '11px', color: 'var(--grey-5)', flex: '1 1 auto', marginRight: '4px',
    }});
    const children: HTMLElement[] = [hintEl];
    if (status) children.push(status);
    if (subsetBtn) children.push(subsetBtn);
    return ui.div(children, {style: {
      display: 'flex', alignItems: 'center', gap: '8px', flex: '0 0 auto',
      padding: '4px 8px 5px', borderBottom: '1px solid var(--grey-2)',
    }});
  };

  // Per-component "Subset by selection" now lives inside each tab's step bar (see makeDataPanel).

  const tabPanel = (header: HTMLElement, gridHost: HTMLElement): HTMLElement => {
    // display:flex on the host turns the grid's inline-flex outer display into a block-level
    // flex item, eliminating the 12px baseline-alignment gap that block+inline-flex produces.
    gridHost.style.cssText += ';position:relative;display:flex;flex-direction:column;flex:1 1 0;min-height:0;overflow:hidden';
    const fade = ui.div([], {style: {
      position: 'absolute', bottom: '0', left: '0', right: '0', height: '48px',
      background: 'linear-gradient(to bottom,transparent,var(--white))', pointerEvents: 'none', zIndex: '1',
    }});
    gridHost.appendChild(fade);
    return ui.div([header, gridHost], {style: {
      height: '100%', display: 'flex', flexDirection: 'column',
      background: 'var(--white)', boxSizing: 'border-box', overflow: 'hidden',
    }});
  };

  // ---- Single-grid per-component panel with a per-step strip ----
  // Each data tab shows ONE grid plus a horizontal step strip: "All steps" shows the full library
  // (with the global "Subset by selection"); "Step k" shows a display-only clone whose row selection
  // is that round's subset. Switching chips swaps what the single grid displays — no second grid.
  interface DataPanelOpts {
    idx: number; noun: string;
    input: DG.InputBase<DG.DataFrame | null>;
    state: SubsetState;
    badge: TabBadge;
    noTableMsg?: string;
    emptyMsg?: string;
    // How to fold this component's per-round work df into a PerRoundOverride fragment.
    apply: (o: PerRoundOverride, work: DG.DataFrame, cfg: EnumeratorConfig) => void;
  }
  interface DataPanel {
    panel: HTMLElement; render: () => void; onTableChanged: () => void; onRoundsChanged: () => void;
    // Applies this component's round-r override (if any) onto `out`; returns whether it did.
    // Lets buildPerRoundOverrides reach each panel's own stepState without a shared indexed array.
    applyOverrideForRound: (r: number, out: PerRoundOverride, cfg: EnumeratorConfig) => boolean;
  }
  function makeDataPanel(o: DataPanelOpts): DataPanel {
    let selStep = 0; // 0 = All steps (full library); 1..rounds = that round's subset
    let filtersOn = false; // funnel toggle: show the standard Datagrok filters panel next to the grid
    let currentDf: DG.DataFrame | null = null;
    // Per-round state (index k-1 = round k), one record per round — a single array instead of
    // parallel ones, so the df/subscription/committed-flag for a given round can never drift out of
    // sync with each other (they're always read and written together).
    // `committed`: explicit "did the user click Subset by selection for this round" flag — NOT
    // inferred from row-count comparison. A step's clone is a full snapshot of the global table
    // taken at whatever moment the user first visits that tab; if the global table is edited in
    // place afterward (row add/delete on the same DataFrame object — no onChanged fires, since the
    // object reference is unchanged), the clone's row count silently drifts from the global's.
    // Inferring "override" from that drift is wrong in both directions (a merely-visited, never-
    // committed clone could look like an override, or a real committed override could look like a
    // stale non-override) and was confirmed to let deleted rows quietly resurrect as part of a
    // round's reactant pool. The flag is set only by subsetStepBySelection and cleared only by
    // useAllForStep/a table swap.
    interface StepState { df: DG.DataFrame | null; sub: Subscription | null; committed: boolean; }
    const stepState: StepState[] = [];
    // Step selector: a real ui.tabControl. Each pane (All steps + Step k) owns its OWN persistent
    // barHost/gridHost, built once via that pane's own addPane content factory — nothing is shared
    // or relocated between panes. The grid itself is still destroyed and rebuilt on every switch
    // (renderGrid), so a shared, relocated node bought no reuse, only the "relocating a live grid
    // corrupts it" hazard the old design had to work around. `paneHosts[selStep]` gives renderBar/
    // renderGrid the current pane's own hosts; rebuilt fresh whenever buildStepTabs rebuilds the tabs.
    const stepTabsHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', display: 'flex',
      flexDirection: 'column', overflow: 'hidden'}});
    let stepTabsSub: Subscription | null = null;
    const stepDots: (HTMLElement | null)[] = []; // index k = 1..rounds; index 0 ("All steps") unused
    let paneHosts: {barHost: HTMLElement; gridHost: HTMLElement}[] = [];

    // A step "has an override" only once stepState[k-1].committed is explicitly set — see the flag's
    // own comment above for why row-count inference was wrong.
    const hasOverride = (k: number): boolean => {
      if (!stepState[k - 1]?.committed) return false;
      // Breadth-first ignores per-round BB overrides entirely (enumerate() draws round-R's reactant
      // pool from the union of all prior products, see `eligibleSmiles`) — don't show the override
      // dot/status as if it did something. Reactions/reagents overrides are unaffected. NOTE:
      // `currentMode()` treats any selected reagents file as "reagents mode"; enumerate()'s actual
      // gate additionally requires at least one reagent to survive canonicalization. A file that's
      // selected but entirely invalid is the one edge case where this dot can read stale — enumerate()
      // itself still warns correctly when that happens, so it's a display-lag risk, not silent data loss.
      if (o.idx === 1 && currentMode() === 'breadth') return false;
      return true;
    };
    // Capped defensively regardless of validation state — an invalid (too-high) input value still
    // blocks Run via `validate()`, but must not make buildStepTabs try to build hundreds of tabs.
    const roundCount = (): number =>
      Math.min(MAX_ROUNDS, Math.max(1, numRoundsInput.value ?? config.enumeration.num_rounds));
    const updateDots = (): void => {
      for (let k = 1; k <= stepDots.length - 1; k++) {
        const dot = stepDots[k];
        if (dot) dot.style.display = hasOverride(k) ? '' : 'none';
      }
    };
    // Backs the "per-step overrides about to be cleared" warning — only worth flagging when there's
    // actually something to lose.
    const hasAnyOverride = (): boolean => {
      for (let k = 1; k <= roundCount(); k++) if (hasOverride(k)) return true;
      return false;
    };

    // Shared bookkeeping behind every place a step's clone changes: unsubscribe the old selection
    // listener, subscribe the new one (if any), and replace the round's state record. Preserves the
    // existing `committed` value — callers that need to change it do so explicitly right after.
    const setStepWork = (k: number, work: DG.DataFrame | null): void => {
      stepState[k - 1]?.sub?.unsubscribe();
      const sub = work ? work.onSelectionChanged.subscribe(() => { updateDots(); renderBar(); }) : null;
      if (sub) view.subs.push(sub);
      stepState[k - 1] = {df: work, sub, committed: stepState[k - 1]?.committed ?? false};
    };

    const stepClone = (k: number): DG.DataFrame | null => {
      const existing = stepState[k - 1]?.df;
      if (existing) return existing;
      const global = o.input.value;
      if (!global) return null;
      // Display-only clone, never registered in the workspace; its selection carries the subset.
      const work = global.clone(null);
      work.name = `${global.name} · step ${k}`;
      detectChemSemTypes(work);
      setStepWork(k, work);
      return work;
    };

    // Per-step mirror of the global "Subset by selection": selecting rows in a step's grid is just
    // staging, the round only narrows once this commits it by swapping in a new clone.
    const subsetStepBySelection = (k: number): void => {
      const w = stepState[k - 1]?.df;
      if (!w) return; // grid must be mounted (via stepClone) before this button is reachable
      const subset = cloneSubsetByRows(w,
        `Select rows (Ctrl/Shift+click) to use only those ${o.noun} in step ${k} first.`);
      if (!subset) return;
      setStepWork(k, subset);
      stepState[k - 1].committed = true;
      renderGrid(); renderBar(); updateDots();
    };

    // Undo: drop the clone entirely so the step falls back to (re-derives from) "All steps" lazily.
    const useAllForStep = (k: number): void => {
      setStepWork(k, null);
      stepState[k - 1].committed = false;
      renderGrid(); renderBar(); updateDots();
    };

    // (Re)build the step tab strip, landing on `initialStep` (clamped). TabControl has no
    // add/remove-pane API, only clear()+re-addPane(); `tc.panes` insertion order lines up with selStep.
    function buildStepTabs(initialStep = 0): void {
      stepTabsSub?.unsubscribe();
      stepTabsHost.innerHTML = '';
      stepDots.length = 0;
      paneHosts = [];
      const tc = ui.tabControl(null, false);
      // tc.root must fill available height, not size to its header-strip content — each pane's own
      // content div (built below) needs the real space to lay its grid out in.
      tc.root.style.cssText += ';width:100%;flex:1 1 0;min-height:0;overflow:hidden';
      // Builds one pane's persistent content: its own barHost + gridHost, recorded in paneHosts at
      // its OWN fixed index k (not push-order) so renderBar/renderGrid can find "the currently
      // selected pane's own hosts" regardless of whether TabControl calls addPane's factory eagerly
      // for every pane up front or lazily on first visit — either way, position k always lands at
      // paneHosts[k], never at "whatever order the factories happened to fire in."
      const makePaneContent = (k: number): () => HTMLElement => () => {
        const barHost = ui.div([], {style: {display: 'flex', alignItems: 'center', gap: '8px', flex: '0 0 auto',
          padding: '4px 8px 5px', borderBottom: '1px solid var(--grey-2)'}});
        const gridHost = ui.div([], {style: {display: 'flex', flexDirection: 'column', flex: '1 1 0',
          minHeight: '0', overflow: 'hidden'}});
        paneHosts[k] = {barHost, gridHost};
        return ui.div([barHost, gridHost], {style: {
          height: '100%', display: 'flex', flexDirection: 'column', overflow: 'hidden'}});
      };
      const allPane = tc.addPane('All steps', makePaneContent(0));
      ui.tooltip.bind(allPane.header, 'Narrow this component for one round only. "All steps" is the full ' +
        'library used by every round; pick a step to restrict just that round. Per-step building-block ' +
        'subsets apply in depth-first / reagents mode; in breadth-first mode a round draws from all earlier ' +
        'products, so a BB subset has no effect. Resets when you swap the input file (in-range overrides ' +
        'survive a round-count change).');
      stepDots.push(null);
      for (let k = 1; k <= roundCount(); k++) {
        const pane = tc.addPane(`Step ${k}`, makePaneContent(k));
        // Position dot absolutely (header stacks children vertically, so marginLeft lands below the
        // label). Positive left offset only — a negative one bleeds into the adjacent tab (tabs sit
        // flush, live-verified). Fixed top offset, not top:50% — the header's box shrinks ~7px when
        // selected (underline indicator), so a percentage-based center drifts on selection.
        const dot = ui.div([], {style: {position: 'absolute', left: '5px', top: '12px',
          width: '6px', height: '6px', borderRadius: '50%',
          background: 'var(--orange-2, #c98a1b)', display: 'none'}});
        pane.header.style.cssText += ';position:relative';
        pane.header.appendChild(dot);
        stepDots.push(dot);
      }
      const renderCurrent = (): void => { renderGrid(); renderBar(); updateDots(); };
      stepTabsSub = tc.onTabChanged.subscribe(() => {
        // Index into tc.panes IS selStep (0 = "All steps", k = "Step k", by insertion order above)
        // — no need to parse the pane's label, which would break if it were ever reworded.
        const idx = tc.currentPane ? tc.panes.indexOf(tc.currentPane) : -1;
        selStep = idx < 0 ? 0 : idx;
        renderCurrent();
      });
      // Each rebuild's OWN unsubscribe (above, at the top of buildStepTabs) handles every prior
      // instance — but the LAST one built before the view closes has no later rebuild to retire it.
      // Track it in view.subs too so it isn't the one subscription in this panel that view close
      // never reaches (already-unsubscribed entries left behind here are inert, same tradeoff as
      // the per-step selection subs in setStepWork).
      view.subs.push(stepTabsSub);
      stepTabsHost.appendChild(tc.root);
      // Select explicitly — onTabChanged may not fire if the target is already the control's default.
      const clamped = Math.min(Math.max(0, initialStep), roundCount());
      selStep = clamped;
      const target = tc.panes[clamped] ?? allPane;
      if (target !== tc.currentPane) tc.currentPane = target;
      renderCurrent();
    }

    function renderBar(): void {
      const barHost = paneHosts[selStep]?.barHost;
      if (!barHost) return; // pane not built yet (shouldn't happen once buildStepTabs has run once)
      barHost.innerHTML = '';
      const hintEl = (t: string): HTMLElement =>
        ui.divText(t, {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '1 1 auto', marginRight: '4px'}});
      // Funnel toggle — shows the standard Datagrok filters panel for the visible grid. Off by default.
      const filterIcon = ui.iconFA('filter',
        () => { filtersOn = !filtersOn; renderGrid(); renderBar(); },
        filtersOn ? 'Hide filters' : 'Show filters');
      filterIcon.style.cssText += `;cursor:pointer;padding:2px 5px;flex:0 0 auto;` +
        `color:${filtersOn ? 'var(--blue-2)' : 'var(--grey-5)'}`;
      if (selStep === 0) {
        // Warn only when the click actually swaps the table AND overrides existed to lose — a no-op
        // click already gets its own info toast from subsetBySelection/restoreFullTable.
        const doGlobalAction = (action: () => void, clearedSuffix: string): void => {
          const hadOverride = hasAnyOverride();
          const prevValue = o.input.value;
          action();
          if (hadOverride && o.input.value !== prevValue)
            grok.shell.info(`Per-step ${o.noun} overrides were cleared — every step now uses the ${clearedSuffix}.`);
        };
        const btn = ui.button('Subset by selection', () => doGlobalAction(
          () => subsetBySelection(o.input, o.state, renderGrid, `${o.noun} grid`, o.noTableMsg),
          `new ${o.noun} subset`));
        ui.tooltip.bind(btn, `Replace the ${o.noun} library with only the selected rows (applies to every ` +
          `step). Click "Use all" to restore the full set.`);
        const useAll = ui.button('Use all', () => doGlobalAction(
          () => restoreFullTable(o.input, o.state, renderGrid, o.noun),
          `full ${o.noun} library again`));
        ui.tooltip.bind(useAll, `Restore the full ${o.noun} library (undo "Subset by selection").`);
        barHost.append(hintEl(`Full ${o.noun} library — used by every step unless a step overrides it.`),
          filterIcon, btn, useAll);
      } else {
        const w = stepState[selStep - 1]?.df;
        const status = ui.divText(
          w ? (hasOverride(selStep) ? `using ${w.rowCount} / ${o.input.value?.rowCount ?? w.rowCount}` : `all ${w.rowCount}`) : '',
          {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
        const btn = ui.button('Subset by selection', () => subsetStepBySelection(selStep));
        ui.tooltip.bind(btn, `Narrow step ${selStep} to only the selected rows (select rows with ` +
          `Ctrl/Shift+click first). Click "Use all" to go back to the full ${o.noun} library.`);
        const useAll = ui.button('Use all', () => useAllForStep(selStep));
        ui.tooltip.bind(useAll, `Undo "Subset by selection" so step ${selStep} uses the full ${o.noun} library ` +
          `(same as "All steps").`);
        barHost.append(
          hintEl(`Select rows (Ctrl/Shift+click), then "Subset by selection" to use only those ${o.noun} in ` +
            `step ${selStep}.`),
          status, filterIcon, btn, useAll);
      }
    }

    function renderGrid(): void {
      const gridHost = paneHosts[selStep]?.gridHost;
      if (!gridHost) return; // pane not built yet (shouldn't happen once buildStepTabs has run once)
      gridHost.innerHTML = '';
      currentDf = selStep === 0 ? o.input.value : stepClone(selStep);
      if (!currentDf) {
        closeLiveFilterView(gridHost); // this branch bypasses mountDf, so it must close it directly
        gridHost.appendChild(ui.divText(o.emptyMsg ?? `No ${o.noun} table selected.`,
          {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}}));
        if (selStep === 0) o.badge.refresh(null);
        return;
      }
      mountDf(gridHost, currentDf, filtersOn);
      if (selStep === 0) o.badge.refresh(currentDf.rowCount);
    }

    function render(): void { renderGrid(); renderBar(); updateDots(); }

    const onTableChanged = (): void => {
      if (o.input.value) detectChemSemTypes(o.input.value);
      // Guarded: subsetBySelection/restoreFullTable reassign o.input.value under suppress=true so
      // THIS handler won't treat their own programmatic swap as a genuine user file change and erase
      // "Use all"'s remembered original.
      if (!o.state.suppress) o.state.original = null;
      for (const s of stepState) s?.sub?.unsubscribe();
      stepState.length = 0;
      buildStepTabs(0);
      refreshValidation();
    };
    // Deliberately does NOT drop stepState entries for rounds beyond the new count. Dart int inputs
    // fire onChanged per keystroke, so typing "10" over "5" transiently passes through "1" — eagerly
    // truncating there would irreversibly destroy committed overrides on steps 2-5 before the user
    // finishes typing. buildStepTabs already only builds tabs up to roundCount() (capped at
    // MAX_ROUNDS), so out-of-range entries are simply invisible/unused, not leaked — nothing ever
    // writes past index MAX_ROUNDS-1 in the first place.
    const onRoundsChanged = (): void => buildStepTabs(selStep);

    const applyOverrideForRound = (r: number, out: PerRoundOverride, cfg: EnumeratorConfig): boolean => {
      const entry = stepState[r - 1];
      if (!entry?.committed || !entry.df) return false; // !entry.df shouldn't happen if committed
      const work = entry.df;
      o.apply(out, work, cfg);
      return true;
    };

    const panel = ui.div([stepTabsHost], {style: {
      height: '100%', display: 'flex', flexDirection: 'column', background: 'var(--white)', overflow: 'hidden'}});
    buildStepTabs(0); // also mounts the grid once hosts are in the DOM
    return {panel, render, onTableChanged, onRoundsChanged, applyOverrideForRound};
  }

  const templatesCtl = makeDataPanel({idx: 0, noun: 'reaction templates',
    input: templatesInput, state: templatesState, badge: templatesBadge,
    apply: (o, work, cfg) => { o.templates = extractTemplates(cfg, work); }});
  const bbsCtl = makeDataPanel({idx: 1, noun: 'building blocks',
    input: bbsInput, state: bbsState, badge: bbsBadge,
    apply: (o, work, cfg) => { o.buildingBlocks = extractBuildingBlocks(cfg, work); }});
  const reagentsCtl = makeDataPanel({idx: 2, noun: 'reagents',
    input: reagentsInput, state: reagentsState, badge: reagentsBadge,
    apply: (o, work, cfg) => { o.reagents = extractReagents(cfg, work); },
    noTableMsg: 'No reagents file selected.', emptyMsg: 'No reagents file selected. Add one in the Extras ' +
      'section to subset reagents per step.'});
  const templatesPanel = templatesCtl.panel;
  const bbsPanel = bbsCtl.panel;
  const reagentsPanel = reagentsCtl.panel;
  const dataCtls = [templatesCtl, bbsCtl, reagentsCtl];

  // Rebuilds step strips on round-count change or when a component's "All steps" table changes (file
  // load, Subset by selection, Use all) — narrowing "All steps" intentionally discards every step's
  // committed override, since "All steps" is the fallback source of truth. Runs unconditionally, not
  // gated by `*State.suppress` (that only guards `state.original` inside onTableChanged).
  //
  // `max` on an int input only shows a tooltip, it doesn't clamp the value — an over-max entry is
  // instead caught by validate() (which disables Run), same as every other invalid-input case in
  // this file; roundCount() separately caps at MAX_ROUNDS so an invalid value can't blow up tab count.
  view.subs.push(numRoundsInput.onChanged.subscribe(() => {
    refreshValidation();
    dataCtls.forEach((c) => c.onRoundsChanged());
  }));
  view.subs.push(templatesInput.onChanged.subscribe(() => templatesCtl.onTableChanged()));
  view.subs.push(bbsInput.onChanged.subscribe(() => bbsCtl.onTableChanged()));
  // BB override dot/status is mode-aware (hasOverride hides it in breadth-first) — re-render the
  // OTHER panels on any mode switch so they don't show stale state. reagentsCtl itself already gets
  // a full rebuild from onTableChanged below (its own table changed) — re-rendering it a second time
  // here would be redundant work, so it's excluded from the mode-refresh pass.
  view.subs.push(depthFirstInput.onChanged.subscribe(() => dataCtls.forEach((c) => c.render())));
  view.subs.push(reagentsInput.onChanged.subscribe(() => {
    reagentsCtl.onTableChanged();
    templatesCtl.render();
    bbsCtl.render();
  }));

  // A step's clone IS the subset once committed via "Subset by selection" — no deriving from a live
  // .selection bitset at run time. A round with no narrowed component falls back to the global set;
  // the whole result is undefined when nothing is overridden.
  function buildPerRoundOverrides(cfg: EnumeratorConfig): PerRoundOverride[] | undefined {
    const overrides: PerRoundOverride[] = [];
    let any = false;
    for (let r = 0; r < cfg.enumeration.num_rounds; r++) {
      const o: PerRoundOverride = {};
      for (const ctl of dataCtls) {
        if (ctl.applyOverrideForRound(r + 1, o, cfg)) any = true;
      }
      overrides.push(o);
    }
    return any ? overrides : undefined;
  }

  // ---- Preview tab (lazy) ----
  // Same idea as the old persistent preview: build inputs, run a budgeted enumerate, show the
  // result in a grid. Re-renders whenever the user switches to this tab. Subsequent activations
  // bump previewRunId so any in-flight preview is short-circuited via isCancelled.
  const PREVIEW_TARGET_ROWS = 20;
  const PREVIEW_MAX_COMBOS_PER_TEMPLATE = 3;
  const PREVIEW_MAX_ROUNDS = 2;
  const previewHost = ui.div([]);
  const previewStatus = ui.divText('',
    {style: {fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto'}});
  const previewHeader = panelHeader(
    'Quick preview — runs a small subset (≤ 2 rounds, ≤ 3 combos / template) to give a flavour of ' +
    'products. Uses the global reactions / building blocks / reagents; per-step subsets are not applied here.',
    undefined,
    previewStatus);
  const previewPanel = tabPanel(
    previewHeader,
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

    let rdkit: ReturnType<typeof getRdKitModule>;
    try {rdkit = await getRdKitModule();} catch (e) {
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
    } catch (e) {
      console.warn('Preview grid styling failed:', e);
    }
    showInPreview(grid.root);
    applyGridColumnSizing(grid, false); // route is not the last column — skip extendLastColumn
    previewStatus.textContent =
      `${samples.length} samples of ${rows.length} preview rows (≤ ${previewConfig.enumeration.num_rounds} rounds, ≤ ${PREVIEW_MAX_COMBOS_PER_TEMPLATE} combos / template)`;
  }

  // ---- Right pane: TabControl with lazy panes ----
  // addPane's factory is called once when the tab is first activated; we then keep the panel DOM
  // alive across switches.
  const tabs = ui.tabControl(null, false);
  tabs.root.style.cssText += ';width:100%;flex:1 1 0;min-height:0;overflow:hidden';
  templatesPane = tabs.addPane('Reaction templates', () => templatesPanel);
  bbsPane = tabs.addPane('Building blocks', () => bbsPanel);
  reagentsPane = tabs.addPane('Reagents', () => reagentsPanel);
  tabs.addPane('Preview', () => previewPanel);
  view.subs.push(tabs.onTabChanged.subscribe(() => {
    if (tabs.currentPane?.name === 'Preview') refreshPreview();
    // Bump previewRunId on tab-away too, or an in-flight preview keeps running unattended
    // (isCancelled only checked previewRunId bumps from starting a NEW preview before).
    else previewRunId++;
  }));
  // Row-count badges on each data tab header — updated whenever the underlying grid changes.
  templatesPane.header.appendChild(templatesBadge.el);
  bbsPane.header.appendChild(bbsBadge.el);
  reagentsPane.header.appendChild(reagentsBadge.el);

  // Right pane: a single component-tab control (each tab has its own step strip + one grid).
  const rightPane = ui.divV([tabs.root], {style: {height: '100%', overflow: 'hidden'}});

  // Resizable horizontal split — drag the divider to rebalance inputs vs side grids.
  const mainRow = ui.splitH([leftPane, rightPane],
    {style: {flex: '1 1 0', minHeight: '0', width: '100%'}}, true);
  // Initial split ~38/62: right pane gets more space for the grids.
  // Must target the wrapper boxes (mainRow.children[0/2]) — direct flex children of the split
  // container — because spliterResize() reads their flexGrow during drag to compute sumFlexGrow.
  const splitLeft = mainRow.children[0] as HTMLElement;
  const splitRight = mainRow.children[2] as HTMLElement;
  if (splitLeft && splitRight)
    sizeSplitOnceLaidOut(splitLeft, splitRight, (total) => Math.round(total * 0.38));

  const root = ui.divV([
    mainRow,
    validationDiv,
  ], {style: {padding: '0 0 0 16px', height: '100%', boxSizing: 'border-box', overflow: 'hidden'},
    classes: 'chem-enumerator'});

  view.setRibbonPanels([[appInfoIcon, runGroup, chipReactions, ribbonArrow, chipBbs, chipCombine, cfgEstEl]]);
  view.append(root);

  // Render each component panel (step strip + "All steps" grid) and run validation once wired up.
  templatesCtl.render();
  bbsCtl.render();
  reagentsCtl.render();
  refreshValidation();
  return view;
}
