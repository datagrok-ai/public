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
  const numRoundsInput = ui.input.int('Number of rounds', {value: config.enumeration.num_rounds, min: 1});
  numRoundsInput.setTooltip(
    'Number of consecutive enumeration rounds. Round 1 reacts BBs only; round 2 takes round-1 ' +
    'products and (in depth-first mode) reacts each one with original BBs. Increase for deeper libraries.');

  const depthFirstInput = ui.input.bool('Depth first', {value: config.enumeration.depth_first});
  depthFirstInput.setTooltip(
    'When checked, each round-r > 1 step must combine EXACTLY ONE round-(r-1) product with original ' +
    'BBs (linear chain extension, no merging two complex products). Off (breadth-first) allows any ' +
    'combination from rounds 0..r-1 — typically explodes the search space and produces convergent routes.');

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
    const mode = reagentsInput.value != null ? 'reagents' : (en.depth_first ? 'depth' : 'breadth');
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
  // Each side grid renders the DataFrame currently held in the corresponding table input. The
  // user picks rows in the grid (Ctrl/Shift+click), then presses "Subset by selection": we clone
  // the current DataFrame with the selection mask, register the clone with the workspace (the
  // table input is a choice widget backed by the workspace tables list — values not registered
  // are silently rejected and the input snaps back), set it as the input value, and remount the
  // grid against the subset. Previous subsets we created are closed to keep the workspace tidy.
  type SubsetState = {prev: DG.DataFrame | null; suppress: boolean};
  const templatesState: SubsetState = {prev: null, suppress: false};
  const bbsState: SubsetState = {prev: null, suppress: false};
  const reagentsState: SubsetState = {prev: null, suppress: false};

  const templatesGridHost = ui.div([]);
  const bbsGridHost = ui.div([]);
  const reagentsGridHost = ui.div([]);

  // Nullable refs to tab row-count badges — null until tabs are built; closures use optional chaining.
  type TabBadge = {el: HTMLSpanElement; refresh: (n: number | null) => void};
  let templatesBadge: TabBadge | null = null;
  let bbsBadge: TabBadge | null = null;
  let reagentsBadge: TabBadge | null = null;

  const GRID_ROW_HEIGHT = 75;
  const reagentsEmptyEl = ui.divText(
    'No reagents file selected. Pick one in the Data section to enable reagents-mode ' +
    'enumeration (every step uses exactly one BB / earlier product plus reagents in the ' +
    'remaining slots).',
    {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}},
  );

  function applyGridColumnSizing(grid: DG.Grid, extendLast = true): void {
    try {
      grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
      if (extendLast && grid.props.hasProperty('extendLastColumn'))
        (grid.props as any).extendLastColumn = true;
    } catch { /* setColumnsWidthType not available on older Dart builds */ }
  }

  function mountGrid(
    host: HTMLElement, input: DG.InputBase<DG.DataFrame | null>,
    badge: TabBadge | null, emptyEl?: HTMLElement,
  ): void {
    host.innerHTML = '';
    const df = input.value;
    if (!df) {
      if (emptyEl) host.appendChild(emptyEl);
      badge?.refresh(null);
      return;
    }
    const grid = DG.Viewer.grid(df);
    grid.props.rowHeight = GRID_ROW_HEIGHT;
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    host.appendChild(grid.root);
    applyGridColumnSizing(grid);
    badge?.refresh(df.rowCount);
  }

  const mountTemplates = (): void => mountGrid(templatesGridHost, templatesInput, templatesBadge);
  const mountBbs = (): void => mountGrid(bbsGridHost, bbsInput, bbsBadge);
  const mountReagents = (): void => mountGrid(reagentsGridHost, reagentsInput, reagentsBadge, reagentsEmptyEl);

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
      try {input.value = null;} catch {/* nullable: false rejects */}
      input.value = df;
    } finally {setSuppress(false);}
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
    const sel = df.selection;
    if (sel.trueCount === 0) {
      grok.shell.info(`Select rows in the ${gridLabel} first (Ctrl/Shift+click).`);
      return;
    }
    if (sel.trueCount === df.rowCount) {
      grok.shell.info('All rows are selected — nothing to subset.');
      return;
    }
    const subset = df.clone(sel);
    subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
    detectChemSemTypes(subset);
    const prev = state.prev;
    state.prev = subset;
    assignTableInput(input, subset, (v) => { state.suppress = v; });
    mountFn();
    refreshValidation();
    // Close the previous subset only after the input has switched away from it.
    if (prev && prev !== subset && prev !== df)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  // When the user picks a different table or uploads a new CSV through the table input control,
  // detect chem semtypes (if needed) and re-mount the side grid. Programmatic value sets (our own
  // subset path) are guarded by `templatesState.suppress` / `bbsState.suppress`.
  const wireTableInput = (
    input: DG.InputBase<DG.DataFrame | null>, state: SubsetState, mountFn: () => void,
  ): void => {
    view.subs.push(input.onChanged.subscribe(() => {
      if (state.suppress) return;
      const df = input.value;
      if (df) detectChemSemTypes(df);
      mountFn();
      refreshValidation();
    }));
  };
  wireTableInput(templatesInput, templatesState, mountTemplates);
  wireTableInput(bbsInput, bbsState, mountBbs);
  wireTableInput(reagentsInput, reagentsState, mountReagents);

  // Re-validate on every input change so the Run button stays accurate.
  const wireValidation = (input: DG.InputBase<unknown>): void => {
    view.subs.push(input.onChanged.subscribe(() => {
      refreshValidation();
    }));
  };
  [smartsColInput, blockingColInput, rxnNameColInput, bbColInput, reagentsColInput,
    exclusionInput, exclusionColInput, numRoundsInput, depthFirstInput,
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
  const runBtn = ui.bigButton('Enumerate', async () => {
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
  ui.tooltip.bind(runBtn, 'Run enumeration with the current config and add the result to the workspace.');

  const cancelBtn = ui.button('Cancel', () => {cancelled = true;});
  cancelBtn.style.display = 'none';

  async function runFullEnumeration(): Promise<void> {
    progressLabel.textContent = 'Loading RDKit…';
    const rdkit = await getRdKitModule();
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

    if (cancelled) {
      grok.shell.warning(`Enumeration cancelled. Partial results: ${rows.length} rows.`);
    } else {
      grok.shell.info(`Enumeration done in ${elapsed}s — ${rows.length} rows.`);
    }
    if (warnings.length > 0) {
      console.warn('Enumeration warnings:', warnings);
      grok.shell.warning(`${warnings.length} warning(s); see console for details.`);
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

  // Strategy cards replace the depth-first checkbox: three explicit choices. depth/breadth drive the
  // hidden `depthFirstInput` (so all existing sync/validation keeps working); the reagents card is
  // active only when a reagents file is selected in Extras (reagents-mode follows its presence).
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
    ui.form([numRoundsInput]),
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

  const templatesSubsetBtn = ui.button('Subset by selection', () =>
    subsetBySelection(templatesInput, templatesState, mountTemplates, 'reaction templates grid'));
  ui.tooltip.bind(templatesSubsetBtn, 'Replace the reaction templates with only the rows currently ' +
    'selected in the grid. To restore the full set, reload the table from the input on the left.');

  const bbsSubsetBtn = ui.button('Subset by selection', () =>
    subsetBySelection(bbsInput, bbsState, mountBbs, 'building blocks grid'));
  ui.tooltip.bind(bbsSubsetBtn, 'Replace the building blocks with only the rows currently selected ' +
    'in the grid. To restore the full set, reload the table from the input on the left.');

  const reagentsSubsetBtn = ui.button('Subset by selection', () =>
    subsetBySelection(reagentsInput, reagentsState, mountReagents, 'reagents grid', 'No reagents file selected.'));
  ui.tooltip.bind(reagentsSubsetBtn, 'Replace the reagents with only the rows currently selected ' +
    'in the grid. To restore the full set, reload the table from the input on the left.');

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

  const subsetHint = (what: string) =>
    `Select rows (Ctrl/Shift+click) and click "Subset by selection" to enumerate only chosen ${what}.`;

  const templatesPanel = tabPanel(panelHeader(subsetHint('templates'), templatesSubsetBtn), templatesGridHost);
  const bbsPanel = tabPanel(panelHeader(subsetHint('building blocks'), bbsSubsetBtn), bbsGridHost);
  const reagentsPanel = tabPanel(panelHeader(subsetHint('reagents'), reagentsSubsetBtn), reagentsGridHost);

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
    'Quick preview — runs a small subset (≤ 2 rounds, ≤ 3 combos / template) to give a flavour of products.',
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
  // alive across switches. The preview tab additionally listens on onTabChanged to re-run whenever
  // the user comes back to it (so it picks up edits made on the left while another tab was shown).
  const tabs = ui.tabControl(null, false);
  tabs.root.style.cssText += ';width:100%;flex:1 1 0;min-height:0;overflow:hidden';
  templatesPane = tabs.addPane('Reaction templates', () => templatesPanel);
  bbsPane = tabs.addPane('Building blocks', () => bbsPanel);
  reagentsPane = tabs.addPane('Reagents', () => reagentsPanel);
  tabs.addPane('Preview', () => previewPanel);
  // Row-count badges on each data tab header — updated whenever the underlying grid changes.
  const makeTabBadge = (): TabBadge => {
    const el = document.createElement('span');
    el.className = 'chem-enum-tab-badge';
    return {el, refresh: (n: number | null) => {
      el.textContent = n != null ? String(n) : '';
      el.style.display = n != null ? '' : 'none';
    }};
  };
  templatesBadge = makeTabBadge();
  bbsBadge = makeTabBadge();
  reagentsBadge = makeTabBadge();
  templatesPane.header.appendChild(templatesBadge.el);
  bbsPane.header.appendChild(bbsBadge.el);
  reagentsPane.header.appendChild(reagentsBadge.el);

  view.subs.push(tabs.onTabChanged.subscribe(() => {
    if (tabs.currentPane?.name === 'Preview') refreshPreview();
  }));

  const rightPane = ui.divV([tabs.root], {style: {height: '100%', overflow: 'hidden'}});

  // Resizable horizontal split — drag the divider to rebalance inputs vs side grids.
  const mainRow = ui.splitH([leftPane, rightPane],
    {style: {flex: '1 1 0', minHeight: '0', width: '100%'}}, true);
  // Initial split ~38/62: right pane gets more space for the grids.
  // Must target the wrapper boxes (mainRow.children[0/2]) — direct flex children of the split
  // container — because spliterResize() reads their flexGrow during drag to compute sumFlexGrow.
  requestAnimationFrame(() => {
    const splitLeft = mainRow.children[0] as HTMLElement;
    const splitRight = mainRow.children[2] as HTMLElement;
    if (splitLeft && splitRight) {
      const total = splitLeft.clientWidth + splitRight.clientWidth;
      if (total > 0) {
        splitLeft.style.width = Math.round(total * 0.38) + 'px';
        splitLeft.style.flexGrow = '0.38';
        splitRight.style.width = Math.round(total * 0.62) + 'px';
        splitRight.style.flexGrow = '0.62';
      }
    }
  });

  const root = ui.divV([
    mainRow,
    validationDiv,
  ], {style: {padding: '0 0 0 16px', height: '100%', boxSizing: 'border-box', overflow: 'hidden'},
    classes: 'chem-enumerator'});

  view.setRibbonPanels([[appInfoIcon, runGroup, chipReactions, ribbonArrow, chipBbs, chipCombine, cfgEstEl]]);
  view.append(root);

  // Mount the initial grids and run validation once everything is wired up.
  mountTemplates();
  mountBbs();
  mountReagents();
  refreshValidation();
  return view;
}
