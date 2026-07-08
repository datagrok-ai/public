/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import {Subscription} from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package, PackageFunctions} from '../../package';
import {cloneConfig, configFromYaml, configToYaml, DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {openConfigDialog} from './config-form';
import {getRdKitModule} from '../chem-common-rdkit';
import {enumerate, EnumerationProgress, OutputRow, PerRoundOverride, TemplateInput, tryGetRxn} from './enumerate';

const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
const BUNDLED_BBS = 'enumerations/bb.csv';
const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

// A freshly-(re)mounted Filters viewer can asynchronously reapply a stale per-column categorical
// selection over a DataFrame's .filter BitSet shortly after construction (a Datagrok platform
// behavior, not something this app triggers) — this is how long every reset-after-remount in this
// file waits before re-asserting the intended filter state, past that clobber window.
const FILTER_REMOUNT_SETTLE_MS = 200;

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

// A column being SELECTED doesn't mean it holds the right data — e.g. the user picks "Name" where
// "SMILES" was auto-detected. DG.Detector.sampleCategories is the same sampler Chem's own detectors
// use (detectSmiles/detectReactions in detectors.js) — reuse it instead of a bespoke row scanner.
// validate() reruns on every tracked input's onChanged, including ones unrelated to file/column
// choice (rounds, depth-first, ...); a WeakMap keyed on the Column object itself means an RDKit
// re-parse only happens the first time a given column is checked, not on every unrelated keystroke,
// and it never needs manual invalidation — picking a different column or file yields a new Column.
const columnContentValidCache = new WeakMap<DG.Column, boolean>();
function cachedSampleValid(col: DG.Column, isValid: (s: string) => boolean): boolean {
  let result = columnContentValidCache.get(col);
  if (result === undefined) {
    result = DG.Detector.sampleCategories(col, isValid, 1, 5, 0.8);
    columnContentValidCache.set(col, result);
  }
  return result;
}

function isValidReactionSmarts(s: string): boolean {
  const rxn = tryGetRxn(getRdKitModule(), s);
  const ok = !!rxn;
  rxn?.delete();
  return ok;
}

// validateMolecule returns '' (not null) for a valid molecule — only a non-empty string is an error.
function isValidSmiles(s: string): boolean {
  return !PackageFunctions.validateMolecule(s);
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
  // `min`/`max` above only affect the tooltip/spinner, not validation — add that separately.
  numRoundsInput.addValidator((v) => {
    const n = Number(v);
    if (!Number.isFinite(n) || n < 1) return 'Must be at least 1.';
    if (n > MAX_ROUNDS) return `Must be at most ${MAX_ROUNDS} — the step strip shows the first ${MAX_ROUNDS} rounds only.`;
    return null;
  });

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
    // A table with no rows (e.g. an unparseable/empty upload that still produced a valid, if
    // degenerate, DataFrame) is truthy and passes a plain null check — Enumerate stayed clickable
    // and silently did nothing. Row-count checks close that gap.
    const tDf = templatesInput.value;
    if (!tDf) return 'Select a reaction templates file.';
    if (tDf.rowCount === 0) return 'Reaction templates file has no rows.';
    if (!smartsColInput.value) return 'Select a reaction SMARTS column.';
    if (!cachedSampleValid(smartsColInput.value, isValidReactionSmarts))
      return 'Selected column does not contain valid reaction SMARTS — pick a different column.';

    const bDf = bbsInput.value;
    if (!bDf) return 'Select a building blocks file.';
    if (bDf.rowCount === 0) return 'Building blocks file has no rows.';
    if (!bbColInput.value)
      return 'Building blocks are missing, or the wrong column is selected.';
    if (!cachedSampleValid(bbColInput.value, isValidSmiles))
      return 'Selected column does not contain valid SMILES — pick a different column.';

    const rDf = reagentsInput.value;
    if (rDf && rDf.rowCount === 0) return 'Reagents file has no rows — clear it or pick a non-empty one.';
    if (rDf && !reagentsColInput.value)
      return 'Select a reagent SMILES column or clear the reagents file.';

    const xDf = exclusionInput.value;
    if (xDf && xDf.rowCount === 0)
      return 'Exclusion substructures file has no rows — clear it or pick a non-empty one.';
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
    bindRunTooltip(err);
  }
  if (!pushingConfigToInputs) syncQuickInputsToConfig();

  // ---- Side grids with explicit selection-driven subsetting ----
  // "Subset by selection" clones by the selection mask and registers the clone with the workspace,
  // since the table input's choice widget silently rejects unregistered values. `original` restores
  // the user's own table for "Use all".
  type SubsetState = {prev: DG.DataFrame | null; original: DG.DataFrame | null};
  const templatesState: SubsetState = {prev: null, original: null};
  const bbsState: SubsetState = {prev: null, original: null};
  const reagentsState: SubsetState = {prev: null, original: null};

  // ---- Per-step (per-round) subsetting state ----
  // Each component (reactions/BBs/reagents) can be narrowed per round via a display-only clone
  // whose row selection is that round's subset. `stepState` lives inside makeDataPanel (one array
  // per panel) — see DataPanel.applyOverrideForRound.

  // Tab row-count badge. Reactions/BBs already show their row count via the always-visible ribbon
  // chips (chipReactions/chipBbs) and the accordion pane subtitles — a tab badge there would just
  // repeat the same number a third time. Reagents has neither, so it keeps a badge.
  type TabBadge = {el: HTMLSpanElement; refresh: (n: number | null) => void};
  const makeTabBadge = (): TabBadge => {
    const el = document.createElement('span');
    el.className = 'chem-enum-tab-badge';
    return {el, refresh: (n: number | null) => {
      el.textContent = n != null ? String(n) : '';
      el.style.display = n != null ? '' : 'none';
    }};
  };
  const reagentsBadge = makeTabBadge();

  const GRID_ROW_HEIGHT = 75;

  function applyGridColumnSizing(grid: DG.Grid, extendLast = true): void {
    try {
      grid.setColumnsWidthType(DG.ColumnWidthType.Optimal);
      if (extendLast) grid.props.extendLastColumn = true;
    } catch { /* setColumnsWidthType not available on older Dart builds */ }
  }

  // A single RAF after appending a resizable ui.splitH isn't reliably enough for a real clientWidth
  // yet — ui.onSizeChanged (ResizeObserver-backed) fires once real layout completes; size once,
  // then unsubscribe. Tracked in view.subs in case the view closes before that first size arrives.
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

  // Plain viewers (grid, filters) leak (host.innerHTML='' only drops the DOM node) unless closed
  // explicitly — track the live ones per host and close them before mounting a replacement.
  const mountedViewers = new Map<HTMLElement, DG.Viewer[]>();
  function closeMountedViewers(host: HTMLElement): void {
    const prev = mountedViewers.get(host);
    if (!prev) return;
    mountedViewers.delete(host);
    for (const v of prev) {
      // A viewer mid-Dart-side-construction throws TypeError on close — expected on initial
      // load, not actionable.
      try {v.close();} catch (e) {
        if (!(e instanceof TypeError)) console.warn('Could not close previous viewer:', e);
      }
    }
  }
  // Neither of the above runs on VIEW close — only on remount/no-table. Close every still-mounted
  // viewer when the view itself tears down, or their Dart-side resources leak for the session.
  view.subs.push(new Subscription(() => {
    for (const host of mountedViewers.keys()) closeMountedViewers(host);
  }));

  // ChemicalReaction has no meaningful substructure-filter semantics for a whole reaction template,
  // so it's simply left out of the filters below (not a workaround for anything — DG.Viewer.filters
  // only builds filters for the columns it's explicitly given).
  function mountDf(host: HTMLElement, df: DG.DataFrame, withFilters: boolean): void {
    closeMountedViewers(host);
    host.innerHTML = '';
    const grid = DG.Viewer.grid(df);
    grid.props.rowHeight = GRID_ROW_HEIGHT;
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    if (!withFilters) {
      mountedViewers.set(host, [grid]);
      host.appendChild(grid.root);
      applyGridColumnSizing(grid);
      return;
    }
    const filterStates = df.columns.toList()
      // `~`-prefixed columns are internal/technical (e.g. `~SMILES.Pattern`, added by the
      // substructure filter's own fingerprint cache) — never surface a filter for those.
      .filter((col) => col.semType !== 'ChemicalReaction' && !col.name.startsWith('~'))
      .map((col) => ({
        // A numeric column with zero variance (every value equal — common right after subsetting to
        // the exact value it was filtered on) makes the histogram filter widget clobber the whole
        // dataframe's .filter to all-false, both on construction and again if the user later toggles
        // it histogram<->categorical (live-verified — unrelated to and not fixed by the subset-clone
        // filter reset above). A categorical filter has no degenerate-range code path, so route
        // zero-variance numeric columns there instead of ever constructing the buggy histogram widget.
        type: col.isNumerical ? (col.stats.min === col.stats.max ? DG.FILTER_TYPE.CATEGORICAL : DG.FILTER_TYPE.HISTOGRAM) :
          col.semType === 'Molecule' ? DG.FILTER_TYPE.SUBSTRUCTURE : DG.FILTER_TYPE.CATEGORICAL,
        column: col.name,
      }));
    const filtersViewer = DG.Viewer.filters(df, {filters: filterStates});
    mountedViewers.set(host, [grid, filtersViewer]);
    filtersViewer.root.style.width = '100%';
    filtersViewer.root.style.height = '100%';
    filtersViewer.root.style.overflow = 'auto';
    const gridBox = ui.div([grid.root], {style: {flex: '1 1 0', minWidth: '0', height: '100%', overflow: 'hidden'}});
    // width must track the wrapper sizeSplitOnceLaidOut resizes below (split.children[0]), not a fixed
    // px value — a fixed width here left a dead gap between this box and the divider whenever the
    // wrapper's computed width (total * 0.25, capped 260) exceeded it (live-verified).
    const filtersBox = ui.div([filtersViewer.root], {style: {flex: '0 0 auto', width: '100%', height: '100%',
      overflow: 'auto', borderRight: '1px solid var(--grey-2)'}});
    const split = ui.splitH([filtersBox, gridBox], {style: {width: '100%', height: '100%', minHeight: '0'}}, true);
    host.appendChild(split);
    applyGridColumnSizing(grid);
    // ui.splitH ignores child width/flex style on first layout; its resize handler reads flexGrow
    // off the wrapper boxes it creates (children[0]/[2], skipping the divider at [1]) — set size there.
    const filtersWrap = split.children[0] as HTMLElement;
    const gridWrap = split.children[2] as HTMLElement;
    if (gridWrap && filtersWrap)
      sizeSplitOnceLaidOut(filtersWrap, gridWrap, (total) => Math.min(260, Math.round(total * 0.25)));
  }

  // ui.input.table is a ChoiceInput over the workspace tables list; setting .value to a DataFrame
  // not in that list is silently rejected. Register via grok.shell.addTable first, and dip through
  // null to drop any cached pointer-equality with the previous selection.
  function assignTableInput(input: DG.InputBase<DG.DataFrame | null>, df: DG.DataFrame): void {
    grok.shell.addTable(df);
    try {input.value = null;} catch {/* nullable: false rejects */}
    input.value = df;
  }

  // Shared wording for the two identically-phrased "how to subset" prompts below (per-step's empty-
  // selection message and its hint text) — kept as one constant so the two stay in sync.
  const SELECT_ROWS_OR_FILTER = 'Select rows (Ctrl/Shift+click) or apply a filter';

  // Shared clone/rename/detect core for both the global ("All steps") and per-step "Subset by
  // selection" actions — returns null (after an info toast) when there's nothing to subset.
  // Falls back to the active filter when nothing is explicitly selected, so filtering alone is
  // enough to subset without an extra manual "select all" step — explicit selection still wins
  // when both are present, since it's the more deliberate action.
  function cloneSubsetByRows(df: DG.DataFrame, emptyMsg: string): DG.DataFrame | null {
    const sel = df.selection;
    const mask = sel.trueCount === 0 ? df.filter : sel;
    if (mask.trueCount === 0) {
      grok.shell.info(emptyMsg);
      return null;
    }
    // Also covers the trivial "filter matches every row" case — nothing left to narrow either way.
    if (mask.trueCount === df.rowCount) {
      grok.shell.info('All rows are selected — nothing to subset.');
      return null;
    }
    const subset = df.clone(mask);
    subset.name = `${df.name} (subset, ${subset.rowCount}/${df.rowCount} rows)`;
    detectChemSemTypes(subset);
    // df.clone(mask) leaves the SOURCE df's own .selection set to `mask` as a side effect
    // (confirmed live) — left uncleared, a later filter-only "Subset by selection" on the restored
    // original table picks up this stale selection instead of the new filter, since selection is
    // checked first above.
    df.selection.setAll(false, false);
    // Same clone-carries-BitSet-state issue on the other side: the new subset's OWN .filter comes
    // back all-false (confirmed live) when cloning off a filtered source — reset it so the subset
    // isn't born with every row hidden.
    subset.filter.setAll(true, false);
    return subset;
  }

  // NOTE: assignTableInput's null-then-real two-step assignment does NOT reliably re-render via the
  // input's own onChanged -> onTableChanged reaction alone — live-verified: relying solely on that
  // path left the grid empty after "Subset by selection" (reproduced with filters off too, so it
  // isn't a filters-panel issue). The explicit calls below are intentionally redundant with
  // onTableChanged's own re-render — correctness over the extra render.
  function subsetBySelection(
    input: DG.InputBase<DG.DataFrame | null>, state: SubsetState,
    mountFn: () => void, gridLabel: string, noTableMsg?: string,
  ): void {
    const df = input.value;
    if (!df) {
      if (noTableMsg) grok.shell.info(noTableMsg);
      return;
    }
    const subset = cloneSubsetByRows(df,
      `Select rows in the ${gridLabel} (Ctrl/Shift+click) or apply a filter first.`);
    if (!subset) return;
    // `input` is a ChoiceInput whose `.value` getter re-wraps the underlying table on every read, so
    // `df` is never reference-equal to `state.prev` even when it's the same table — compare by name
    // instead (live-verified via diagnostic logging: same `.name`, but `!==` by reference).
    if (df.name !== state.prev?.name) state.original = df; // remember the user's own table for "Use all"
    const prev = state.prev;
    state.prev = subset;
    assignTableInput(input, subset);
    mountFn();
    // The just-mounted Filters viewer clobbers subset.filter back to all-false when a filtered
    // column has zero variance in the subset (e.g. the very HBD value it was filtered to) — its
    // histogram widgets recompute their range shortly after construction and reset the bitset in the
    // process. Live-verified via staged timing diagnostics: the clobber lands within ~100ms of mount
    // and never recurs after that, so a single reset at 200ms (past that window) sticks permanently —
    // an immediate or same-tick reset does not survive it.
    setTimeout(() => subset.filter.setAll(true, true), FILTER_REMOUNT_SETTLE_MS);
    refreshValidation();
    // Close the previous subset only after the input has switched away from it.
    if (prev && prev.name !== subset.name && prev.name !== df.name)
      try {grok.shell.closeTable(prev);} catch (e) {console.warn(`Could not close prev subset: ${e}`);}
  }

  // Undo "Subset by selection": put the remembered full table back into the input and close the
  // subset clone. No-op (with a hint) when the full library is already in use.
  function restoreFullTable(
    input: DG.InputBase<DG.DataFrame | null>, state: SubsetState, mountFn: () => void, noun: string,
  ): void {
    const orig = state.original;
    if (!orig || input.value?.name === orig.name) {
      // No pending subset swap to undo, but the user may still have an active filter (via the
      // filters panel, with no "Subset by selection" click in between) — "Use all" reads as "show
      // every row", not just "undo the last subset", so clear it too instead of no-op'ing.
      const current = input.value;
      if (current && current.filter.trueCount < current.rowCount) {
        current.filter.setAll(true, true);
        return;
      }
      grok.shell.info(`The full ${noun} library is already in use.`);
      return;
    }
    const prev = state.prev;
    state.prev = null;
    assignTableInput(input, orig);
    mountFn();
    // `orig` can carry a stale filter from before it was swapped out (e.g. the user filtered it,
    // then clicked "Subset by selection" without ever clearing that filter) — "Use all" should mean
    // every row is visible, not "whatever was filtered last time this table was active". Same
    // deferred reset as subsetBySelection: mountFn() just remounted a fresh Filters viewer, whose
    // async histogram-widget setup can clobber an immediate reset (live-verified elsewhere this
    // session), so defer past that window.
    // NOTE: a categorical filter's own checkboxes can still show the pre-reset selection after this
    // (they read a per-column filter memory, not the live bitset) — live-verified that re-mounting
    // to force them to redraw actually re-applies that same remembered selection and re-clobbers the
    // bitset right back, which is worse. Leaving the checkbox display stale is the lesser bug: the
    // grid and row counts are correct, only the filter panel's own checkmarks lag until the user
    // next touches that widget.
    setTimeout(() => orig.filter.setAll(true, true), FILTER_REMOUNT_SETTLE_MS);
    refreshValidation();
    if (prev && prev.name !== orig.name)
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

  // `syncConfigToQuickInputs()` sets values on the quick inputs (rounds, depth-first, etc.), whose
  // own onChanged listeners cascade into a full remount of the Reactions/BBs/Reagents grids — and
  // each freshly-mounted Filters viewer can silently reapply a stale, much-earlier categorical
  // selection instead of the table's actual current filter (live-verified via mount-order tracing:
  // two full remounts fire per call, each seeing the correct filter at mount time, so the reapply
  // happens in the filter widget's own async post-construction step, same root cause fixed elsewhere
  // in this file for subsetBySelection/restoreFullTable). Snapshot each table's filter before the
  // cascade and restore it after, once that async reapply has had time to fire.
  function withPreservedFilters(fn: () => void): void {
    const saved = [templatesInput.value, bbsInput.value, reagentsInput.value]
      .filter((d): d is DG.DataFrame => d != null)
      .map((d) => ({df: d, mask: d.filter.clone()}));
    fn();
    setTimeout(() => {
      for (const {df, mask} of saved) df.filter.copyFrom(mask, true);
    }, FILTER_REMOUNT_SETTLE_MS);
  }

  // ---- Buttons ----
  const editConfigBtn = ui.button('Advanced limits & product filters', async () => {
    syncQuickInputsToConfig();
    const updated = await openConfigDialog(config);
    if (updated) {
      config = updated;
      withPreservedFilters(() => syncConfigToQuickInputs());
      refreshValidation();
    }
  });
  ui.tooltip.bind(editConfigBtn,
    'Caps (components, combinations, routes per compound) and product filters — every config field.');

  const loadYamlBtn = ui.button('Load YAML…', async () => {
    const f = await pickFile('.yaml,.yml');
    if (!f) return;
    try {
      const text = await f.text();
      config = configFromYaml(text);
      withPreservedFilters(() => syncConfigToQuickInputs());
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
  // Disabled buttons get pointer-events:none, so hover/click never reaches them — bind the
  // tooltip to the ribbon-item ancestor instead (must stay its direct, unwrapped child).
  let runBtnRibbonItem: HTMLElement | null = null;
  let lastValidationMsg: string | null = null;
  // WeakSet survives the platform replacing the ribbon-item node mid-session (e.g. after
  // "Subset by selection") — re-attaches on a new node, no-ops on an unchanged one.
  const ribbonItemsWithClickListener = new WeakSet<HTMLElement>();
  const bindRunTooltip = (msg: string | null): void => {
    lastValidationMsg = msg;
    const el = runBtn.closest<HTMLElement>('.d4-ribbon-item');
    if (el) {
      runBtnRibbonItem = el;
      if (!ribbonItemsWithClickListener.has(el)) {
        ribbonItemsWithClickListener.add(el);
        el.addEventListener('click', (e) => {
          if (runBtn.disabled && lastValidationMsg)
            ui.tooltip.show(lastValidationMsg, (e as MouseEvent).clientX, (e as MouseEvent).clientY);
        });
      }
    }
    ui.tooltip.bind(runBtnRibbonItem ?? runBtn,
      msg ?? 'Run enumeration with the current config and add the result to the workspace.');
  };

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
  // No custom wrapping div and no manual gap/margin — passed straight into their own ribbon group
  // (see setRibbonPanels below), same as the Markush Enumerator; the native ribbon group's own flex
  // layout handles spacing between runBtn/cancelBtn/progressLabel.

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
    badge?: TabBadge;
    noTableMsg?: string;
    emptyMsg?: string;
    // How to fold this component's per-round work df into a PerRoundOverride fragment.
    apply: (o: PerRoundOverride, work: DG.DataFrame, cfg: EnumeratorConfig) => void;
  }
  interface DataPanel {
    panel: HTMLElement; render: () => void; onTableChanged: () => void; onRoundsChanged: () => void;
    // Refreshes the bar + override dots only — for a mode/reagents change that can affect whether an
    // override applies (see hasOverride's breadth-first check) without touching the grid itself.
    refreshDisplay: () => void;
    // Applies this component's round-r override (if any) onto `out`; returns whether it did.
    // Lets buildPerRoundOverrides reach each panel's own stepState without a shared indexed array.
    applyOverrideForRound: (r: number, out: PerRoundOverride, cfg: EnumeratorConfig) => boolean;
  }
  function makeDataPanel(o: DataPanelOpts): DataPanel {
    let selStep = 0; // 0 = All steps (full library); 1..rounds = that round's subset
    let filtersOn = false; // funnel toggle: show the standard Datagrok filters panel next to the grid
    let currentDf: DG.DataFrame | null = null;
    // Per-round state (index k-1 = round k) — one array instead of parallel ones, so df/sub/committed
    // can't drift out of sync. `committed` is an explicit flag, NOT inferred from row-count: a step's
    // clone can silently drift from the global table if rows are added/removed in place (no onChanged
    // fires), and inferring "override" from that drift was confirmed to let deleted rows resurrect in
    // a round's reactant pool. Set only by subsetStepBySelection, cleared only by useAllForStep.
    interface StepState { df: DG.DataFrame | null; sub: Subscription | null; committed: boolean; }
    const stepState: StepState[] = [];
    // Step selector: a real ui.tabControl. Each pane owns its own persistent barHost/gridHost (built
    // once via its addPane factory) — relocating a live grid between panes corrupts it, so nothing is
    // shared. `paneHosts[selStep]` gives renderBar/renderGrid the current pane's own hosts.
    const stepTabsHost = ui.div([], {style: {flex: '1 1 0', minHeight: '0', display: 'flex',
      flexDirection: 'column', overflow: 'hidden'}});
    let stepTabsSub: Subscription | null = null;
    const stepDots: (HTMLElement | null)[] = []; // index k = 1..rounds; index 0 ("All steps") unused
    let paneHosts: {barHost: HTMLElement; gridHost: HTMLElement}[] = [];

    // A step "has an override" only once stepState[k-1].committed is explicitly set — see the flag's
    // own comment above for why row-count inference was wrong.
    const hasOverride = (k: number): boolean => {
      if (!stepState[k - 1]?.committed) return false;
      // Breadth-first ignores per-round BB overrides entirely (a round draws from all earlier
      // products, see `eligibleSmiles`) — don't show the dot as if it did something. Reactions/
      // reagents overrides are unaffected.
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
      const entry = stepState[k - 1];
      if (entry) { entry.df = work; entry.sub = sub; }
      else stepState[k - 1] = {df: work, sub, committed: false};
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
        `${SELECT_ROWS_OR_FILTER} to use only those ${o.noun} in step ${k}.`);
      if (!subset) return;
      setStepWork(k, subset);
      stepState[k - 1].committed = true;
      renderGrid(); renderBar(); updateDots();
      // See subsetBySelection's identical fix and comment on the 200ms delay.
      setTimeout(() => subset.filter.setAll(true, true), FILTER_REMOUNT_SETTLE_MS);
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
      // Close each pane's mounted viewer(s) BEFORE wiping stepTabsHost — otherwise their gridHost
      // divs are dropped from the DOM while still registered in `mountedViewers`, orphaning the
      // Viewer instances (never closed) instead of releasing their Dart-side resources.
      for (const ph of paneHosts) if (ph) closeMountedViewers(ph.gridHost);
      stepTabsHost.innerHTML = '';
      stepDots.length = 0;
      paneHosts = [];
      const tc = ui.tabControl(null, false);
      // tc.root must fill available height, not size to its header-strip content — each pane's own
      // content div (built below) needs the real space to lay its grid out in.
      tc.root.style.width = '100%';
      tc.root.style.flex = '1 1 0';
      tc.root.style.minHeight = '0';
      tc.root.style.overflow = 'hidden';
      // Builds one pane's persistent content (barHost+gridHost), recorded in paneHosts at its OWN
      // fixed index k (not push-order) — works whether TabControl's addPane factory runs eagerly or
      // lazily, since position k always lands at paneHosts[k] regardless of firing order.
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
        pane.header.style.position = 'relative';
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
      // Each rebuild unsubscribes the prior instance, but the LAST one has no later rebuild to retire
      // it — track it in view.subs too so view close still reaches it.
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
      filterIcon.style.cursor = 'pointer';
      filterIcon.style.padding = '2px 5px';
      filterIcon.style.flex = '0 0 auto';
      filterIcon.style.color = filtersOn ? 'var(--blue-2)' : 'var(--grey-5)';
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
        ui.tooltip.bind(btn, `Replace the ${o.noun} library with only the selected rows, or — if nothing is ` +
          `selected — the currently filtered rows (applies to every step). Click "Use all" to restore the ` +
          `full set.`);
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
        ui.tooltip.bind(btn, `Narrow step ${selStep} to only the selected rows (Ctrl/Shift+click), or — if ` +
          `nothing is selected — the currently filtered rows. Click "Use all" to go back to the full ` +
          `${o.noun} library.`);
        const useAll = ui.button('Use all', () => useAllForStep(selStep));
        ui.tooltip.bind(useAll, `Undo "Subset by selection" so step ${selStep} uses the full ${o.noun} library ` +
          `(same as "All steps").`);
        barHost.append(
          hintEl(`${SELECT_ROWS_OR_FILTER}, then "Subset by selection" to use only those ${o.noun} in step ` +
            `${selStep}.`),
          status, filterIcon, btn, useAll);
      }
    }

    function renderGrid(): void {
      const gridHost = paneHosts[selStep]?.gridHost;
      if (!gridHost) return; // pane not built yet (shouldn't happen once buildStepTabs has run once)
      currentDf = selStep === 0 ? o.input.value : stepClone(selStep);
      if (!currentDf) {
        // Close mounted viewers BEFORE wiping the DOM — closing after innerHTML='' hands the viewer
        // a detached container, which throws ("Cannot read properties of null") deep in the Dart-side
        // close path. Under rapid re-triggering (e.g. the filter icon clicked several times in quick
        // succession) that cascades badly enough to crash the tab's renderer.
        closeMountedViewers(gridHost);
        gridHost.innerHTML = '';
        gridHost.appendChild(ui.divText(o.emptyMsg ?? `No ${o.noun} table selected.`,
          {style: {color: 'var(--grey-5)', padding: '20px', textAlign: 'center'}}));
        if (selStep === 0) o.badge?.refresh(null);
        return;
      }
      mountDf(gridHost, currentDf, filtersOn); // mountDf itself closes-then-clears the host
      if (selStep === 0) o.badge?.refresh(currentDf.rowCount);
    }

    function render(): void { renderGrid(); renderBar(); updateDots(); }
    function refreshDisplay(): void { renderBar(); updateDots(); }

    const onTableChanged = (): void => {
      if (o.input.value) detectChemSemTypes(o.input.value);
      // `state.original` is deliberately NOT touched here. onChanged's exact timing relative to
      // subsetBySelection's own synchronous assignment isn't reliable enough to gate a reset on
      // (two earlier attempts — a suppress flag, then a direct o.state.prev comparison here — were
      // both live-verified to occasionally wipe `original` mid-chain, making "Use all" undo only the
      // last "Subset by selection" instead of the whole chain back to the true original). Ownership
      // of `original`'s lifecycle lives entirely in subsetBySelection/restoreFullTable's own
      // synchronous, user-triggered logic instead — see subsetBySelection's `df.name !== state.prev?.name`
      // check, which self-corrects even after an unrelated file swap. Trade-off: if the user swaps
      // to a different file WITHOUT ever subsetting it, then clicks "Use all", it restores the
      // previous file's original instead of a no-op — an edge case outside "Use all"'s intended use.
      for (const s of stepState) s?.sub?.unsubscribe();
      stepState.length = 0;
      buildStepTabs(0);
      refreshValidation();
    };
    // Deliberately does NOT drop stepState entries beyond the new round count — Dart int inputs fire
    // onChanged per keystroke, so typing "10" over "5" transiently passes through "1", which would
    // eagerly destroy committed overrides on steps 2-5 before the user finishes typing.
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
    return {panel, render, onTableChanged, onRoundsChanged, refreshDisplay, applyOverrideForRound};
  }

  const templatesCtl = makeDataPanel({idx: 0, noun: 'reaction templates',
    input: templatesInput, state: templatesState,
    apply: (o, work, cfg) => { o.templates = extractTemplates(cfg, work); }});
  const bbsCtl = makeDataPanel({idx: 1, noun: 'building blocks',
    input: bbsInput, state: bbsState,
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

  // Rebuilds step strips on round-count change or when a component's "All steps" table changes —
  // narrowing "All steps" intentionally discards every step's committed override.
  //
  // `max` on an int input only shows a tooltip, it doesn't clamp — an over-max value is instead
  // caught by validate(); roundCount() separately caps at MAX_ROUNDS so tab count can't blow up.
  view.subs.push(numRoundsInput.onChanged.subscribe(() => {
    refreshValidation();
    dataCtls.forEach((c) => c.onRoundsChanged());
  }));
  view.subs.push(templatesInput.onChanged.subscribe(() => templatesCtl.onTableChanged()));
  view.subs.push(bbsInput.onChanged.subscribe(() => bbsCtl.onTableChanged()));
  // BB override dot/status is mode-aware (hasOverride hides it in breadth-first) — refresh the bar +
  // dots on any mode switch so they don't show stale state. Neither the underlying table nor the
  // grid changes here, so refreshDisplay (not render) is enough — no grid rebuild needed.
  view.subs.push(depthFirstInput.onChanged.subscribe(() => dataCtls.forEach((c) => c.refreshDisplay())));
  view.subs.push(reagentsInput.onChanged.subscribe(() => {
    reagentsCtl.onTableChanged(); // its own table changed — needs the full rebuild
    templatesCtl.refreshDisplay();
    bbsCtl.refreshDisplay();
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
  // Reagents row-count badge — the only data tab without a count shown elsewhere (reactions/BBs
  // already have it in the ribbon chips and accordion pane subtitles).
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

  // Enumerate gets its own ribbon group, matching the Markush Enumerator (pt-chem-enum-dialog.ts) —
  // its `appActionHost` (the run button) is also its own group, never bundled with other items into
  // one custom flex div. A run button sharing a group with unrelated items is what picks up the
  // platform's own group-level background/shadow oddly.
  view.setRibbonPanels([
    [appInfoIcon],
    [runBtn, cancelBtn, progressLabel],
    [chipReactions, ribbonArrow, chipBbs, chipCombine, cfgEstEl],
  ]);
  view.append(root);

  // Strips the platform's default ribbon-group shadow/background (Markush Enumerator does the
  // same). Walk up from `runBtn`, not down from `.d4-root` — that's the grid viewer's own
  // wrapper class, not a ribbon ancestor. Style gets wiped by platform startup, so reassert on
  // every tick rather than stopping at first find; click-listener attach lives in bindRunTooltip.
  function applyRibbonFixup(attempt = 0): void {
    const el = runBtn.closest<HTMLElement>('.d4-ribbon-item');
    if (el) {
      document.querySelectorAll<HTMLElement>('.d4-ribbon-group, .d4-ribbon-item').forEach((g) => {
        g.style.background = 'transparent';
        g.style.boxShadow = 'none';
        g.style.border = 'none';
      });
    }
    if (attempt < 200) setTimeout(() => applyRibbonFixup(attempt + 1), 50);
  }
  applyRibbonFixup();
  bindRunTooltip(validate());

  // Looks redundant with buildStepTabs(0)'s own render, but isn't: removing it left the
  // initial grid rendered into a not-yet-sized host and empty on some loads.
  templatesCtl.render();
  bbsCtl.render();
  reagentsCtl.render();
  refreshValidation();
  return view;
}
