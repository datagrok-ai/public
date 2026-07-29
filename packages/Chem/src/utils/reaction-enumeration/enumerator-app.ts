/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import {Subscription} from 'rxjs';
import {debounceTime} from 'rxjs/operators';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package, PackageFunctions} from '../../package';
import {cloneConfig, configFromYaml, configToYaml, DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {buildCombinationLimitFields, buildProductFilterFields} from './config-form';
import {getRdKitModule} from '../chem-common-rdkit';
import {enumerate, EnumerationProgress, OutputRow, PerRoundOverride, TemplateInput, tryGetRxn} from './enumerate';
import {MountedViewerRegistry} from './viewer-mount';
import {StrategySummary} from './strategy-summary';
import {PreviewPanel} from './preview-panel';
import {
  buildPerRoundOverrides as buildPerRoundOverridesForPanels, DataPanel, DataPanelDeps, makeTabBadge,
  overrideCountFor,
} from './data-panel';
import {RunControls} from './run-controls';

const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
const BUNDLED_BBS = 'enumerations/bb.csv';
const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

// Shared "custom subset" indicator color — round tabs' dot and the Strategy summary's dot.
export const OVERRIDE_DOT_COLOR = 'var(--orange-2, #c98a1b)';
// Shared look for the small "changed/custom" dots; call sites add their own display mode and spacing.
export const CHANGED_DOT_STYLE = {width: '6px', height: '6px', borderRadius: '50%', background: OVERRIDE_DOT_COLOR};

// Shared mode label/rounds text — used by the ribbon chip, the Strategy summary, and the Preview
// recap so they can't drift out of sync (e.g. one saying "2 rounds", another "1 rounds").
export const MODE_LABEL = {depth: 'Depth-first', breadth: 'Breadth-first', reagents: 'Reagents'} as const;
export const roundsLabel = (n: number): string => `${n} round${n === 1 ? '' : 's'}`;

// Shared "differs from platform defaults" checks — drive both the Combination limits/Product
// filters toggle dots and the Strategy summary's "changed from defaults" caveat.
export function combinationLimitsChanged(cfg: EnumeratorConfig): boolean {
  return cfg.max_num_combinations_per_template !== DEFAULT_CONFIG.max_num_combinations_per_template ||
    cfg.keep_building_blocks_in_final_output !== DEFAULT_CONFIG.keep_building_blocks_in_final_output;
}

// Shared panel chrome — a hint/status header bar plus a content host — used by every right-pane
// tab (data grids, Strategy summary, Preview).
export const panelHeader = (hint: string, subsetBtn?: HTMLElement, status?: HTMLElement): HTMLElement => {
  // flex:0 0 auto (not 1 1 auto) so hint and status sit side by side, both left-aligned, instead
  // of hint growing to push status to the far right of the row.
  const hintEl = ui.divText(hint, {style: {
    fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto', marginRight: '4px',
  }});
  const children: HTMLElement[] = [hintEl];
  if (status) children.push(status);
  if (subsetBtn) children.push(subsetBtn);
  return ui.div(children, {style: {
    display: 'flex', alignItems: 'center', gap: '8px', flex: '0 0 auto',
    padding: '4px 8px 5px', borderBottom: '1px solid var(--grey-2)',
  }});
};

// `scrollable` is for plain content hosts (e.g. the Strategy summary card) that can outgrow a
// short window — unlike a grid, which manages its own internal scroll, so its host stays
// overflow:hidden with the bottom fade. min-height:0 is unconditional: without it this wrapper
// (a flex child of the platform's own .d4-tab-content, itself fixed up to min-height:0 in
// chem.css) refuses to shrink below its content's natural height, so a narrow window makes
// .d4-tab-content overflow and get silently clipped further up instead of ever asking this pane
// — or gridHost's own overflow below — to actually engage.
export const tabPanel = (header: HTMLElement, gridHost: HTMLElement, scrollable = false): HTMLElement => {
  // display:flex on the host turns the grid's inline-flex outer display into a block-level
  // flex item, eliminating the 12px baseline-alignment gap that block+inline-flex produces.
  gridHost.style.display = 'flex';
  gridHost.style.flexDirection = 'column';
  gridHost.style.flex = '1 1 0';
  gridHost.style.minHeight = '0';
  if (scrollable) {
    // Vertical only, matching the left "How to combine" pane — fully native scroll.
    gridHost.style.overflowY = 'auto';
    gridHost.style.overflowX = 'hidden';
  } else {
    gridHost.style.position = 'relative';
    gridHost.style.overflow = 'hidden';
    const fade = ui.div([], {style: {
      position: 'absolute', bottom: '0', left: '0', right: '0', height: '48px',
      background: 'linear-gradient(to bottom,transparent,var(--white))', pointerEvents: 'none', zIndex: '1',
    }});
    gridHost.appendChild(fade);
  }
  return ui.div([header, gridHost], {style: {
    height: '100%', display: 'flex', flexDirection: 'column', minHeight: '0',
    background: 'var(--white)', boxSizing: 'border-box', overflow: 'hidden',
  }});
};

export function productFiltersChangedCount(cfg: EnumeratorConfig): number {
  const ps = cfg.products_specs;
  const dps = DEFAULT_CONFIG.products_specs;
  return [
    ps.max_num_heavy_atoms !== dps.max_num_heavy_atoms,
    ps.min_num_carbon_atoms !== dps.min_num_carbon_atoms,
    ps.max_num_carbon_atoms !== dps.max_num_carbon_atoms,
    ps.max_num_hetero_atoms !== dps.max_num_hetero_atoms,
    ps.max_num_nitrogen !== dps.max_num_nitrogen,
    ps.max_num_sulfur !== dps.max_num_sulfur,
    ps.max_num_oxygen !== dps.max_num_oxygen,
    ps.max_num_metals !== dps.max_num_metals,
    ps.max_num_halogens !== dps.max_num_halogens,
    ps.max_num_aromatic_atoms !== dps.max_num_aromatic_atoms,
    ps.max_num_unsaturated_nonaromatic_bonds !== dps.max_num_unsaturated_nonaromatic_bonds,
    ps.only_these_atoms_allowed.join(',') !== dps.only_these_atoms_allowed.join(','),
    ps.remove_radicals !== dps.remove_radicals,
    ps.remove_isotope_information !== dps.remove_isotope_information,
    ps.remove_charged_species !== dps.remove_charged_species,
  ].filter(Boolean).length;
}

// Dart int inputs fire onChanged per keystroke — debounce the expensive step-tab rebuild.
const ROUNDS_INPUT_DEBOUNCE_MS = 300;
// A round tab is built for every round, and product counts grow combinatorially with each one.
export const MAX_ROUNDS = 10;

// Sniff string columns and set semType so the grid renders reactions and molecules: presence of
// `>>` in sampled values wins as ChemicalReaction, else auto-detection handles Molecule etc.
export function detectChemSemTypes(df: DG.DataFrame): void {
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
    let done = false;
    const cleanup = (f: File | null): void => {
      if (done) return;
      done = true;
      window.removeEventListener('focus', onFocus);
      input.remove();
      resolve(f);
    };
    // Cancelling the OS dialog never fires 'change' — only returns focus to the window. Give
    // 'change' a moment to win the race if a file WAS picked, then resolve null otherwise.
    const onFocus = (): void => {
      setTimeout(() => cleanup(input.files?.[0] ?? null), 300);
    };
    window.addEventListener('focus', onFocus);
    input.onchange = () => cleanup(input.files?.[0] ?? null);
    // Newer browsers fire 'cancel' directly on an actual dismiss, no timing guesswork needed.
    // Older browsers simply never dispatch it, leaving the focus fallback above as the only path.
    input.addEventListener('cancel', () => cleanup(null));
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

export function buildResultDataFrame(rows: OutputRow[], name = 'Enumeration result'): DG.DataFrame {
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

export interface BuiltInputs {
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
  // `nullable` must apply even when `table` is already set, or the platform defaults to required.
  const opts: any = {filter, nullable};
  if (table) {
    opts.table = table;
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

export function buildInputs(
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
  view.name = 'Reaction Enumerator';
  view.box = true;

  let config: EnumeratorConfig = cloneConfig(DEFAULT_CONFIG);

  // Fields are built eagerly (syncQuickInputsToConfig needs .inputs before the accordion exists);
  // their ui.form() wrapper is built lazily instead (see lazyFilterForm below).
  let combinationLimitFields = buildCombinationLimitFields(config);
  let productFilterFields = buildProductFilterFields(config);

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

  const blockingColInput = makeColInput('Blocking SMARTS (optional)', templatesDf,
    config.enumeration.reactant_blocking_groups_per_template_column, isStringCol,
    'Optional column whose values are SMARTS patterns (separated by ";" or "|"). Excludes building ' +
    'blocks with functional groups incompatible with this template — a building block matching any ' +
    'of them is skipped for this template only and stays available to all others.', true);

  const rxnNameColInput = makeColInput('Reaction name (optional)', templatesDf,
    config.enumeration.reaction_name_col, isStringCol,
    'Optional column with a friendly name for each reaction template. Surfaces in the output ' +
    '"reaction_name" column.', true);

  const bbsInput = ui.input.table('Building blocks file', {
    value: bbsDf ?? undefined, nullable: false,
    tooltipText: 'Table with the building-block library (SMILES). Pick a workspace table or upload a CSV.',
  });
  const bbColInput = makeColInput('SMILES column', bbsDf, config.enumeration.bb_smiles_column, isStringCol,
    'Column in the building blocks file that contains SMILES.', false);

  const reagentsInput = ui.input.table('Reagents file (optional)', {
    value: undefined, nullable: true,
    tooltipText: 'Optional table of reagent SMILES. When set, switches to reagents mode: every ' +
      'round uses exactly one building block (or product of an earlier round) and fills every ' +
      'remaining slot with reagents from this file — produces derivatives of each BB across rounds.',
  });
  const REAGENTS_COL_TOOLTIP = 'Column in the reagents file that contains the reagent SMILES.';
  const reagentsColInput = makeColInput('Reagent SMILES column', null,
    config.enumeration.reagent_smiles_column, isStringCol, REAGENTS_COL_TOOLTIP, true);
  // The platform marks this invalid because it has no table to pick a column from yet, not because
  // it's empty (nullable is true) — disable it and say why, instead of fighting the invalid style.
  const syncReagentsColEnabled = (): void => {
    const hasTable = reagentsInput.value != null;
    reagentsColInput.enabled = hasTable;
    reagentsColInput.setTooltip(hasTable ? REAGENTS_COL_TOOLTIP : 'Select a reagents file first.');
  };
  syncReagentsColEnabled();

  const exclusionInput = ui.input.table('Exclusion SMARTS file (optional)', {
    value: exclusionDf ?? undefined, nullable: true,
    tooltipText: 'Optional table of SMARTS patterns. Any product matching one of these is rejected.',
  });
  const exclusionColInput = makeColInput('Exclusion SMARTS column', exclusionDf,
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
    reagentsInput.onChanged.subscribe(syncReagentsColEnabled),
  );

  // ---- CONFIG inputs ----
  const numRoundsInput = ui.input.int('Number of rounds',
    {value: config.enumeration.num_rounds, min: 1, max: MAX_ROUNDS, showPlusMinus: true});
  numRoundsInput.setTooltip(
    'Number of consecutive enumeration rounds. Round 1 reacts BBs only; round 2 takes round-1 ' +
    `products and (in depth-first mode) reacts each one with original BBs. Increase for deeper ` +
    `libraries (capped at ${MAX_ROUNDS} — a round tab is built for every round, and product counts ` +
    `grow combinatorially with each one).`);
  // `min`/`max` above only affect the tooltip/spinner, not validation — add that separately.
  numRoundsInput.addValidator((v) => {
    const n = Number(v);
    if (!Number.isFinite(n) || n < 1) return 'Must be at least 1.';
    if (n > MAX_ROUNDS) return `Must be at most ${MAX_ROUNDS} — the round strip shows the first ${MAX_ROUNDS} rounds only.`;
    return null;
  });

  const depthFirstInput = ui.input.bool('Depth first', {value: config.enumeration.depth_first});
  depthFirstInput.setTooltip(
    'When checked, each round r > 1 must combine EXACTLY ONE round-(r-1) product with original ' +
    'BBs (linear chain extension, no merging two complex products). Off (breadth-first) allows any ' +
    'combination from rounds 0..r-1 — typically explodes the search space and produces convergent routes.');

  // Promoted out of "Advanced limits & product filters" — used often enough to live at top level.
  const maxComponentsInput = ui.input.int('Max # components',
    {value: config.max_num_components, min: 1, showPlusMinus: true});
  maxComponentsInput.setTooltip('Max number of reactant components a template may have.');
  // -1 is the config's own "no cap" sentinel — showing it as a literal number reads as a developer
  // detail, not a value a user would type. Blank means the same thing and is shown/read as such.
  const maxRoutesInput = ui.input.int('Max routes per compound', {
    value: config.max_num_routes_per_compound < 0 ? undefined : config.max_num_routes_per_compound,
    nullable: true, showPlusMinus: true,
  });
  maxRoutesInput.setTooltip('Cap on the number of routes saved per product. Leave blank for no cap.');
  // Decrementing the stepper past 0 floors through to -1 instead of blank.
  maxRoutesInput.onChanged.subscribe(() => {
    if (maxRoutesInput.value === -1) maxRoutesInput.value = null;
  });

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
        <li><b>Reagents file</b> (optional) — switches to <i>reagents mode</i>: every round uses exactly one BB or earlier-round product, with reagents in the remaining slots. Yields derivatives of each BB across rounds.</li>
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
  const MODE_ABBR = {depth: 'DF', breadth: 'BF', reagents: 'RM'} as const; // ribbon chip only, MODE_LABEL elsewhere
  // The raw round count as currently displayed/edited (not the defensively-clamped one DataPanel's
  // roundCount() uses for building round tabs) — shared by the ribbon chip, the Strategy summary, and
  // the Preview recap so they can't drift out of sync (e.g. one saying "2 rounds", another "1 rounds").
  function currentRounds(): number {
    return numRoundsInput.value ?? config.enumeration.num_rounds;
  }

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
    // No-op while syncConfigToQuickInputs() is pushing config -> inputs: each setAndFire() fires
    // onChanged, which reaches here via refreshValidation()/validate() before every input has been
    // updated — reading them now would write stale (not-yet-pushed) values back into config.
    if (pushingConfigToInputs) return;
    config.enumeration.num_rounds = numRoundsInput.value ?? config.enumeration.num_rounds;
    config.enumeration.depth_first = !!depthFirstInput.value;
    config.max_num_components = maxComponentsInput.value ?? config.max_num_components;
    config.max_num_routes_per_compound = maxRoutesInput.value ?? -1; // blank == no cap
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
    combinationLimitFields.syncToConfig(config);
    productFilterFields.syncToConfig(config);
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
      setAndFire(maxRoutesInput, config.max_num_routes_per_compound < 0 ? null : config.max_num_routes_per_compound);
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
      // Neither field group has a "set value" hook — rebuild from the loaded config; the rebuilt
      // inputs are brand new objects, so they need their own revalidation wiring too.
      combinationLimitFields = buildCombinationLimitFields(config);
      combinationLimitsForm.invalidate();
      combinationLimitFields.inputs.forEach((inp) => wireValidation(inp));
      productFilterFields = buildProductFilterFields(config);
      productFilterForm.invalidate();
      productFilterFields.inputs.forEach((inp) => wireValidation(inp));
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
    // Sync before reading config below — validate() syncs too, but only after refreshCfgRibbon(),
    // which would otherwise read one refresh behind.
    syncQuickInputsToConfig();
    refreshCfgRibbon();
    refreshStrategyCards();
    const err = validate();
    validationDiv.textContent = err ?? '';
    runControls.setValidation(err);
  }
  syncQuickInputsToConfig();

  // Tab row-count badge. Reactions/BBs already show their row count via the always-visible ribbon
  // chips (chipReactionsC/chipBbsC) and the accordion pane subtitles — a tab badge there would just
  // repeat the same number a third time. Reagents has neither, so it keeps a badge.
  const reagentsBadge = makeTabBadge();

  const viewerHost = new MountedViewerRegistry(view);

  // Re-validate on every input change so the Run button stays accurate.
  const wireValidation = (input: DG.InputBase<unknown>): void => {
    view.subs.push(input.onChanged.subscribe(() => {
      refreshValidation();
    }));
  };
  [smartsColInput, blockingColInput, rxnNameColInput, bbColInput, reagentsColInput,
    exclusionInput, exclusionColInput, numRoundsInput, depthFirstInput,
    maxComponentsInput, maxRoutesInput,
    ...combinationLimitFields.inputs, ...productFilterFields.inputs,
  ].forEach((inp) => wireValidation(inp));

  // ---- Buttons ----
  // Icon buttons in the ribbon: 'folder-open' for import, 'arrow-to-bottom' for export.
  const loadYamlBtn = ui.iconFA('folder-open', async () => {
    const f = await pickFile('.yaml,.yml');
    if (!f) return;
    try {
      const text = await f.text();
      config = configFromYaml(text);
      viewerHost.withPreservedFilters(
        [templatesInput.value, bbsInput.value, reagentsInput.value].filter((d): d is DG.DataFrame => d != null),
        () => syncConfigToQuickInputs());
      refreshValidation();
      grok.shell.info(`Loaded config from ${f.name}.`);
    } catch (e) {
      grok.shell.error(`Could not load YAML: ${e instanceof Error ? e.message : String(e)}`);
    }
  }, 'Load a YAML config file from disk and apply it to the form.');

  const saveYamlBtn = ui.iconFA('arrow-to-bottom', () => {
    syncQuickInputsToConfig();
    DG.Utils.download('enumerator-config.yaml', configToYaml(config), 'text/yaml');
  }, 'Download the current config as a YAML file.');

  // ---- Run / Cancel ----
  const runControls = new RunControls({
    getConfig: () => config,
    templatesInput, bbsInput, reagentsInput, exclusionInput,
    validate, syncQuickInputsToConfig, buildPerRoundOverrides, refreshValidation,
  });

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
  // Reagents mode is auto-derived (set via a file in Extras), never manually selectable here.
  const reagentsModeNote = ui.divH([ui.iconFA('info-circle'),
    ui.span([' Reagents mode active — set via a reagents file in Extras.'])],
  {style: {fontSize: '11px', color: 'var(--grey-6)', gap: '6px', alignItems: 'center',
    padding: '6px 10px', display: 'none'}});
  ui.tooltip.bind(reagentsModeNote, 'Every round uses exactly one building block or earlier-round product, ' +
    'with reagents filling the remaining slots. Automatically selected as soon as a reagents file is loaded ' +
    'in the Extras tab — to go back to Depth-first/Breadth-first, clear the reagents file there.');

  function refreshStrategyCards(): void {
    const cur = currentMode();
    const hasReagents = cur === 'reagents';
    applyStratCardStyle(stratDepthCard, 'depth', cur, !hasReagents);
    applyStratCardStyle(stratBreadthCard, 'breadth', cur, !hasReagents);
    reagentsModeNote.style.display = hasReagents ? '' : 'none';
  }
  stratDepthCard.root.onclick = (): void => {
    if (reagentsInput.value != null) return;
    setAndFire(depthFirstInput, true);
  };
  stratBreadthCard.root.onclick = (): void => {
    if (reagentsInput.value != null) return;
    setAndFire(depthFirstInput, false);
  };

  // Right-pane tab references — assigned when tabs are built; used by section-open handlers for
  // context-sensitive tab switching. Declared here so openAccPaneAndSyncTab can close over them.
  let templatesPane: DG.TabPane | undefined;
  let bbsPane: DG.TabPane | undefined;
  let reagentsPane: DG.TabPane | undefined;

  // Target pane resolved lazily via a thunk: Reactions pane is added expanded, so its factory runs
  // synchronously inside addPane, before later panes exist — capturing directly would hit the TDZ.
  const mkNavBtn = (kind: 'next' | 'back', getTarget: () => DG.AccordionPane, label: string): HTMLElement => {
    const btn = ui.button(`${kind === 'next' ? 'Next' : 'Back'}: ${label}`, () => openAccPaneAndSyncTab(getTarget()));
    btn.classList.add(`chem-enum-${kind}-btn`);
    return btn;
  };
  const mkNextBtn = (getTarget: () => DG.AccordionPane, label: string): HTMLElement =>
    mkNavBtn('next', getTarget, label);
  const mkBackBtn = (getTarget: () => DG.AccordionPane, label: string): HTMLElement =>
    mkNavBtn('back', getTarget, label);
  // Back (if any) on the left, Next/action (if any) on the right — one consistent row per pane.
  const navRow = (back: HTMLElement | null, next: HTMLElement | null): HTMLElement =>
    ui.divH([back ?? ui.div([]), next ?? ui.div([])], {classes: 'chem-enum-nav-row'});

  const accordion = ui.accordion();
  accordion.root.classList.add('chem-enum-accordion');
  // allowDragOut (5th arg) defaults to true; panes shouldn't be draggable out of this panel.
  // One shared form so all four fields' labels align (two forms would size independently).
  const accReactionsPane = accordion.addPane('Reactions', () =>
    ui.divV([ui.form([templatesInput, smartsColInput, blockingColInput, rxnNameColInput]),
      navRow(mkBackBtn(() => accCombinePane, 'How to combine'), mkNextBtn(() => accBbsPane, 'Building blocks'))]),
  true, null, false);
  const accBbsPane = accordion.addPane('Building blocks', () =>
    ui.divV([ui.form([bbsInput, bbColInput]),
      navRow(mkBackBtn(() => accReactionsPane, 'Reactions'), mkNextBtn(() => accExtrasPane, 'Extras'))]),
  false, null, false);
  const extrasForm = ui.form([reagentsInput, reagentsColInput, exclusionInput, exclusionColInput]);
  const accExtrasPane = accordion.addPane('Extras', () =>
    ui.divV([extrasForm,
      navRow(mkBackBtn(() => accBbsPane, 'Building blocks'), mkNextBtn(() => accPreviewPane, 'Preview'))]),
  false, null, false);
  // Flags "differs from platform defaults" without expanding the toggle.
  const mkChangedDot = (tooltip: string): HTMLElement => {
    const dot = ui.div([], {style: {...CHANGED_DOT_STYLE, display: 'none'}});
    ui.tooltip.bind(dot, tooltip);
    return dot;
  };
  // Attaches a changed-dot to a pane's own header, spaced off the label text.
  const attachChangedDot = (pane: DG.AccordionPane, tooltip: string): HTMLElement => {
    const dot = mkChangedDot(tooltip);
    dot.style.marginLeft = '6px';
    pane.root.querySelector('.d4-accordion-pane-header')?.appendChild(dot);
    return dot;
  };
  // Cancels the nested accordion's own chevron indent so its rows start flush with the main fields —
  // measured against the real gap instead of a hand-picked constant, and re-measured on every real
  // size change (not just once at mount) since a sibling reflow can still be in flight when the pane
  // first expands. Observes a stable ancestor, not el itself, since el's own margin-left change would
  // otherwise perturb its resolved width and could re-trigger a ResizeObserver watching el directly.
  const flushIndent = (el: HTMLElement, reference: HTMLElement): Subscription => {
    const apply = (): void => {
      const extra = el.getBoundingClientRect().left - reference.getBoundingClientRect().left;
      const current = parseFloat(getComputedStyle(el).marginLeft) || 0;
      const next = current - extra;
      if (next === current) return;
      if (next === 0) el.style.removeProperty('margin-left');
      else el.style.setProperty('margin-left', `${next}px`, 'important');
    };
    const sub = new Subscription();
    sub.add(ui.onSizeChanged(reference).subscribe(apply));
    sub.add(ui.onSizeChanged(accCombinePane.root).subscribe(apply));
    setTimeout(apply, 0);
    view.subs.push(sub);
    return sub;
  };
  // Builds ui.form() lazily, only once the pane's content factory first runs — building it while
  // still detached would measure it at 0 width and get it marked .ui-form-condensed regardless of
  // label width. invalidate() rebuilds in place after a config reload swaps in fresh inputs — the
  // stale build's flushIndent subscriptions are disposed first, or they'd keep firing against a
  // detached form on every future resize.
  const lazyFilterForm = (getInputs: () => DG.InputBase<unknown>[]):
  {getRoot: () => HTMLElement; invalidate: () => void} => {
    let root: HTMLElement | null = null;
    let indentSub: Subscription | null = null;
    const build = (): HTMLElement => {
      const form = ui.form(getInputs());
      indentSub = flushIndent(form, numRoundsInput.root);
      return form;
    };
    return {
      getRoot: (): HTMLElement => root ??= build(),
      invalidate: (): void => {
        if (!root) return;
        indentSub?.unsubscribe();
        const fresh = build();
        root.replaceWith(fresh);
        root = fresh;
      },
    };
  };
  // Shared label-column width across all three forms in this section (rounds, combination limits,
  // product filters) — so Product filters renders with the exact same column Combination limits
  // does, not its own wider one. Measured from the actual widest caption across all three forms
  // (via a hidden offscreen probe in the label's own font) rather than a hardcoded constant, so it
  // never goes stale if a caption is added, removed, or reworded later. Set via a CSS custom
  // property + stylesheet !important rule, not a plain inline style: the platform's own per-form
  // label auto-sizing runs asynchronously after mount and would silently overwrite a plain inline
  // style outright.
  const measureLabelWidths = (labels: HTMLElement[]): number[] => {
    const probe = document.createElement('span');
    probe.style.position = 'absolute';
    probe.style.visibility = 'hidden';
    probe.style.whiteSpace = 'nowrap';
    probe.style.left = '-9999px';
    document.body.appendChild(probe);
    const widths = labels.map((label) => {
      probe.style.font = getComputedStyle(label).font;
      probe.textContent = label.textContent;
      return probe.getBoundingClientRect().width;
    });
    probe.remove();
    return widths;
  };
  const sectionInputs = [
    numRoundsInput, maxComponentsInput, maxRoutesInput,
    ...combinationLimitFields.inputs, ...productFilterFields.inputs,
  ];
  const sharedLabelWidth = Math.ceil(Math.max(...measureLabelWidths(sectionInputs.map((inp) => inp.captionLabel))));
  // Editors need pinning too, not just labels: the platform widens an editor from its normal ~150px
  // to ~176px on whichever of the three forms it currently considers .ui-form-condensed, so without
  // this the Product filters column can end up visibly wider than Combination limits' even with
  // labels aligned. 150px matches the platform's own un-condensed default.
  const sharedEditorWidth = 150;
  // Mirrors the per-input-type minimum widths the platform itself uses to decide .ui-form-condensed
  // (js-api ui.ts's getInputsMinWidths — text/table inputs need 200px, float/date 140px, everything
  // else 100px). Sizing the form to fit the widest one up front, the same way the platform already
  // does for its own dialog forms (a min-width computed from label + input minimums instead of
  // reacting to condensed after the fact), means these forms never need condensed layout at all.
  const platformInputMinWidth = (input: DG.InputBase<unknown>): number => {
    const el = input.root;
    if (el.classList.contains('ui-input-text') || el.classList.contains('ui-input-table')) return 200;
    if (el.classList.contains('ui-input-float') || el.classList.contains('ui-input-date')) return 140;
    return 100;
  };
  const formMinInputWidth = Math.ceil(Math.max(...sectionInputs.map(platformInputMinWidth)));
  // Independently-collapsible sub-sections within "How to combine" (no forced exclusivity, unlike
  // the outer wizard-navigation accordion).
  const limitsAccordion = ui.accordion();
  const combinationLimitsForm = lazyFilterForm(() => combinationLimitFields.inputs);
  const productFilterForm = lazyFilterForm(() => productFilterFields.inputs);
  const combinationLimitsPane = limitsAccordion.addPane('Combination limits (optional)',
    () => combinationLimitsForm.getRoot(), false, null, false);
  const productFilterPane = limitsAccordion.addPane('Product filters (optional)',
    () => productFilterForm.getRoot(), false, null, false);
  const combinationLimitsDot = attachChangedDot(combinationLimitsPane, 'Changed from platform defaults.');
  const productFiltersDot = attachChangedDot(productFilterPane, 'Changed from platform defaults.');
  const accCombinePane = accordion.addPane('How to combine', () => ui.divV([
    ui.divH([
      ui.divText('Strategy', {style: {fontSize: '11px', color: 'var(--grey-6)', marginBottom: '2px'}}),
      configInfoIcon,
    ], {style: {alignItems: 'center', gap: '4px'}}),
    ui.divV([stratDepthCard.root, stratBreadthCard.root, reagentsModeNote], {style: {gap: '6px'}}),
    ui.form([numRoundsInput, maxComponentsInput, maxRoutesInput]),
    limitsAccordion.root,
    // First pane in the chain — no Back target.
    navRow(null, mkNextBtn(() => accReactionsPane, 'Reactions')),
  ], {style: {gap: '8px'}}), false, null, false);
  accCombinePane.root.classList.add('chem-enum-combine-pane');
  accCombinePane.root.style.setProperty('--chem-enum-label-width', `${sharedLabelWidth}px`);
  accCombinePane.root.style.setProperty('--chem-enum-editor-width', `${sharedEditorWidth}px`);
  // +40 matches the buffer the platform's own d4-dialog-contents branch adds on top of
  // label + input minimums (js-api ui.ts handleFormResize).
  accCombinePane.root.style.setProperty('--chem-enum-form-min-width', `${sharedLabelWidth + formMinInputWidth + 40}px`);
  // Left panel for Preview — content built once previewPanel exists below; this factory itself
  // only runs lazily when the user actually opens the pane, well after that assignment.
  const accPreviewPane = accordion.addPane('Preview', () => ui.divV([
    ui.divText('Samples a small subset of products.',
      {style: {fontSize: '12px', color: 'var(--grey-6)'}}),
    previewPanel.buildRecapCard(),
    // Last pane in the chain — the run action takes the Next slot instead of a target pane.
    navRow(mkBackBtn(() => accExtrasPane, 'Extras'), runControls.previewEnumerateBtnWrap),
  ], {style: {gap: '10px'}}), false, null, false);
  const accPanes = [accReactionsPane, accBbsPane, accExtrasPane, accCombinePane, accPreviewPane];

  // Navigation is chip/button-driven; native click-to-collapse on the header would just empty the
  // panel, so it's disabled at the source.
  for (const p of accPanes) {
    const header = p.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    if (header) {
      header.style.pointerEvents = 'none';
      header.style.cursor = 'default';
    }
  }

  // Subtitle spans injected into each pane's native header — updated by refreshCfgRibbon().
  const injectPaneSub = (pane: DG.AccordionPane): HTMLElement => {
    const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
    const sub = document.createElement('span');
    sub.className = 'chem-enum-pane-subtitle';
    header?.appendChild(sub);
    return sub;
  };
  const subReactions = injectPaneSub(accReactionsPane);
  const subBbs = injectPaneSub(accBbsPane);
  const subExtras = injectPaneSub(accExtrasPane);
  const subCombine = injectPaneSub(accCombinePane);

  // Must fit --chem-enum-form-min-width (the "How to combine" forms' own min-width, computed above)
  // plus the padding/accordion chrome between this pane and those forms — otherwise this pane's
  // content overflows it on every load, not just at a narrow splitter, since the forms now refuse to
  // shrink below their own minimum. The +80 buffer covers this pane's own padding-right plus the
  // nested accordion/pane padding.
  const leftPane = ui.divV([accordion.root],
    {style: {minWidth: `${sharedLabelWidth + formMinInputWidth + 80}px`, overflowY: 'auto', overflowX: 'hidden',
      paddingRight: '8px'}});

  // === Config-summary ribbon (shown above the right-pane tabs) ===
  // Each chip's dot flags "customized": Reactions/Building blocks/Reagents show it when any round
  // has a custom subset (same check behind the round tabs' own dot); Strategy shows it when
  // Combination limits or Product filters differ from platform defaults. Text lives in a child
  // span so refreshing it doesn't wipe the dot out (chip.textContent replaces all children).
  const cfgChipEl = (tooltip: string): {root: HTMLElement; textEl: HTMLElement; dot: HTMLElement} => {
    const dot = mkChangedDot(tooltip);
    const textEl = ui.span([], {style: {overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', minWidth: '0'}});
    const root = ui.divH([dot, textEl], {classes: 'chem-enum-chip', style: {alignItems: 'center', gap: '4px'}});
    return {root, textEl, dot};
  };
  const chipReactionsC = cfgChipEl('One or more rounds use a custom subset of reaction templates.');
  const chipBbsC = cfgChipEl('One or more rounds use a custom subset of building blocks.');
  const chipExtrasC = cfgChipEl('One or more rounds use a custom subset of reagents.');
  const chipCombineC = cfgChipEl('Changed from platform defaults.');
  const cfgEstEl = ui.divText('');
  cfgEstEl.className = 'chem-enum-chip';
  // A separate arrow node per gap — one DOM node can't be reused in three places at once.
  const mkRibbonArrow = (): HTMLElement => {
    const a = ui.iconFA('arrow-right');
    a.classList.add('chem-enum-ribbon-arrow');
    return a;
  };
  // No custom wrapping div and no manual gap/margin — passed straight into their own ribbon group
  // (see setRibbonPanels below); the native ribbon group's own flex layout handles spacing between
  // runControls.runBtn/cancelBtn/progressLabel.

  function refreshCfgRibbon(): void {
    const tDf = templatesInput.value; const bDf = bbsInput.value; const rDf = reagentsInput.value;
    const mode = currentMode();
    const roundsText = roundsLabel(currentRounds());
    const combineText = `${MODE_LABEL[mode]} · ${roundsText}`;
    const setChip = (chip: {root: HTMLElement; textEl: HTMLElement}, text: string, err: boolean): void => {
      chip.textEl.textContent = text;
      chip.root.classList.toggle('chem-enum-chip--err', err);
    };
    setChip(chipReactionsC, tDf ? `${tDf.rowCount} reactions` : 'No reaction table',
      !tDf || !smartsColInput.value);
    setChip(chipBbsC, bDf ? `${bDf.rowCount} BBs` : 'No BBs table',
      !bDf || !bbColInput.value);
    // Extras is fully optional — never flagged as an error state.
    chipExtrasC.textEl.textContent = rDf ? `${rDf.rowCount} reagents` : 'Reagents (None)';
    // "Strategy:" prefix only on the ribbon chip — the accordion pane itself already says "How to combine".
    chipCombineC.textEl.textContent = `Strategy: ${MODE_ABBR[mode]} · ${roundsText}`;
    const n = (tDf && bDf) ? tDf.rowCount * bDf.rowCount : 0;
    cfgEstEl.textContent = n > 0 ? `≈ ${n.toLocaleString('en-US')} products` : '';
    const combChanged = combinationLimitsChanged(config);
    const prodChangedCount = productFiltersChangedCount(config);
    combinationLimitsDot.style.display = combChanged ? '' : 'none';
    productFiltersDot.style.display = prodChangedCount > 0 ? '' : 'none';
    chipCombineC.dot.style.display = (combChanged || prodChangedCount > 0) ? '' : 'none';
    chipReactionsC.dot.style.display = templatesCtl.hasAnyOverride() ? '' : 'none';
    chipBbsC.dot.style.display = bbsCtl.hasAnyOverride() ? '' : 'none';
    chipExtrasC.dot.style.display = reagentsCtl.hasAnyOverride() ? '' : 'none';
    // Re-render Strategy/Preview even when already the visible tab, so in-tab edits stay current.
    // Pass the values just computed above instead of having strategySummary.render() re-derive them.
    if (tabs.currentPane === strategyPane) strategySummary.render(combChanged, prodChangedCount);
    if (tabs.currentPane === previewPane) previewPanel.renderRecap();
    subReactions.textContent = tDf ? `${tDf.rowCount} reactions` : 'No table selected';
    subBbs.textContent = bDf ? `${bDf.rowCount} building blocks` : 'No table selected';
    subExtras.textContent = rDf ? `${rDf.rowCount} reagents` : 'Optional';
    subCombine.textContent = combineText;
  }
  // === Exclusive accordion — only the selected section is shown, on either side ===

  function switchTabForAccPane(pane: DG.AccordionPane): void {
    if (pane === accReactionsPane && templatesPane) {
      tabs.currentPane = templatesPane;
    } else if (pane === accBbsPane && bbsPane) {
      tabs.currentPane = bbsPane;
    } else if (pane === accExtrasPane && reagentsPane) {
      tabs.currentPane = reagentsPane;
    } else if (pane === accCombinePane) {
      tabs.currentPane = strategyPane;
      strategySummary.render();
    } else if (pane === accPreviewPane) {
      tabs.currentPane = previewPane;
      previewPanel.renderRecap();
      previewPanel.refresh();
    }
  }

  // Maps a left-navigator pane to its ribbon chip, so the chip for the currently shown section can
  // be marked active — same `--<state>` modifier convention as chem-enum-chip--err.
  const chipForPane = (pane: DG.AccordionPane): HTMLElement | undefined => {
    if (pane === accReactionsPane) return chipReactionsC.root;
    if (pane === accBbsPane) return chipBbsC.root;
    if (pane === accExtrasPane) return chipExtrasC.root;
    if (pane === accCombinePane) return chipCombineC.root;
    if (pane === accPreviewPane) return cfgEstEl;
    return undefined;
  };

  // Expands `pane` and hides every other section's header entirely (not just collapsed), syncing
  // the right-side tab to match — the ribbon chips are the only navigator now.
  function openAccPaneAndSyncTab(pane: DG.AccordionPane): void {
    accPanes.forEach((p) => {
      p.expanded = (p === pane);
      const header = p.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header) header.style.display = (p === pane) ? '' : 'none';
    });
    const activeChip = chipForPane(pane);
    [chipReactionsC.root, chipBbsC.root, chipExtrasC.root, chipCombineC.root, cfgEstEl].forEach((c) =>
      c.classList.toggle('chem-enum-chip--active', c === activeChip));
    switchTabForAccPane(pane);
  }

  chipReactionsC.root.onclick = () => openAccPaneAndSyncTab(accReactionsPane);
  chipBbsC.root.onclick = () => openAccPaneAndSyncTab(accBbsPane);
  chipExtrasC.root.onclick = () => openAccPaneAndSyncTab(accExtrasPane);
  chipCombineC.root.onclick = () => openAccPaneAndSyncTab(accCombinePane);
  // Initial pane selection must happen near the end of this function, after `tabs`/`strategyPane`
  // exist (switchTabForAccPane reads them) — calling it earlier crashes with a TDZ error.

  // Per-component "Subset by selection" now lives inside each tab's step bar (see DataPanel).

  // ---- Single-grid per-component panel with a per-step strip ----
  // Each data tab shows ONE grid plus a horizontal step strip: "All steps" shows the full library
  // (with the global "Subset by selection"); "Step k" shows a display-only clone whose row selection
  // is that round's subset. Switching chips swaps what the single grid displays — no second grid.
  const dataPanelDeps: DataPanelDeps = {
    view, viewerHost, getConfig: () => config, currentMode, currentRounds,
    refreshValidation, refreshCfgRibbon,
  };
  const templatesCtl = new DataPanel({idx: 0, noun: 'reaction templates',
    input: templatesInput,
    apply: (o, work, cfg) => { o.templates = extractTemplates(cfg, work); }}, dataPanelDeps);
  const bbsCtl = new DataPanel({idx: 1, noun: 'building blocks',
    input: bbsInput,
    apply: (o, work, cfg) => { o.buildingBlocks = extractBuildingBlocks(cfg, work); }}, dataPanelDeps);
  const reagentsCtl = new DataPanel({idx: 2, noun: 'reagents',
    input: reagentsInput, badge: reagentsBadge,
    apply: (o, work, cfg) => { o.reagents = extractReagents(cfg, work); },
    noTableMsg: 'No reagents file selected.', emptyMsg: 'No reagents file selected. Add one in the Extras ' +
      'section to subset reagents per round.'}, dataPanelDeps);
  const templatesPanel = templatesCtl.panel;
  const bbsPanel = bbsCtl.panel;
  const reagentsPanel = reagentsCtl.panel;
  const dataCtls = [templatesCtl, bbsCtl, reagentsCtl];

  // Rebuilds step strips on round-count change or when a component's "All steps" table changes —
  // narrowing "All steps" intentionally discards every step's committed override.
  //
  // `max` on an int input only shows a tooltip, it doesn't clamp — an over-max value is instead
  // caught by validate(); roundCount() separately caps at MAX_ROUNDS so tab count can't blow up.
  view.subs.push(numRoundsInput.onChanged.subscribe(() => refreshValidation()));
  view.subs.push(numRoundsInput.onChanged.pipe(debounceTime(ROUNDS_INPUT_DEBOUNCE_MS)).subscribe(() => {
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

  // buildPerRoundOverrides needs the panels array, which is only known at this call site — the real
  // implementation lives in data-panel.ts (shared with the Strategy summary and Preview recap).
  function buildPerRoundOverrides(cfg: EnumeratorConfig): PerRoundOverride[] | undefined {
    return buildPerRoundOverridesForPanels(dataCtls, cfg);
  }

  // ---- Strategy summary (right-pane view when "How to combine" is the active left section) ----
  const strategySummary = new StrategySummary({
    getConfig: () => config,
    currentMode, currentRounds,
    templatesInput, bbsInput, reagentsInput, exclusionInput,
    buildPerRoundOverrides, overrideCountFor,
  });

  // ---- Preview tab (lazy) ----
  const previewPanel = new PreviewPanel({
    getConfig: () => config,
    currentMode, currentRounds,
    templatesInput, bbsInput, reagentsInput, exclusionInput,
    viewerHost, buildPerRoundOverrides, overrideCountFor, validate,
  });

  // ---- Right pane: TabControl with lazy panes ----
  // addPane's factory is called once when the tab is first activated; we then keep the panel DOM
  // alive across switches.
  const tabs = ui.tabControl(null, false);
  tabs.root.style.width = '100%';
  tabs.root.style.flex = '1 1 0';
  tabs.root.style.minHeight = '0';
  tabs.root.style.overflow = 'hidden';
  templatesPane = tabs.addPane('Reaction templates', () => templatesPanel);
  bbsPane = tabs.addPane('Building blocks', () => bbsPanel);
  reagentsPane = tabs.addPane('Reagents', () => reagentsPanel);
  const strategyPane = tabs.addPane('Strategy', () => strategySummary.panel);
  const previewPane = tabs.addPane('Preview', () => previewPanel.panel);
  view.subs.push(tabs.onTabChanged.subscribe(() => {
    if (tabs.currentPane?.name === 'Preview') previewPanel.refresh();
    // Cancel on tab-away too, or an in-flight preview keeps running unattended (isCancelled only
    // checked at points that read the bumped run id from starting a NEW preview before).
    else previewPanel.cancelPendingRun();
  }));
  // Reagents row-count badge — the only data tab without a count shown elsewhere (reactions/BBs
  // already have it in the ribbon chips and accordion pane subtitles).
  reagentsPane.header.appendChild(reagentsBadge.el);
  // Right-pane tabs mirror the left accordion's sections — hide every tab header so the accordion
  // is the only way to navigate; `tabs.currentPane` is driven entirely by code.
  [templatesPane, bbsPane, reagentsPane, strategyPane, previewPane].forEach((p) => {p.header.style.display = 'none';});
  // "≈ N products" doubles as a shortcut into Preview, through the same navigation path as the chips.
  cfgEstEl.style.cursor = 'pointer';
  cfgEstEl.onclick = () => {if (cfgEstEl.textContent) openAccPaneAndSyncTab(accPreviewPane);};
  ui.tooltip.bind(cfgEstEl, 'Open Preview to sample products before running the full enumeration.');

  // Right pane: a single component-tab control (each tab has its own step strip + one grid).
  // min-width matches the left pane's approach: below it, the grid+filters split squishes into an
  // unreadable sliver rather than shrinking usefully, so it holds a floor and scrolls instead.
  const rightPane = ui.divV([tabs.root], {style: {height: '100%', minWidth: '400px', overflow: 'hidden'}});

  // Resizable horizontal split — drag the divider to rebalance inputs vs side grids.
  const mainRow = ui.splitH([leftPane, rightPane],
    {style: {flex: '1 1 0', minHeight: '0', width: '100%'}}, true);

  // Redundant with the ribbon title and the app-info icon's tooltip — dropped rather than moved.
  const root = ui.divV([
    mainRow,
    validationDiv,
  ], {style: {padding: '0 0 0 16px', height: '100%', boxSizing: 'border-box', overflow: 'hidden'},
    classes: 'chem-enumerator'});

  // Enumerate gets its own ribbon group, never bundled with unrelated items into one custom flex
  // div — a run button sharing a group with unrelated items is what picks up the platform's own
  // group-level background/shadow oddly.
  view.setRibbonPanels([
    [appInfoIcon],
    [runControls.runBtn, runControls.cancelBtn, runControls.progressLabel],
    // Strategy first (the "how"), then reactions/BBs/extras (the "what"), then the resulting estimate.
    [chipCombineC.root, mkRibbonArrow(), chipReactionsC.root, mkRibbonArrow(), chipBbsC.root, mkRibbonArrow(),
      chipExtrasC.root, mkRibbonArrow(), cfgEstEl],
    [loadYamlBtn, saveYamlBtn],
  ]);
  view.append(root);

  // The ribbon's own group/item shadow-background and the chips panel's scroll are both handled
  // declaratively in chem.css via :has(.chem-enum-chip) — no JS reacting to the platform's ribbon
  // re-renders (that pattern leaked onto other views' ribbons and had to keep "winning a race"
  // against re-renders; a plain CSS rule matches by class name regardless of how many times the
  // platform recreates the DOM, and never needs reasserting).
  runControls.setValidation(validate());

  // Looks redundant with buildStepTabs(0)'s own render, but isn't: removing it left the
  // initial grid rendered into a not-yet-sized host and empty on some loads.
  templatesCtl.render();
  bbsCtl.render();
  reagentsCtl.render();
  refreshValidation();
  // Strategy opens first. Must run here, after `tabs`/`strategyPane` exist (see the TDZ note above).
  openAccPaneAndSyncTab(accCombinePane);
  return view;
}
