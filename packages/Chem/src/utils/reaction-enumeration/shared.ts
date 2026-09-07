import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {OutputRow, OutputStep, TemplateInput} from './enumerate';

export type Mode = 'depth' | 'breadth' | 'reagents';
export type DataKey = 'templates' | 'buildingBlocks' | 'reagents';

export const MODE_LABEL = {depth: 'Depth-first', breadth: 'Breadth-first', reagents: 'Reagents'} as const;
export const roundsLabel = (n: number): string => `${n} step${n === 1 ? '' : 's'}`;

export const OVERRIDE_DOT_COLOR = 'var(--orange-2, #c98a1b)';
export const CHANGED_DOT_STYLE = {width: '6px', height: '6px', borderRadius: '50%', background: OVERRIDE_DOT_COLOR};

export const MAX_ROUNDS = 10;

/** Every per-round loop that builds DOM must clamp: "Number of steps" can transiently hold a much
 * larger typed value, which would freeze the tab building thousands of rows per keystroke. */
export function clampRounds(rounds: number): number {
  return Math.min(MAX_ROUNDS, Math.max(1, rounds));
}

export function combinationLimitsChanged(cfg: EnumeratorConfig): boolean {
  return cfg.max_num_components !== DEFAULT_CONFIG.max_num_components ||
    cfg.max_num_routes_per_compound !== DEFAULT_CONFIG.max_num_routes_per_compound ||
    cfg.max_num_combinations_per_template !== DEFAULT_CONFIG.max_num_combinations_per_template ||
    cfg.keep_building_blocks_in_final_output !== DEFAULT_CONFIG.keep_building_blocks_in_final_output;
}

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

/** Naive upper bound shown in the ribbon and the Strategy summary. */
export function estimateProductCount(tDf: DG.DataFrame | null, bDf: DG.DataFrame | null): number {
  return (tDf && bDf) ? tDf.rowCount * bDf.rowCount : 0;
}

/** Hint/status header bar shared by every right-pane tab. */
export const panelHeader = (hint: string, status?: HTMLElement): HTMLElement => {
  const hintEl = ui.divText(hint, {style: {
    fontSize: '11px', color: 'var(--grey-5)', flex: '0 0 auto', marginRight: '4px',
  }});
  const children: HTMLElement[] = [hintEl];
  if (status) children.push(status);
  return ui.div(children, {style: {
    display: 'flex', alignItems: 'center', gap: '8px', flex: '0 0 auto',
    padding: '4px 8px 5px', borderBottom: '1px solid var(--grey-2)',
  }});
};

/** `scrollable` is for plain content hosts; grids scroll internally instead and get a bottom fade.
 * min-height:0 lets this flex child shrink below its content height. */
export const tabPanel = (header: HTMLElement, gridHost: HTMLElement, scrollable = false): HTMLElement => {
  gridHost.style.display = 'flex';
  gridHost.style.flexDirection = 'column';
  gridHost.style.flex = '1 1 0';
  gridHost.style.minHeight = '0';
  if (scrollable) {
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

/** Tags reaction columns by sampling for `>>`, then runs one auto-detection pass over the frame.
 * Returns that pass's promise so callers can await it instead of scanning twice. */
export function detectChemSemTypes(df: DG.DataFrame): Promise<void> {
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
  return df.meta.detectSemanticTypes();
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

const STEP_RXN_PREFIX = '~reaction_name_';

export function buildResultDataFrame(rows: OutputRow[], name = 'Enumeration result'): DG.DataFrame {
  const maxSteps = rows.reduce((m, r) => Math.max(m, r.steps.length), 0);
  // maxSteps is the longest route; shorter routes leave their tail cells blank.
  const perStep = (colName: (k: number) => string, pick: (s: OutputStep) => string): DG.Column[] =>
    Array.from({length: maxSteps}, (_, k) =>
      DG.Column.fromStrings(colName(k + 1), rows.map((r) => (r.steps[k] ? pick(r.steps[k]) : ''))));
  // Per-step templates stay visible (a SMARTS is only readable drawn); per-step reaction names are
  // hidden and exist to feed the step hierarchy filter.
  const templateCols = perStep((k) => `template_${k}`, (s) => s.template);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('product', rows.map((r) => r.product)),
    DG.Column.fromStrings('route', rows.map((r) => r.route)),
    DG.Column.fromStrings('reaction_names',
      rows.map((r) => r.steps.map((s, i) => `Step ${i + 1}: ${s.reactionName}`).join('\n'))),
    ...templateCols,
    DG.Column.fromInt32Array('round', new Int32Array(rows.map((r) => r.round))),
    DG.Column.fromInt32Array('n_routes', new Int32Array(rows.map((r) => r.n_routes))),
    ...perStep((k) => `${STEP_RXN_PREFIX}${k}`, (s) => s.reactionName),
  ]);
  df.name = name;
  df.col('product')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('route')!.semType = 'ChemicalReaction';
  for (const c of templateCols) c.semType = 'ChemicalReaction';
  // A filter on this column then offers each "Step k: …" line as its own category.
  df.col('reaction_names')!.meta.multiValueSeparator = '\n';
  return df;
}

/** The filter type the platform picks for a column, except that a zero-variance numeric column gets
 * the categorical filter: the histogram one clobbers the whole frame's .filter to all-false. */
export function defaultFilterState(col: DG.Column): DG.FilterState {
  const type = col.isNumerical ?
    (col.stats.min === col.stats.max ? DG.FILTER_TYPE.CATEGORICAL : DG.FILTER_TYPE.HISTOGRAM) :
    col.semType === DG.SEMTYPE.MOLECULE ? DG.FILTER_TYPE.SUBSTRUCTURE : DG.FILTER_TYPE.CATEGORICAL;
  return {type, column: col.name};
}

/** Per-column filters for the result view, minus `route` (every value is distinct) and, when there is
 * a step tree, minus the summary column it duplicates; plus a "step 1 reaction → step 2 reaction → …"
 * tree over the hidden per-step columns, which the hierarchical filter's own column picker (visible
 * columns only) could not build. */
export function addResultFilters(tv: DG.TableView): void {
  const df = tv.dataFrame;
  const stepCols = df.columns.names().filter((c) => c.startsWith(STEP_RXN_PREFIX));
  const tree = stepCols.length > 1;
  const skip = new Set(tree ? ['route', 'reaction_names'] : ['route']);
  const fg = tv.getFiltersGroup({createDefaultFilters: false});
  // The newest filter shows on top, so add bottom-up: the panel then follows column order, tree first.
  const cols = df.columns.toList().filter((c) => !c.name.startsWith('~') && !skip.has(c.name)).reverse();
  for (const c of cols) fg.updateOrAdd(defaultFilterState(c), false);
  // 'hierarchical' isn't in DG.FILTER_TYPE; the state shape is what the group itself serializes.
  if (tree) fg.updateOrAdd({type: 'hierarchical', colNames: stepCols}, false);
}

export interface BuiltInputs {
  templates: TemplateInput[];
  buildingBlocks: string[];
  exclusionSmarts: string[];
  reagents: string[];
}

/** Also runs against a per-round subset table, hence the separate export. */
export function extractTemplates(config: EnumeratorConfig, tDf: DG.DataFrame): TemplateInput[] {
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

export function extractBuildingBlocks(config: EnumeratorConfig, bDf: DG.DataFrame): string[] {
  return getStringColumn(bDf, config.enumeration.bb_smiles_column).filter((s) => s.trim().length > 0);
}

export function extractReagents(config: EnumeratorConfig, rDf: DG.DataFrame): string[] {
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
