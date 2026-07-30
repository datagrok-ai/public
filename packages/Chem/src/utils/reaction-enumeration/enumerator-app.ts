/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import {debounceTime} from 'rxjs/operators';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DEFAULT_CONFIG, EnumeratorConfig} from './config';
import {OutputRow, PerRoundOverride, TemplateInput} from './enumerate';
import {MountedViewerRegistry} from './viewer-mount';
import {StrategySummary} from './strategy-summary';
import {PreviewPanel} from './preview-panel';
import {
  buildPerRoundOverrides as buildPerRoundOverridesForPanels, DataPanel, DataPanelDeps, makeTabBadge,
  overrideCountFor,
} from './data-panel';
import {RunControls} from './run-controls';
import {EnumeratorConfigForm} from './enumerator-config-form';

// Shared "custom subset" indicator color — round tabs' dot and the Strategy summary's dot.
export const OVERRIDE_DOT_COLOR = 'var(--orange-2, #c98a1b)';
// Shared look for the small "changed/custom" dots; call sites add their own display mode and spacing.
export const CHANGED_DOT_STYLE = {width: '6px', height: '6px', borderRadius: '50%', background: OVERRIDE_DOT_COLOR};

// Single source of truth for the enumeration-mode/data-key literal unions — every module that
// needs them imports from here instead of redeclaring the same literals.
export type Mode = 'depth' | 'breadth' | 'reagents';
export type DataKey = 'templates' | 'buildingBlocks' | 'reagents';

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

// Shared naive product-count estimate — used by the ribbon chip and the Strategy summary so they
// can't drift out of sync (e.g. one estimate rounding differently from the other).
export function estimateProductCount(tDf: DG.DataFrame | null, bDf: DG.DataFrame | null): number {
  return (tDf && bDf) ? tDf.rowCount * bDf.rowCount : 0;
}

// Shared panel chrome — a hint/status header bar plus a content host — used by every right-pane
// tab (data grids, Strategy summary, Preview).
export const panelHeader = (hint: string, status?: HTMLElement): HTMLElement => {
  // flex:0 0 auto (not 1 1 auto) so hint and status sit side by side, both left-aligned, instead
  // of hint growing to push status to the far right of the row.
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
export function detectChemSemTypes(df: DG.DataFrame): Promise<void> {
  // detectSemanticTypes() scans the WHOLE dataframe; calling it per-column made this O(columns²)
  // and it ran on every step-clone. Tag ChemicalReaction columns first, then auto-detect once.
  // Returns the scan's own promise so a caller that needs to wait for it (e.g. loadBundledCsv) can
  // await this call directly instead of triggering a second, redundant full-table scan.
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

  const viewerHost = new MountedViewerRegistry(view);
  const validationDiv = ui.divText('', {style: {color: 'var(--red-3)', fontSize: '12px', flex: '0 0 auto'}});

  // ---- Late-bound refresh mediators ----
  // refreshValidation/refreshCfgRibbon are threaded into every class below as a constructor dep, but
  // their real implementations need templatesCtl/bbsCtl/reagentsCtl/strategySummary/previewPanel/
  // runControls/tabs/strategyPane/previewPane — all constructed AFTER those deps are handed out.
  // Every class calls `ctx.refreshValidation()`/`ctx.refreshCfgRibbon()` (never captures the
  // function itself at construction time, or it would freeze on today's no-op) — so a premature
  // synchronous call (e.g. a DG input auto-selecting a value during its own construction, firing
  // onChanged inline) is now a harmless no-op instead of a "Cannot access 'X' before
  // initialization" crash. Reassigned to the real implementations once everything they coordinate
  // exists — see the bottom of this function.
  const ctx: {
    refreshValidation: () => void; refreshCfgRibbon: () => void; hasAnyPerRoundOverride: () => boolean;
  } = {
    refreshValidation: () => {},
    refreshCfgRibbon: () => {},
    hasAnyPerRoundOverride: () => false,
  };

  // switchTabForAccPane/chipForPane/openAccPaneAndSyncTab stay plain hoisted function statements
  // (not part of `ctx`): they're only ever invoked from explicit click handlers, never from a
  // widget's own synchronous construction-time event, so they don't share refreshValidation's risk.
  // Right-pane tab references — assigned when tabs are built; used by section-open handlers for
  // context-sensitive tab switching. Declared here so openAccPaneAndSyncTab can close over them.
  let templatesPane: DG.TabPane | undefined;
  let bbsPane: DG.TabPane | undefined;
  let reagentsPane: DG.TabPane | undefined;

  function switchTabForAccPane(pane: DG.AccordionPane): void {
    if (pane === configForm.accReactionsPane && templatesPane)
      tabs.currentPane = templatesPane;
    else if (pane === configForm.accBbsPane && bbsPane)
      tabs.currentPane = bbsPane;
    else if (pane === configForm.accExtrasPane && reagentsPane)
      tabs.currentPane = reagentsPane;
    else if (pane === configForm.accCombinePane) {
      tabs.currentPane = strategyPane;
      strategySummary.render();
    } else if (pane === configForm.accPreviewPane) {
      tabs.currentPane = previewPane;
      previewPanel.renderRecap();
      previewPanel.refresh();
    }
  }

  // Maps a left-navigator pane to its ribbon chip, so the chip for the currently shown section can
  // be marked active — same `--<state>` modifier convention as chem-enum-chip--err.
  const chipForPane = (pane: DG.AccordionPane): HTMLElement | undefined => {
    if (pane === configForm.accReactionsPane) return configForm.chipReactionsC.root;
    if (pane === configForm.accBbsPane) return configForm.chipBbsC.root;
    if (pane === configForm.accExtrasPane) return configForm.chipExtrasC.root;
    if (pane === configForm.accCombinePane) return configForm.chipCombineC.root;
    if (pane === configForm.accPreviewPane) return configForm.cfgEstEl;
    return undefined;
  };

  // Expands `pane` and hides every other section's header entirely (not just collapsed), syncing
  // the right-side tab to match — the ribbon chips are the only navigator now.
  function openAccPaneAndSyncTab(pane: DG.AccordionPane): void {
    configForm.accPanes.forEach((p) => {
      p.expanded = (p === pane);
      const header = p.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header) header.style.display = (p === pane) ? '' : 'none';
    });
    const activeChip = chipForPane(pane);
    [configForm.chipReactionsC.root, configForm.chipBbsC.root, configForm.chipExtrasC.root,
      configForm.chipCombineC.root, configForm.cfgEstEl].forEach((c) =>
      c.classList.toggle('chem-enum-chip--active', c === activeChip));
    switchTabForAccPane(pane);
  }

  // ---- Config form: config itself, data/quick-config inputs, YAML load/save, validation, the
  // strategy cards, and the left-nav accordion (5 panes) ----
  const configForm = await EnumeratorConfigForm.create({
    view, viewerHost, refreshValidation: () => ctx.refreshValidation(), openAccPaneAndSyncTab,
    getPreviewRecapCard: () => previewPanel.buildRecapCard(),
    getPreviewEnumerateBtnWrap: () => runControls.previewEnumerateBtnWrap,
    hasAnyPerRoundOverride: () => ctx.hasAnyPerRoundOverride(),
  });

  configForm.chipReactionsC.root.onclick = () => openAccPaneAndSyncTab(configForm.accReactionsPane);
  configForm.chipBbsC.root.onclick = () => openAccPaneAndSyncTab(configForm.accBbsPane);
  configForm.chipExtrasC.root.onclick = () => openAccPaneAndSyncTab(configForm.accExtrasPane);
  configForm.chipCombineC.root.onclick = () => openAccPaneAndSyncTab(configForm.accCombinePane);
  // Initial pane selection must happen near the end of this function, after `tabs`/`strategyPane`
  // exist (switchTabForAccPane reads them) — calling it earlier crashes with a TDZ error.

  // Tab row-count badge. Reactions/BBs already show their row count via the always-visible ribbon
  // chips (chipReactionsC/chipBbsC) and the accordion pane subtitles — a tab badge there would just
  // repeat the same number a third time. Reagents has neither, so it keeps a badge.
  const reagentsBadge = makeTabBadge();

  // Per-component "Subset by selection" now lives inside each tab's step bar (see DataPanel).

  // ---- Single-grid per-component panel with a per-step strip (see DataPanel's class doc) ----
  const dataPanelDeps: DataPanelDeps = {
    view, viewerHost, getConfig: configForm.getConfig, currentMode: configForm.currentMode,
    currentRounds: configForm.currentRounds,
    refreshValidation: () => ctx.refreshValidation(), refreshCfgRibbon: () => ctx.refreshCfgRibbon(),
  };
  const templatesCtl = new DataPanel({idx: 0, noun: 'reaction templates',
    input: configForm.templatesInput,
    apply: (o, work, cfg) => {o.templates = extractTemplates(cfg, work);}}, dataPanelDeps);
  const bbsCtl = new DataPanel({idx: 1, noun: 'building blocks',
    input: configForm.bbsInput,
    apply: (o, work, cfg) => {o.buildingBlocks = extractBuildingBlocks(cfg, work);}}, dataPanelDeps);
  const reagentsCtl = new DataPanel({idx: 2, noun: 'reagents',
    input: configForm.reagentsInput, badge: reagentsBadge,
    apply: (o, work, cfg) => {o.reagents = extractReagents(cfg, work);},
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
  view.subs.push(configForm.numRoundsInput.onChanged.subscribe(() => ctx.refreshValidation()));
  view.subs.push(configForm.numRoundsInput.onChanged.pipe(debounceTime(ROUNDS_INPUT_DEBOUNCE_MS)).subscribe(() => {
    dataCtls.forEach((c) => c.onRoundsChanged());
  }));
  view.subs.push(configForm.templatesInput.onChanged.subscribe(() => templatesCtl.onTableChanged()));
  view.subs.push(configForm.bbsInput.onChanged.subscribe(() => bbsCtl.onTableChanged()));
  // BB override dot/status is mode-aware (hasOverride hides it in breadth-first) — refresh the bar +
  // dots on any mode switch so they don't show stale state. Neither the underlying table nor the
  // grid changes here, so refreshDisplay (not render) is enough — no grid rebuild needed.
  view.subs.push(configForm.depthFirstInput.onChanged.subscribe(() => dataCtls.forEach((c) => c.refreshDisplay())));
  view.subs.push(configForm.reagentsInput.onChanged.subscribe(() => {
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
    getConfig: configForm.getConfig,
    currentMode: configForm.currentMode, currentRounds: configForm.currentRounds,
    templatesInput: configForm.templatesInput, bbsInput: configForm.bbsInput,
    reagentsInput: configForm.reagentsInput, exclusionInput: configForm.exclusionInput,
    buildPerRoundOverrides, overrideCountFor,
  });

  // ---- Preview tab (lazy) ----
  const previewPanel = new PreviewPanel({
    getConfig: configForm.getConfig,
    currentMode: configForm.currentMode, currentRounds: configForm.currentRounds,
    templatesInput: configForm.templatesInput, bbsInput: configForm.bbsInput,
    reagentsInput: configForm.reagentsInput, exclusionInput: configForm.exclusionInput,
    viewerHost, buildPerRoundOverrides, overrideCountFor, validate: configForm.validate,
  });

  // ---- Run / Cancel ----
  const runControls = new RunControls({
    getConfig: configForm.getConfig,
    templatesInput: configForm.templatesInput, bbsInput: configForm.bbsInput,
    reagentsInput: configForm.reagentsInput, exclusionInput: configForm.exclusionInput,
    validate: configForm.validate, syncQuickInputsToConfig: configForm.syncQuickInputsToConfig,
    buildPerRoundOverrides, refreshValidation: () => ctx.refreshValidation(),
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
  configForm.cfgEstEl.style.cursor = 'pointer';
  configForm.cfgEstEl.onclick = () => {
    if (configForm.cfgEstEl.textContent) openAccPaneAndSyncTab(configForm.accPreviewPane);
  };
  ui.tooltip.bind(configForm.cfgEstEl, 'Open Preview to sample products before running the full enumeration.');

  // Right pane: a single component-tab control (each tab has its own step strip + one grid).
  // min-width matches the left pane's approach: below it, the grid+filters split squishes into an
  // unreadable sliver rather than shrinking usefully, so it holds a floor and scrolls instead.
  const rightPane = ui.divV([tabs.root], {style: {height: '100%', minWidth: '400px', overflow: 'hidden'}});

  // Resizable horizontal split — drag the divider to rebalance inputs vs side grids.
  const mainRow = ui.splitH([configForm.leftPane, rightPane],
    {style: {flex: '1 1 0', minHeight: '0', width: '100%'}}, true);

  // Redundant with the ribbon title and the app-info icon's tooltip — dropped rather than moved.
  const root = ui.divV([
    mainRow,
    validationDiv,
  ], {style: {padding: '0 0 0 16px', height: '100%', boxSizing: 'border-box', overflow: 'hidden'},
    classes: 'chem-enumerator'});

  // A separate arrow node per gap — one DOM node can't be reused in three places at once.
  const mkRibbonArrow = (): HTMLElement => {
    const a = ui.iconFA('arrow-right');
    a.classList.add('chem-enum-ribbon-arrow');
    return a;
  };
  // Enumerate gets its own ribbon group, never bundled with unrelated items into one custom flex
  // div — a run button sharing a group with unrelated items is what picks up the platform's own
  // group-level background/shadow oddly.
  view.setRibbonPanels([
    [configForm.appInfoIcon],
    [runControls.runBtn, runControls.cancelBtn, runControls.progressLabel],
    // Strategy first (the "how"), then reactions/BBs/extras (the "what"), then the resulting estimate.
    [configForm.chipCombineC.root, mkRibbonArrow(), configForm.chipReactionsC.root, mkRibbonArrow(),
      configForm.chipBbsC.root, mkRibbonArrow(), configForm.chipExtrasC.root, mkRibbonArrow(), configForm.cfgEstEl],
    [configForm.loadYamlBtn, configForm.saveYamlBtn],
  ]);
  view.append(root);

  // ---- Bind the real mediator implementations ----
  // Every class above received a `ctx.refreshValidation()`/`ctx.refreshCfgRibbon()` indirection
  // rather than a direct function reference — reassigning here (now that everything they coordinate
  // exists) takes effect immediately for all of them, with no need to re-thread anything.
  ctx.hasAnyPerRoundOverride = () => dataCtls.some((p) => p.hasAnyOverride());
  ctx.refreshValidation = (): void => {
    // validate() syncs quick inputs into config as its own first step — call it before
    // refreshCfgRibbon()/refreshStrategyCards() so they read the just-updated config, instead of
    // syncing here too and reading one refresh behind.
    const err = configForm.validate();
    ctx.refreshCfgRibbon();
    configForm.refreshStrategyCards();
    validationDiv.textContent = err ?? '';
    runControls.setValidation(err);
  };
  ctx.refreshCfgRibbon = (): void => {
    // Re-render Strategy/Preview even when already the visible tab, so in-tab edits stay current.
    configForm.refreshRibbonChips({
      templatesOverride: templatesCtl.hasAnyOverride(),
      bbsOverride: bbsCtl.hasAnyOverride(),
      reagentsOverride: reagentsCtl.hasAnyOverride(),
    });
    if (tabs.currentPane === strategyPane) strategySummary.render();
    if (tabs.currentPane === previewPane) previewPanel.renderRecap();
  };

  // The ribbon's own group/item shadow-background and the chips panel's scroll are both handled
  // declaratively in chem.css via :has(.chem-enum-chip) — no JS reacting to the platform's ribbon
  // re-renders (that pattern leaked onto other views' ribbons and had to keep "winning a race"
  // against re-renders; a plain CSS rule matches by class name regardless of how many times the
  // platform recreates the DOM, and never needs reasserting).
  runControls.setValidation(configForm.validate());

  // Looks redundant with buildStepTabs(0)'s own render, but isn't: removing it left the
  // initial grid rendered into a not-yet-sized host and empty on some loads.
  templatesCtl.render();
  bbsCtl.render();
  reagentsCtl.render();
  ctx.refreshValidation();
  // Strategy opens first. Must run here, after `tabs`/`strategyPane` exist (see the TDZ note above).
  openAccPaneAndSyncTab(configForm.accCombinePane);
  return view;
}
