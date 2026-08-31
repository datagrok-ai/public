import {debounceTime} from 'rxjs/operators';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EnumeratorConfig} from './config';
import {PerRoundOverride} from './enumerate';
import {extractBuildingBlocks, extractReagents, extractTemplates} from './shared';
import {MountedViewerRegistry} from './viewer-mount';
import {StrategySummary} from './strategy-summary';
import {PreviewPanel} from './preview-panel';
import {
  buildPerRoundOverrides as buildPerRoundOverridesForPanels, DataPanel, DataPanelDeps, makeTabBadge,
  overrideCountFor,
} from './data-panel';
import {RunControls} from './run-controls';
import {EnumeratorConfigForm} from './enumerator-config-form';

// Dart int inputs fire onChanged per keystroke — debounce the expensive round-tab rebuild.
const ROUNDS_INPUT_DEBOUNCE_MS = 300;

export async function buildEnumeratorView(): Promise<DG.ViewBase> {
  const view = DG.View.create();
  view.name = 'Reaction Enumerator';
  view.box = true;

  const viewerHost = new MountedViewerRegistry(view);
  const validationDiv = ui.divText('', {style: {color: 'var(--red-3)', fontSize: '12px', flex: '0 0 auto'}});

  // Refresh mediators start as no-ops and are reassigned at the bottom of this function, once
  // everything they coordinate exists. Without the indirection, an input firing onChanged during
  // its own widget construction would reach a not-yet-initialized const.
  const ctx: {
    refreshValidation: () => void; refreshCfgRibbon: () => void; hasAnyPerRoundOverride: () => boolean;
  } = {
    refreshValidation: () => {},
    refreshCfgRibbon: () => {},
    hasAnyPerRoundOverride: () => false,
  };

  // The three functions below are only ever invoked from click handlers, well after construction,
  // so they can reference the tabs and panes declared further down.
  function switchTabForAccPane(pane: DG.AccordionPane): void {
    if (pane === configForm.accReactionsPane)
      tabs.currentPane = templatesPane;
    else if (pane === configForm.accBbsPane)
      tabs.currentPane = bbsPane;
    else if (pane === configForm.accExtrasPane)
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

  const chipForPane = (pane: DG.AccordionPane): HTMLElement | undefined => {
    if (pane === configForm.accReactionsPane) return configForm.chipReactionsC.root;
    if (pane === configForm.accBbsPane) return configForm.chipBbsC.root;
    if (pane === configForm.accExtrasPane) return configForm.chipExtrasC.root;
    if (pane === configForm.accCombinePane) return configForm.chipCombineC.root;
    if (pane === configForm.accPreviewPane) return configForm.cfgEstEl;
    return undefined;
  };

  /** Expands `pane`, hides every other section's header, and syncs the right-side tab to match. */
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

  // Reactions/BBs show their row count in the ribbon chips and pane subtitles already; reagents
  // has neither, so it keeps a tab badge.
  const reagentsBadge = makeTabBadge();

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
  const dataCtls = [templatesCtl, bbsCtl, reagentsCtl];

  view.subs.push(configForm.numRoundsInput.onChanged.subscribe(() => ctx.refreshValidation()));
  view.subs.push(configForm.numRoundsInput.onChanged.pipe(debounceTime(ROUNDS_INPUT_DEBOUNCE_MS)).subscribe(() => {
    dataCtls.forEach((c) => c.onRoundsChanged());
  }));
  view.subs.push(configForm.templatesInput.onChanged.subscribe(() => templatesCtl.onTableChanged()));
  view.subs.push(configForm.bbsInput.onChanged.subscribe(() => bbsCtl.onTableChanged()));
  // A mode switch changes whether a BB override applies at all, so the bars and dots need a refresh
  // even though no table or grid changed.
  view.subs.push(configForm.depthFirstInput.onChanged.subscribe(() => dataCtls.forEach((c) => c.refreshDisplay())));
  view.subs.push(configForm.reagentsInput.onChanged.subscribe(() => {
    reagentsCtl.onTableChanged();
    templatesCtl.refreshDisplay();
    bbsCtl.refreshDisplay();
  }));

  function buildPerRoundOverrides(cfg: EnumeratorConfig): PerRoundOverride[] | undefined {
    return buildPerRoundOverridesForPanels(dataCtls, cfg);
  }

  const strategySummary = new StrategySummary({
    getConfig: configForm.getConfig,
    currentMode: configForm.currentMode, currentRounds: configForm.currentRounds,
    templatesInput: configForm.templatesInput, bbsInput: configForm.bbsInput,
    reagentsInput: configForm.reagentsInput, exclusionInput: configForm.exclusionInput,
    buildPerRoundOverrides, overrideCountFor,
  });

  const previewPanel = new PreviewPanel({
    getConfig: configForm.getConfig,
    currentMode: configForm.currentMode, currentRounds: configForm.currentRounds,
    templatesInput: configForm.templatesInput, bbsInput: configForm.bbsInput,
    reagentsInput: configForm.reagentsInput, exclusionInput: configForm.exclusionInput,
    viewerHost, buildPerRoundOverrides, overrideCountFor, validate: configForm.validate,
  });

  const runControls = new RunControls({
    getConfig: configForm.getConfig,
    templatesInput: configForm.templatesInput, bbsInput: configForm.bbsInput,
    reagentsInput: configForm.reagentsInput, exclusionInput: configForm.exclusionInput,
    validate: configForm.validate, syncQuickInputsToConfig: configForm.syncQuickInputsToConfig,
    buildPerRoundOverrides, refreshValidation: () => ctx.refreshValidation(),
  });

  const tabs = ui.tabControl(null, false);
  tabs.root.style.width = '100%';
  tabs.root.style.flex = '1 1 0';
  tabs.root.style.minHeight = '0';
  tabs.root.style.overflow = 'hidden';
  const templatesPane = tabs.addPane('Reaction templates', () => templatesCtl.panel);
  const bbsPane = tabs.addPane('Building blocks', () => bbsCtl.panel);
  const reagentsPane = tabs.addPane('Reagents', () => reagentsCtl.panel);
  const strategyPane = tabs.addPane('Strategy', () => strategySummary.panel);
  const previewPane = tabs.addPane('Preview', () => previewPanel.panel);
  view.subs.push(tabs.onTabChanged.subscribe(() => {
    if (tabs.currentPane?.name === 'Preview') previewPanel.refresh();
    else previewPanel.cancelPendingRun();
  }));
  reagentsPane.header.appendChild(reagentsBadge.el);
  // The left accordion is the only navigator; tabs.currentPane is driven entirely by code.
  [templatesPane, bbsPane, reagentsPane, strategyPane, previewPane].forEach((p) => {p.header.style.display = 'none';});
  configForm.cfgEstEl.style.cursor = 'pointer';
  configForm.cfgEstEl.onclick = () => {
    if (configForm.cfgEstEl.textContent) openAccPaneAndSyncTab(configForm.accPreviewPane);
  };
  ui.tooltip.bind(configForm.cfgEstEl, 'Open Preview to sample products before running the full enumeration.');

  // Below this min-width the grid+filters split squishes into an unreadable sliver, so it holds a
  // floor and the root scrolls horizontally instead.
  const rightPane = ui.divV([tabs.root], {style: {height: '100%', minWidth: '400px', overflow: 'hidden'}});

  const mainRow = ui.splitH([configForm.leftPane, rightPane],
    {style: {flex: '1 1 0', minHeight: '0', width: '100%'}}, true);

  const root = ui.divV([
    mainRow,
    validationDiv,
  ], {style: {padding: '0 0 0 16px', height: '100%', boxSizing: 'border-box', overflowX: 'auto', overflowY: 'hidden'},
    classes: 'chem-enumerator'});

  const mkRibbonArrow = (): HTMLElement => {
    const a = ui.iconFA('arrow-right');
    a.classList.add('chem-enum-ribbon-arrow');
    return a;
  };
  view.setRibbonPanels([
    [configForm.appInfoIcon],
    [runControls.runBtn, runControls.cancelBtn, runControls.progressLabel],
    [configForm.chipCombineC.root, mkRibbonArrow(), configForm.chipReactionsC.root, mkRibbonArrow(),
      configForm.chipBbsC.root, mkRibbonArrow(), configForm.chipExtrasC.root, mkRibbonArrow(), configForm.cfgEstEl],
    [configForm.loadYamlBtn, configForm.saveYamlBtn],
  ]);
  view.append(root);

  ctx.hasAnyPerRoundOverride = () => dataCtls.some((p) => p.hasAnyOverride());
  ctx.refreshValidation = (): void => {
    // validate() syncs the quick inputs into config first, so the refreshes below read current values.
    const err = configForm.validate();
    ctx.refreshCfgRibbon();
    configForm.refreshStrategyCards();
    validationDiv.textContent = err ?? '';
    runControls.setValidation(err);
  };
  ctx.refreshCfgRibbon = (): void => {
    configForm.refreshRibbonChips({
      templatesOverride: templatesCtl.hasAnyOverride(),
      bbsOverride: bbsCtl.hasAnyOverride(),
      reagentsOverride: reagentsCtl.hasAnyOverride(),
    });
    // Re-render even when already visible, so in-tab edits stay current.
    if (tabs.currentPane === strategyPane) strategySummary.render();
    if (tabs.currentPane === previewPane) previewPanel.renderRecap();
  };

  runControls.setValidation(configForm.validate());

  // Must run after view.append(root): a grid rendered into an unsized host comes out empty.
  templatesCtl.render();
  bbsCtl.render();
  reagentsCtl.render();
  ctx.refreshValidation();
  openAccPaneAndSyncTab(configForm.accCombinePane);
  return view;
}
