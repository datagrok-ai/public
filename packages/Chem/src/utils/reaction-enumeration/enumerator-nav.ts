/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {CHANGED_DOT_STYLE, Mode} from './enumerator-app';

export type ChipEl = {root: HTMLElement; textEl: HTMLElement; dot: HTMLElement};
type StratCard = {root: HTMLElement; icon: HTMLElement};

// Plain-data snapshot of everything the ribbon chips + accordion-pane subtitles need to show —
// computed by EnumeratorConfigForm (which owns EnumeratorConfig) and handed in wholesale, so this
// class never needs to read EnumeratorConfig itself.
export interface RibbonChipState {
  reactionsText: string; reactionsOverride: boolean; reactionsSubtitle: string;
  bbsText: string; bbsOverride: boolean; bbsSubtitle: string;
  extrasText: string; extrasOverride: boolean; extrasSubtitle: string;
  combineChipText: string; combineSubtitle: string;
  estimateText: string;
  limitsChanged: boolean; filtersChanged: boolean;
}

export interface EnumeratorNavDeps {
  view: DG.View;
  templatesInput: DG.InputBase<DG.DataFrame | null>;
  smartsColInput: DG.InputBase<DG.Column | null>;
  blockingColInput: DG.InputBase<DG.Column | null>;
  rxnNameColInput: DG.InputBase<DG.Column | null>;
  bbsInput: DG.InputBase<DG.DataFrame | null>;
  bbColInput: DG.InputBase<DG.Column | null>;
  reagentsInput: DG.InputBase<DG.DataFrame | null>;
  reagentsColInput: DG.InputBase<DG.Column | null>;
  exclusionInput: DG.InputBase<DG.DataFrame | null>;
  exclusionColInput: DG.InputBase<DG.Column | null>;
  numRoundsInput: DG.InputBase<number | null>;
  maxComponentsInput: DG.InputBase<number | null>;
  maxRoutesInput: DG.InputBase<number | null>;
  depthFirstInput: DG.InputBase<boolean>;
  configInfoIcon: HTMLElement;
  // Getters, not snapshots — combinationLimitFields/productFilterFields are rebuilt wholesale on
  // every YAML load (see EnumeratorConfigForm.syncConfigToQuickInputs), so the lazy form builder
  // below must re-read them at build/invalidate time, not capture the array that existed at
  // construction.
  getCombinationLimitInputs: () => DG.InputBase<unknown>[];
  getProductFilterInputs: () => DG.InputBase<unknown>[];
  setAndFire: <T>(input: DG.InputBase<T>, v: T) => void;
  // Orchestrator-level nav mediator — see EnumeratorConfigFormDeps for the same TDZ-safety note;
  // referenced (not called) by nav buttons at construction time, only invoked later on click.
  openAccPaneAndSyncTab: (pane: DG.AccordionPane) => void;
  getPreviewRecapCard: () => HTMLElement;
  getPreviewEnumerateBtnWrap: () => HTMLElement;
}

/** Owns the strategy picker cards and the left-nav accordion (5 panes forming one navigational
 * ring), plus the ribbon chip DOM and accordion-pane subtitles — bundled together because both are
 * just two surfaces for the same "which section is active / what does it currently hold" state,
 * kept in sync by one applyRibbonState() call. Pure layout: never reads EnumeratorConfig itself,
 * only the already-built inputs and plain values EnumeratorConfigForm computes and hands in. */
export class EnumeratorNav {
  readonly accReactionsPane: DG.AccordionPane;
  readonly accBbsPane: DG.AccordionPane;
  readonly accExtrasPane: DG.AccordionPane;
  readonly accCombinePane: DG.AccordionPane;
  readonly accPreviewPane: DG.AccordionPane;
  readonly accPanes: DG.AccordionPane[];
  readonly leftPane: HTMLElement;

  readonly chipReactionsC: ChipEl;
  readonly chipBbsC: ChipEl;
  readonly chipExtrasC: ChipEl;
  readonly chipCombineC: ChipEl;
  readonly cfgEstEl: HTMLElement;

  private readonly stratDepthCard: StratCard;
  private readonly stratBreadthCard: StratCard;
  private readonly reagentsModeNote: HTMLElement;

  private readonly combinationLimitsForm: {getRoot: () => HTMLElement; invalidate: () => void};
  private readonly productFilterForm: {getRoot: () => HTMLElement; invalidate: () => void};
  private readonly combinationLimitsDot: HTMLElement;
  private readonly productFiltersDot: HTMLElement;

  private readonly subReactions: HTMLElement;
  private readonly subBbs: HTMLElement;
  private readonly subExtras: HTMLElement;
  private readonly subCombine: HTMLElement;

  constructor(private readonly deps: EnumeratorNavDeps) {
    // Strategy cards replace the depth-first checkbox: depth/breadth drive the hidden `depthFirstInput`
    // (existing sync/validation keeps working). Depth/breadth are disabled (reagentsModeNote shown
    // instead) whenever a reagents file is loaded, since reagents mode makes strategy irrelevant.
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
    this.stratDepthCard = buildStratCard('arrow-right', 'Depth-first',
      'Grow each product with original blocks — linear chains.');
    this.stratBreadthCard = buildStratCard('sitemap', 'Breadth-first',
      'Combine any earlier products — convergent routes.');
    // Reagents mode is auto-derived (set via a file in Extras), never manually selectable here.
    this.reagentsModeNote = ui.divH([ui.iconFA('info-circle'),
      ui.span([' Reagents mode active — set via a reagents file in Extras.'])],
    {style: {fontSize: '11px', color: 'var(--grey-6)', gap: '6px', alignItems: 'center',
      padding: '6px 10px', display: 'none'}});
    ui.tooltip.bind(this.reagentsModeNote, 'Every round uses exactly one building block or earlier-round product, ' +
      'with reagents filling the remaining slots. Automatically selected as soon as a reagents file is loaded ' +
      'in the Extras tab — to go back to Depth-first/Breadth-first, clear the reagents file there.');

    this.stratDepthCard.root.onclick = (): void => {
      if (this.deps.reagentsInput.value != null) return;
      this.deps.setAndFire(this.deps.depthFirstInput, true);
    };
    this.stratBreadthCard.root.onclick = (): void => {
      if (this.deps.reagentsInput.value != null) return;
      this.deps.setAndFire(this.deps.depthFirstInput, false);
    };

    const accordion = ui.accordion();
    accordion.root.classList.add('chem-enum-accordion');

    // Target pane resolved lazily via a thunk: Reactions pane is added expanded, so its factory runs
    // synchronously inside addPane, before later panes exist — capturing directly would hit the TDZ.
    const mkNavBtn = (kind: 'next' | 'back', getTarget: () => DG.AccordionPane, label: string): HTMLElement => {
      const btn = ui.button(`${kind === 'next' ? 'Next' : 'Back'}: ${label}`,
        () => this.deps.openAccPaneAndSyncTab(getTarget()));
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

    // allowDragOut (5th arg) defaults to true; panes shouldn't be draggable out of this panel.
    // One shared form so all four fields' labels align (two forms would size independently).
    this.accReactionsPane = accordion.addPane('Reactions', () =>
      ui.divV([ui.form([this.deps.templatesInput, this.deps.smartsColInput, this.deps.blockingColInput, this.deps.rxnNameColInput]),
        navRow(mkBackBtn(() => this.accCombinePane, 'How to combine'), mkNextBtn(() => this.accBbsPane, 'Building blocks'))]),
    true, null, false);
    this.accBbsPane = accordion.addPane('Building blocks', () =>
      ui.divV([ui.form([this.deps.bbsInput, this.deps.bbColInput]),
        navRow(mkBackBtn(() => this.accReactionsPane, 'Reactions'), mkNextBtn(() => this.accExtrasPane, 'Extras'))]),
    false, null, false);
    const extrasForm = ui.form([this.deps.reagentsInput, this.deps.reagentsColInput, this.deps.exclusionInput, this.deps.exclusionColInput]);
    this.accExtrasPane = accordion.addPane('Extras', () =>
      ui.divV([extrasForm,
        navRow(mkBackBtn(() => this.accBbsPane, 'Building blocks'), mkNextBtn(() => this.accPreviewPane, 'Preview'))]),
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
    // Builds ui.form() lazily, only once the pane's content factory first runs — building it while
    // still detached would measure it at 0 width and get it marked .ui-form-condensed regardless of
    // label width. invalidate() rebuilds in place after a config reload swaps in fresh inputs.
    const lazyFilterForm = (getInputs: () => DG.InputBase<unknown>[]):
    {getRoot: () => HTMLElement; invalidate: () => void} => {
      let root: HTMLElement | null = null;
      const build = (): HTMLElement => ui.form(getInputs());
      return {
        getRoot: (): HTMLElement => root ??= build(),
        invalidate: (): void => {
          if (!root) return;
          const fresh = build();
          root.replaceWith(fresh);
          root = fresh;
        },
      };
    };
    // Shared label-column width across all three "How to combine" forms so their columns align.
    // Set via a CSS custom property (not inline style) — the platform's own async label
    // auto-sizing would otherwise overwrite a plain inline style after mount.
    const CHEM_ENUM_LABEL_WIDTH = 220;
    // Editors need pinning too, not just labels: the platform widens an editor from its normal ~150px
    // to ~176px on whichever of the three forms it currently considers .ui-form-condensed, so without
    // this the Product filters column can end up visibly wider than Combination limits' even with
    // labels aligned. 150px matches the platform's own un-condensed default.
    const sharedEditorWidth = 150;
    // Fixed floor for the whole form's width, comfortably above what the platform's own
    // .ui-form-condensed decision (js-api ui.ts's getInputsMinWidths) would compute for the widest
    // input here (a string input like "Only these atoms allowed", whose own minimum is ~200px) plus
    // the label column above — reserving this much up front means these forms never cross that
    // threshold and never condense, regardless of where the left/right splitter sits when a section
    // is first expanded.
    const CHEM_ENUM_FORM_MIN_WIDTH = 480;
    // Combination limits/Product filters are panes of limitsAccordion, one accordion level below this
    // one — the platform's own .d4-accordion (6px padding) + .d4-accordion-pane-content (20px margin)
    // indent that extra level, so without a matching offset here the rounds/components form's label
    // column starts 26px to the left of the otherwise-identical 220px column below it.
    const CHEM_ENUM_NESTED_ACCORDION_INDENT = 26;
    // Independently-collapsible sub-sections within "How to combine" (no forced exclusivity, unlike
    // the outer wizard-navigation accordion).
    const limitsAccordion = ui.accordion();
    this.combinationLimitsForm = lazyFilterForm(() => this.deps.getCombinationLimitInputs());
    this.productFilterForm = lazyFilterForm(() => this.deps.getProductFilterInputs());
    const combinationLimitsPane = limitsAccordion.addPane('Combination limits (optional)',
      () => this.combinationLimitsForm.getRoot(), false, null, false);
    const productFilterPane = limitsAccordion.addPane('Product filters (optional)',
      () => this.productFilterForm.getRoot(), false, null, false);
    this.combinationLimitsDot = attachChangedDot(combinationLimitsPane, 'Changed from platform defaults.');
    this.productFiltersDot = attachChangedDot(productFilterPane, 'Changed from platform defaults.');
    this.accCombinePane = accordion.addPane('How to combine', () => ui.divV([
      ui.divH([
        ui.divText('Strategy', {style: {fontSize: '11px', color: 'var(--grey-6)', marginBottom: '2px'}}),
        this.deps.configInfoIcon,
      ], {style: {alignItems: 'center', gap: '4px'}}),
      ui.divV([this.stratDepthCard.root, this.stratBreadthCard.root, this.reagentsModeNote], {style: {gap: '6px'}}),
      ui.div([ui.form([this.deps.numRoundsInput, this.deps.maxComponentsInput, this.deps.maxRoutesInput])],
        {style: {marginLeft: `${CHEM_ENUM_NESTED_ACCORDION_INDENT}px`}}),
      limitsAccordion.root,
      // First pane in the chain — no Back target.
      navRow(null, mkNextBtn(() => this.accReactionsPane, 'Reactions')),
    ], {style: {gap: '8px'}}), false, null, false);
    this.accCombinePane.root.classList.add('chem-enum-combine-pane');
    this.accCombinePane.root.style.setProperty('--chem-enum-label-width', `${CHEM_ENUM_LABEL_WIDTH}px`);
    this.accCombinePane.root.style.setProperty('--chem-enum-editor-width', `${sharedEditorWidth}px`);
    this.accCombinePane.root.style.setProperty('--chem-enum-form-min-width', `${CHEM_ENUM_FORM_MIN_WIDTH}px`);
    // Left panel for Preview — content built once previewPanel exists; this factory itself only runs
    // lazily when the user actually opens the pane, well after previewPanel/runControls are constructed.
    this.accPreviewPane = accordion.addPane('Preview', () => ui.divV([
      ui.divText('Samples a small subset of products.',
        {style: {fontSize: '12px', color: 'var(--grey-6)'}}),
      this.deps.getPreviewRecapCard(),
      // Last pane in the chain — the run action takes the Next slot instead of a target pane.
      navRow(mkBackBtn(() => this.accExtrasPane, 'Extras'), this.deps.getPreviewEnumerateBtnWrap()),
    ], {style: {gap: '10px'}}), false, null, false);
    this.accPanes = [this.accReactionsPane, this.accBbsPane, this.accExtrasPane, this.accCombinePane, this.accPreviewPane];

    // Navigation is chip/button-driven; native click-to-collapse on the header would just empty the
    // panel, so it's disabled at the source.
    for (const p of this.accPanes) {
      const header = p.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header) {
        header.style.pointerEvents = 'none';
        header.style.cursor = 'default';
      }
    }

    // Subtitle spans injected into each pane's native header — updated by applyRibbonState().
    const injectPaneSub = (pane: DG.AccordionPane): HTMLElement => {
      const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      const sub = document.createElement('span');
      sub.className = 'chem-enum-pane-subtitle';
      header?.appendChild(sub);
      return sub;
    };
    this.subReactions = injectPaneSub(this.accReactionsPane);
    this.subBbs = injectPaneSub(this.accBbsPane);
    this.subExtras = injectPaneSub(this.accExtrasPane);
    this.subCombine = injectPaneSub(this.accCombinePane);

    // Must fit CHEM_ENUM_FORM_MIN_WIDTH (the "How to combine" forms' own min-width, set above) plus
    // the padding/accordion chrome between this pane and those forms — otherwise this pane's content
    // overflows it on every load, not just at a narrow splitter, since the forms now refuse to shrink
    // below their own minimum. The +80 buffer covers this pane's own padding-right plus the nested
    // accordion/pane padding.
    this.leftPane = ui.divV([accordion.root],
      {style: {minWidth: `${CHEM_ENUM_FORM_MIN_WIDTH + 80}px`, overflowY: 'auto', overflowX: 'hidden',
        paddingRight: '8px'}});

    // === Config-summary ribbon chips (shown above the right-pane tabs) ===
    // Each chip's dot flags "customized": Reactions/Building blocks/Reagents show it when any round
    // has a custom subset (passed in via applyRibbonState's override flags — DataPanel-owned, not
    // known here); Strategy shows it when Combination limits or Product filters differ from platform
    // defaults. Text lives in a child span so refreshing it doesn't wipe the dot out (chip.textContent
    // replaces all children).
    const cfgChipEl = (tooltip: string): ChipEl => {
      const dot = mkChangedDot(tooltip);
      const textEl = ui.span([], {style: {overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', minWidth: '0'}});
      const root = ui.divH([dot, textEl], {classes: 'chem-enum-chip', style: {alignItems: 'center', gap: '4px'}});
      return {root, textEl, dot};
    };
    this.chipReactionsC = cfgChipEl('One or more rounds use a custom subset of reaction templates.');
    this.chipBbsC = cfgChipEl('One or more rounds use a custom subset of building blocks.');
    this.chipExtrasC = cfgChipEl('One or more rounds use a custom subset of reagents.');
    this.chipCombineC = cfgChipEl('Changed from platform defaults.');
    this.cfgEstEl = ui.divText('');
    this.cfgEstEl.className = 'chem-enum-chip';
  }

  invalidateLimitForms(): void {
    this.combinationLimitsForm.invalidate();
    this.productFilterForm.invalidate();
  }

  refreshStrategyCards(mode: Mode): void {
    const hasReagents = mode === 'reagents';
    this.applyStratCardStyle(this.stratDepthCard, 'depth', mode, !hasReagents);
    this.applyStratCardStyle(this.stratBreadthCard, 'breadth', mode, !hasReagents);
    this.reagentsModeNote.style.display = hasReagents ? '' : 'none';
  }

  private applyStratCardStyle(card: StratCard, mode: string, cur: string, enabled: boolean): void {
    const sel = cur === mode;
    card.root.style.border = sel ? '2px solid var(--blue-2)' : '1px solid var(--grey-3)';
    card.root.style.opacity = enabled ? '1' : '0.5';
    card.root.style.cursor = enabled ? 'pointer' : 'not-allowed';
    card.icon.style.color = sel ? 'var(--blue-2)' : 'var(--grey-5)';
  }

  applyRibbonState(s: RibbonChipState): void {
    const setChip = (chip: ChipEl, text: string, override: boolean): void => {
      chip.textEl.textContent = text;
      chip.dot.style.display = override ? '' : 'none';
    };
    setChip(this.chipReactionsC, s.reactionsText, s.reactionsOverride);
    setChip(this.chipBbsC, s.bbsText, s.bbsOverride);
    this.chipExtrasC.textEl.textContent = s.extrasText;
    this.chipExtrasC.dot.style.display = s.extrasOverride ? '' : 'none';
    this.chipCombineC.textEl.textContent = s.combineChipText;
    this.chipCombineC.dot.style.display = (s.limitsChanged || s.filtersChanged) ? '' : 'none';
    this.cfgEstEl.textContent = s.estimateText;
    this.combinationLimitsDot.style.display = s.limitsChanged ? '' : 'none';
    this.productFiltersDot.style.display = s.filtersChanged ? '' : 'none';
    this.subReactions.textContent = s.reactionsSubtitle;
    this.subBbs.textContent = s.bbsSubtitle;
    this.subExtras.textContent = s.extrasSubtitle;
    this.subCombine.textContent = s.combineSubtitle;
  }
}
