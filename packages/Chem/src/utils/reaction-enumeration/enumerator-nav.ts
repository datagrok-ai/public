/* eslint-disable max-len */
import {Subscription} from 'rxjs';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {CHANGED_DOT_STYLE, Mode} from './enumerator-app';

export type ChipEl = {root: HTMLElement; textEl: HTMLElement; dot: HTMLElement};
type StratCard = {root: HTMLElement; icon: HTMLElement};

// Plain-data snapshot of everything the ribbon chips + accordion-pane subtitles need to show —
// computed by EnumeratorConfigForm (which owns EnumeratorConfig) and handed in wholesale, so this
// class never needs to read EnumeratorConfig itself.
export interface RibbonChipState {
  reactionsText: string; reactionsErr: boolean; reactionsOverride: boolean; reactionsSubtitle: string;
  bbsText: string; bbsErr: boolean; bbsOverride: boolean; bbsSubtitle: string;
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
      sub.add(ui.onSizeChanged(this.accCombinePane.root).subscribe(apply));
      setTimeout(apply, 0);
      this.deps.view.subs.push(sub);
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
        indentSub = flushIndent(form, this.deps.numRoundsInput.root);
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
      this.deps.numRoundsInput, this.deps.maxComponentsInput, this.deps.maxRoutesInput,
      ...this.deps.getCombinationLimitInputs(), ...this.deps.getProductFilterInputs(),
    ];
    const sharedLabelWidth = Math.ceil(Math.max(...measureLabelWidths(sectionInputs.map((inp) => inp.captionLabel))));
    // Editors need pinning too, not just labels: the platform widens an editor from its normal ~150px
    // to ~176px on whichever of the three forms it currently considers .ui-form-condensed, so without
    // this the Product filters column can end up visibly wider than Combination limits' even with
    // labels aligned. 150px matches the platform's own un-condensed default.
    const sharedEditorWidth = 150;
    // Mirrors the per-input-type minimum widths the platform itself uses to decide .ui-form-condensed
    // (js-api ui.ts's getInputsMinWidths). Sizing the form to fit the widest one up front, the same
    // way the platform already does for its own dialog forms (a min-width computed from label + input
    // minimums instead of reacting to condensed after the fact), means these forms never need
    // condensed layout at all. Named constants (not inlined) so a future change to those platform
    // thresholds is easier to spot and re-sync here.
    const PLATFORM_MIN_WIDTH_TEXT_TABLE = 200;
    const PLATFORM_MIN_WIDTH_FLOAT_DATE = 140;
    const PLATFORM_MIN_WIDTH_DEFAULT = 100;
    const platformInputMinWidth = (input: DG.InputBase<unknown>): number => {
      const el = input.root;
      if (el.classList.contains('ui-input-text') || el.classList.contains('ui-input-table')) return PLATFORM_MIN_WIDTH_TEXT_TABLE;
      if (el.classList.contains('ui-input-float') || el.classList.contains('ui-input-date')) return PLATFORM_MIN_WIDTH_FLOAT_DATE;
      return PLATFORM_MIN_WIDTH_DEFAULT;
    };
    const formMinInputWidth = Math.ceil(Math.max(...sectionInputs.map(platformInputMinWidth)));
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
      ui.form([this.deps.numRoundsInput, this.deps.maxComponentsInput, this.deps.maxRoutesInput]),
      limitsAccordion.root,
      // First pane in the chain — no Back target.
      navRow(null, mkNextBtn(() => this.accReactionsPane, 'Reactions')),
    ], {style: {gap: '8px'}}), false, null, false);
    this.accCombinePane.root.classList.add('chem-enum-combine-pane');
    this.accCombinePane.root.style.setProperty('--chem-enum-label-width', `${sharedLabelWidth}px`);
    this.accCombinePane.root.style.setProperty('--chem-enum-editor-width', `${sharedEditorWidth}px`);
    // +40 matches the buffer the platform's own d4-dialog-contents branch adds on top of
    // label + input minimums (js-api ui.ts handleFormResize).
    this.accCombinePane.root.style.setProperty('--chem-enum-form-min-width', `${sharedLabelWidth + formMinInputWidth + 40}px`);
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

    // Must fit --chem-enum-form-min-width (the "How to combine" forms' own min-width, computed above)
    // plus the padding/accordion chrome between this pane and those forms — otherwise this pane's
    // content overflows it on every load, not just at a narrow splitter, since the forms now refuse to
    // shrink below their own minimum. The +80 buffer covers this pane's own padding-right plus the
    // nested accordion/pane padding.
    this.leftPane = ui.divV([accordion.root],
      {style: {minWidth: `${sharedLabelWidth + formMinInputWidth + 80}px`, overflowY: 'auto', overflowX: 'hidden',
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
    const setChip = (chip: ChipEl, text: string, err: boolean, override: boolean): void => {
      chip.textEl.textContent = text;
      chip.root.classList.toggle('chem-enum-chip--err', err);
      chip.dot.style.display = override ? '' : 'none';
    };
    setChip(this.chipReactionsC, s.reactionsText, s.reactionsErr, s.reactionsOverride);
    setChip(this.chipBbsC, s.bbsText, s.bbsErr, s.bbsOverride);
    // Extras is fully optional — never flagged as an error state.
    this.chipExtrasC.textEl.textContent = s.extrasText;
    this.chipExtrasC.dot.style.display = s.extrasOverride ? '' : 'none';
    // "Strategy:" prefix only on the ribbon chip — the accordion pane itself already says "How to combine".
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
