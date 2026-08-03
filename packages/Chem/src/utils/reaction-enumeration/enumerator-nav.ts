/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {CHANGED_DOT_STYLE, Mode} from './shared';

export type ChipEl = {root: HTMLElement; textEl: HTMLElement; dot: HTMLElement};
type StratCard = {root: HTMLElement; icon: HTMLElement};

/** Everything the ribbon chips and pane subtitles show, computed by EnumeratorConfigForm and
 * handed in wholesale so this class never reads EnumeratorConfig itself. */
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
  // Getters, not snapshots: both field groups are rebuilt wholesale on every YAML load, so the
  // lazy form builder must re-read them rather than capture the arrays present at construction.
  getCombinationLimitInputs: () => DG.InputBase<unknown>[];
  getProductFilterInputs: () => DG.InputBase<unknown>[];
  setAndFire: <T>(input: DG.InputBase<T>, v: T) => void;
  openAccPaneAndSyncTab: (pane: DG.AccordionPane) => void;
  getPreviewRecapCard: () => HTMLElement;
  getPreviewEnumerateBtnWrap: () => HTMLElement;
}

/** Strategy picker cards, the left-nav accordion (5 panes forming one navigational ring), and the
 * ribbon chips — two surfaces for the same "which section is active / what does it hold" state,
 * kept in sync by one applyRibbonState() call. Pure layout: it reads no config of its own. */
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
    // Cards replace the depth-first checkbox and drive the (hidden) depthFirstInput. Both are
    // disabled while a reagents file is loaded, since reagents mode makes the strategy irrelevant.
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

    // Target pane resolved via a thunk: the Reactions pane is added expanded, so its factory runs
    // inside addPane before the later panes exist.
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
    const navRow = (back: HTMLElement | null, next: HTMLElement | null): HTMLElement =>
      ui.divH([back ?? ui.div([]), next ?? ui.div([])], {classes: 'chem-enum-nav-row'});

    // Last addPane arg is allowDragOut, which defaults to true. One shared form per pane so the
    // field labels align (separate forms size independently).
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
    const attachChangedDot = (pane: DG.AccordionPane, tooltip: string): HTMLElement => {
      const dot = mkChangedDot(tooltip);
      dot.style.marginLeft = '6px';
      pane.root.querySelector('.d4-accordion-pane-header')?.appendChild(dot);
      return dot;
    };
    // Built lazily on first pane expand: a form built while detached measures 0 width and gets
    // marked .ui-form-condensed. invalidate() rebuilds it after a YAML load swaps in fresh inputs.
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
    // "How to combine" holds three separate forms, which size their label/editor columns
    // independently. These pin all three to one shared column so they line up, and keep the forms
    // wide enough that the platform never marks them .ui-form-condensed. Applied as CSS custom
    // properties, since the platform's async label auto-sizing overwrites plain inline styles.
    const CHEM_ENUM_LABEL_WIDTH = 140;
    const CHEM_ENUM_EDITOR_WIDTH = 150;
    const CHEM_ENUM_FORM_MIN_WIDTH = 280;
    // The limits sub-accordion sits one level deeper, which the platform indents by 6px padding +
    // 20px margin; match it so the rounds form's label column lines up with the ones below it.
    const CHEM_ENUM_NESTED_ACCORDION_INDENT = 26;
    // Independently collapsible, unlike the outer navigation accordion.
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
    this.accCombinePane.root.style.setProperty('--chem-enum-editor-width', `${CHEM_ENUM_EDITOR_WIDTH}px`);
    this.accCombinePane.root.style.setProperty('--chem-enum-form-min-width', `${CHEM_ENUM_FORM_MIN_WIDTH}px`);
    // This factory only runs when the user opens the pane, well after previewPanel/runControls exist.
    this.accPreviewPane = accordion.addPane('Preview', () => ui.divV([
      ui.divText('Samples a small subset of products.',
        {style: {fontSize: '12px', color: 'var(--grey-6)'}}),
      this.deps.getPreviewRecapCard(),
      // Last pane in the chain — the run action takes the Next slot instead of a target pane.
      navRow(mkBackBtn(() => this.accExtrasPane, 'Extras'), this.deps.getPreviewEnumerateBtnWrap()),
    ], {style: {gap: '10px'}}), false, null, false);
    this.accPanes = [this.accReactionsPane, this.accBbsPane, this.accExtrasPane, this.accCombinePane, this.accPreviewPane];

    // Navigation is chip/button-driven; native click-to-collapse would just empty the panel.
    for (const p of this.accPanes) {
      const header = p.root.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
      if (header) {
        header.style.pointerEvents = 'none';
        header.style.cursor = 'default';
      }
    }

    // Subtitle spans injected into each pane header, updated by applyRibbonState().
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

    // The forms above refuse to shrink below CHEM_ENUM_FORM_MIN_WIDTH, so this pane must fit that
    // plus the accordion chrome around them, or its content overflows on every load.
    this.leftPane = ui.divV([accordion.root],
      {style: {minWidth: `${CHEM_ENUM_FORM_MIN_WIDTH + 80}px`, overflowY: 'auto', overflowX: 'hidden',
        paddingRight: '8px'}});

    // Text lives in a child span so refreshing it doesn't wipe out the dot.
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
