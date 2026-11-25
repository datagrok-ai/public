import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {CELL_H, CELL_W, COL_HEADER_H, CORE_W} from '../analysis/sar-matrix/sar-matrix-ui-common';


/** The tour of the SAR Matrix demos: cards that walk the viewer from its Summary to the Make list. It starts by
 *  itself the first time a demo is opened in this browser; the icon in the viewer's header brings it back. */
export function addSarMatrixTour(viewer: DG.Viewer, demoName: string, hints: SarHint[]): void {
  const shownKey = `chem-demo-tour-shown-${demoName}`;
  const start = (): void => {
    running?.end();
    running = new SarMatrixTour(viewer, hints, shownKey);
  };
  addTourIcon(viewer, start);
  if (localStorage.getItem(shownKey) === null) {
    // Either landing surface: the matrix grid does not exist until the SAR Matrix tab is shown.
    awaitCheck(() => viewer.root.querySelector('.chem-sar-sum-overview, .chem-sar-grid-host') !== null, '', 120000)
      .then(start, () => {});
  }
}

function addTourIcon(viewer: DG.Viewer, start: () => void): void {
  const icon = ui.iconFA('graduation-cap', start, 'Show the tutorial');
  icon.classList.add('chem-sar-tour-icon');
  viewer.root.closest('.tab-host')?.querySelector(':scope > .tab-host-header > .panel-titlebar .grok-font-icon-help')
    ?.before(icon);
}

/** One at a time: both demos' cards point at whatever SAR Matrix viewer is on screen. */
let running: SarMatrixTour | null = null;

class SarMatrixTour {
  private step = 0;
  private card: HTMLElement | null = null;
  private timer: ReturnType<typeof setTimeout> | undefined;
  private readonly subs: Subscription[];

  constructor(private readonly viewer: DG.Viewer, private readonly hints: SarHint[],
    private readonly shownKey: string) {
    // A viewer out of sight has no size, whether its view, its tab or the viewer itself went away. The first
    // event comes as the observing starts.
    this.subs = [ui.onSizeChanged(viewer.root).subscribe(() => {
      if (this.viewer.root.clientWidth === 0)
        this.hide();
      else if (this.card === null)
        this.show();
    }), viewer.onDetached.subscribe(() => this.end())];
  }

  end(): void {
    this.hide();
    this.subs.forEach((sub) => sub.unsubscribe());
    if (running === this)
      running = null;
  }

  private hide(): void {
    clearTimeout(this.timer);
    this.card?.remove();
    this.card = null;
    document.querySelector('.chem-sar-hint-cell')?.remove();
  }

  private next(): void {
    this.step++;
    this.show();
  }

  private show(tries = 0): void {
    this.hide();
    if (this.step >= this.hints.length) {
      this.end();
      return;
    }
    const hint = this.hints[this.step];
    if (tries === 0)
      hint.before?.();
    const anchor = hint.anchor();
    if (!(anchor instanceof HTMLElement)) {
      // A pane the step just opened builds on a timer, so its anchor is waited for a few beats before the
      // step is skipped.
      if (tries < 8)
        this.timer = setTimeout(() => this.show(tries + 1), 250);
      else
        this.next();
      return;
    }
    // A card low in the list leaves too little height beside it, and the popup runs past the bottom of the pane.
    if (anchor.closest('.chem-sar-nav-list'))
      anchor.scrollIntoView({block: 'start'});
    const last = this.step === this.hints.length - 1;
    const content = ui.div([ui.h3(hint.title), ui.divText(hint.text),
      ui.divH([ui.button(last ? 'Close' : 'Next', () => this.next())], {style: {justifyContent: 'flex-end'}})]);
    this.card = ui.hints.addHint(anchor, content, hint.position ?? ui.hints.POSITION.BOTTOM);
    this.card.classList.add('chem-sar-tour-card');
    this.card.querySelector('.fa-times')?.addEventListener('click', () => this.end());
    localStorage.setItem(this.shownKey, 'true');
  }
}

export interface SarHint {
  anchor: () => Element | null;
  position?: `${ui.hints.POSITION}`;
  title: string;
  text: string;
  /** Open the tab or section the anchor lives on: an anchor on a pane that is not open is not in the document. */
  before?: () => void;
}

function sarTab(title: string): void {
  sarTabHeader(title)?.click();
}

/** The fallback anchor for hints about a tab the tour has not opened yet, whose markup does not exist yet. */
function sarTabHeader(title: string): HTMLElement | null {
  const found = Array.from(document.querySelectorAll('.d4-tab-header')).find((e) => e.textContent?.trim() === title);
  return found instanceof HTMLElement ? found : null;
}

/** Open a Summary row or fold while it is closed, which its right-pointing chevron shows, so a step shown
 *  again, or one the reader already opened, stays open. */
function openIfClosed(el: Element | null): void {
  if (el instanceof HTMLElement && el.querySelector(':scope > .chem-sar-card-twisty.fa-chevron-right'))
    el.click();
}

/** Show a trigger's evidence in its slot, unless it is already the one shown. */
function showEvidence(trigger: Element | null): void {
  if (trigger instanceof HTMLElement && !trigger.classList.contains('chem-sar-sum-open'))
    trigger.click();
}

/** Open the first analog worth making, then bring the list into view. 'nearest' scrolls only a container the
 *  list is outside of; 'start' would also scroll the two-column Summary, whose own overflow is hidden. */
function openMaking(): void {
  const block = document.querySelector('.chem-sar-sum-making-block');
  openIfClosed(block?.querySelector('.chem-sar-card') ?? null);
  block?.scrollIntoView({block: 'nearest'});
}

/** A marker over the first data cell: the matrix is a canvas grid, so its cells are not elements to anchor
 *  to. Fixed-positioned and outside the grid so it cannot disturb the layout it is measuring. */
function matrixCellAnchor(): Element | null {
  const host = document.querySelector('.chem-sar-grid-host');
  if (!(host instanceof HTMLElement))
    return null;
  const box = host.getBoundingClientRect();
  const marker = ui.div([], 'chem-sar-hint-cell');
  marker.style.position = 'fixed';
  marker.style.pointerEvents = 'none';
  marker.style.left = `${box.left + CORE_W}px`;
  marker.style.top = `${box.top + COL_HEADER_H}px`;
  marker.style.width = `${CELL_W}px`;
  marker.style.height = `${CELL_H}px`;
  document.body.appendChild(marker);
  return marker;
}

/** Where the two tours point. Both walk the same analysis, so a selector belongs in one place: the
 *  fallback is what the step lands on while its pane is still building. */
const at = {
  scale: () => document.querySelector('.chem-sar-sum-scale') ?? sarTabHeader('Summary'),
  coverage: () => document.querySelector('.chem-sar-sum-coverage') ??
    document.querySelector('.chem-sar-sum-overview'),
  components: () => document.querySelector('.chem-sar-sum-comp') ?? document.querySelector('.chem-sar-sum-answers'),
  swapRow: () => document.querySelector('.chem-sar-sum-swap') ?? document.querySelector('.chem-sar-sum-answers'),
  making: () => document.querySelector('.chem-sar-sum-making-block') ?? document.querySelector('.chem-sar-sum-making'),
  firstSeries: () => document.querySelector('.chem-sar-sum-list > .chem-sar-card') ??
    document.querySelector('.chem-sar-sum-overview'),
  analog: () => document.querySelector('.chem-sar-sum-making-block .chem-sar-card') ?? at.making(),
  navCard: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card') ??
    sarTabHeader('SAR Matrix'),
  chips: () => document.querySelector('.chem-sar-chips'),
  contextPanel: () => document.querySelector('.grok-prop-panel'),
  transferTab: () => sarTabHeader('SAR Transfer'),
  transferSources: () => document.querySelector('.chem-sar-xfer-nav .chem-sar-scaffold-card') ??
    document.querySelector('.chem-sar-xfer-list') ?? sarTabHeader('SAR Transfer'),
  transferGrid: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-grid-host') ??
    document.querySelector('.chem-sar-xfer-panel') ?? sarTabHeader('SAR Transfer'),
  makeList: () => sarTabHeader('Make list'),
};

export const SAR_HINTS: SarHint[] = [
  {
    before: () => sarTab('Summary'),
    anchor: at.scale,
    title: 'Start on the Summary',
    text: 'What the analysis concluded. This line is the scale everything rests on: the endpoint, its ' +
      'transform, which direction is better and the observed range.',
  },
  {
    anchor: at.coverage,
    title: 'What reached a matrix',
    text: 'How much of the table reached a series: measured, never assayed, or assayed with no analog ' +
      'to pair with.',
  },
  {
    before: () => openIfClosed(at.firstSeries()),
    anchor: at.firstSeries,
    position: ui.hints.POSITION.TOP,
    title: 'Key series',
    text: 'The series to open first. The squares say why: most compounds, widest range, best ' +
      'compound, most predictions, best-validated fit. A row expands for detail, as the first one here; ' +
      'its grid icon opens the matrix.',
  },
  {
    anchor: at.firstSeries,
    position: ui.hints.POSITION.RIGHT,
    title: 'Reading a row',
    text: 'The series and its level, with badges for whether its fit was checked and what its best ' +
      'prediction would gain. Below: compounds measured, their range and the typical value. On the right, ' +
      'the best one measured. Hover a row for the full counts.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-sum-levels') ?? document.querySelector('.chem-sar-sum-overview'),
    position: ui.hints.POSITION.LEFT,
    title: 'Fold levels',
    text: 'Read the whole Summary at one level: L1 is the finest cut, each level up one cut coarser. All ' +
      'counts every series.',
  },
  {
    before: () => document.querySelector('.chem-sar-sum-answers')?.scrollIntoView({block: 'nearest'}),
    anchor: () => document.querySelector('.chem-sar-sum-answers'),
    position: ui.hints.POSITION.TOP,
    title: 'Three answers across all series',
    text: 'The best substituent at each position, the best core, and the swap that raised potency most ' +
      'among measured pairs. Each says what it rests on, or why it cannot answer.',
  },
  {
    before: () => showEvidence([...document.querySelectorAll('.chem-sar-sum-answer')].pop() ?? null),
    anchor: () => document.querySelector('.chem-sar-sum-answers') ?? document.querySelector('.chem-sar-sum-comp'),
    position: ui.hints.POSITION.TOP,
    title: 'The evidence, in place',
    text: 'Each tile opens the evidence behind its answer right underneath, as this one shows its measured ' +
      'swaps; click it again to close it.',
  },
  {
    anchor: at.making,
    position: ui.hints.POSITION.TOP,
    title: 'Worth making',
    text: 'Analogs not in the table, ranked on their gain over the best compound their series made. The ' +
      'info icon says how to read them. A long list shows 20 and folds the rest; Add all puts every one on ' +
      'the Make list. Compounds you hold but never assayed, if any, sit just above.',
  },
  {
    before: openMaking,
    anchor: at.analog,
    position: ui.hints.POSITION.LEFT,
    title: 'One analog',
    text: 'Its series and cell, n = the measured compounds behind the prediction, its gain and the ' +
      'predicted value. Open, as here, it shows the core and substituent it adds up from. The cart adds ' +
      'it to the Make list; the grid icon opens its matrix.',
  },
  {
    before: () => sarTab('SAR Matrix'),
    anchor: at.navCard,
    position: ui.hints.POSITION.RIGHT,
    title: 'The series themselves',
    text: 'The evidence behind the Summary: the series, their matrices, the transfers and the Make ' +
      'list.',
  },
  {
    // The first card, not the list: the list runs the full height of the pane, so a popup beside it
    // is centred halfway down and lands over the matrix rather than beside what it describes.
    anchor: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Analog series',
    text: 'Each card is an analog series: compounds that share a core and vary at one position. It ' +
      'shows the core, its counts and its score.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-scaffold-card') ??
      document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'One scaffold, several series',
    text: 'A Scaffold row groups the series built on one scaffold. They stay separate because each ' +
      'varies a different position.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .ui-input-root'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank by',
    text: 'Reorders the list by potency, by SAR discontinuity, or by a dominant substituent.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-nav-controls .chem-sar-struct-icon'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Filter the series',
    text: 'Narrows the series by potency, spread, size, level or core substructure.',
  },
  {
    // A card that actually folds: the chevron and the level badge are what the text is about, and a
    // leaf card carries neither. The scaffold row has a chevron but no level, so it is not one either.
    anchor: () => Array.from(document.querySelectorAll(
      '.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card:not(.chem-sar-scaffold-card)'))
      .find((card) => card.querySelector('.chem-sar-card-twisty.fa-chevron-down, ' +
        '.chem-sar-card-twisty.fa-chevron-right')) ??
      document.querySelector('.chem-sar-nav-list:not(.chem-sar-xfer-list) .chem-sar-card'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Select and unfold',
    text: 'Click a card to open its matrix. L2 and L3 group cores that agree one or two cuts deeper ' +
      '("L2·6" holds six L1 matrices); the chevron unfolds them. An orange badge means some ' +
      'compounds sit in groups too thin to show.',
  },
  {
    anchor: at.chips,
    title: 'What the open matrix holds',
    text: 'Compounds, size as cores × substituents, activity range, predicted analogs and compounds ' +
      'with no value. The size follows any filter.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-filter-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Filter the matrix',
    text: 'Narrows the cells by potency, reference points, MW, core or substituent.',
  },
  {
    anchor: matrixCellAnchor,
    position: ui.hints.POSITION.RIGHT,
    title: 'Then click a cell',
    text: 'Cores run down the rows, substituents across. A solid cell is a compound you have, a ' +
      'dashed one an analog nobody has made. Click one to inspect it.',
  },
  {
    anchor: at.contextPanel,
    position: ui.hints.POSITION.LEFT,
    title: 'A compound that exists',
    text: 'Its structure, measured potency and SAR neighbours. A "~" value marks a compound you hold ' +
      'without a measurement: the compound is real, the number estimated.',
  },
  {
    anchor: at.contextPanel,
    position: ui.hints.POSITION.LEFT,
    title: 'A virtual analog',
    text: 'Never made, so its potency is predicted from how its core and substituent do elsewhere, ' +
      'with the compounds behind each term. Add to Make list collects it.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-cart-icon'),
    position: ui.hints.POSITION.LEFT,
    title: 'Add to the Make list',
    text: 'Adds the selected cell, measured or predicted, to the Make list.',
  },
  {
    anchor: at.transferTab,
    title: 'Carry the SAR across scaffolds',
    text: 'Pairs cores whose potency trends run in parallel, so an optimization found on one may ' +
      'carry to the other.',
  },
  {
    anchor: at.transferSources,
    position: ui.hints.POSITION.RIGHT,
    title: 'Sources, nested by series',
    text: 'One card per source core, grouped by series. A source with several targets switches ' +
      'between them in the dropdown above the grid.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-nav .chem-sar-nav-controls .ui-input-root') ??
      sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.RIGHT,
    title: 'Rank the transfers',
    text: 'Order by Potent compounds, Fold match (whether step sizes carry, not just direction) or ' +
      'Analogs gained. The funnel filters on all three.',
  },
  {
    anchor: at.transferGrid,
    position: ui.hints.POSITION.LEFT,
    title: 'Analogs the transfer argues for',
    text: 'Both series over their shared substituents. Dashed "predicted" cells are what the transfer ' +
      'suggests; click one for its breakdown and Add to Make list.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-control-bar .ui-input-root') ??
      document.querySelector('.chem-sar-xfer-panel .chem-sar-main-bar') ?? sarTabHeader('SAR Transfer'),
    title: 'Pick the target core',
    text: 'Switches the target core while the source stays put; the chips follow the target.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-xfer-panel .chem-sar-cart-icon') ??
      document.querySelector('.chem-sar-xfer-panel .chem-sar-main-bar') ?? sarTabHeader('SAR Transfer'),
    position: ui.hints.POSITION.LEFT,
    title: 'Collect from a transfer',
    text: 'Adds the selected cell to the Make list, as on the matrix tab.',
  },
  {
    anchor: at.makeList,
    title: 'The Make list',
    text: 'Everything you collected, with structure, potency, evidence and origin. Status says ' +
      'whether it exists, Method whether its number is measured or predicted.',
  },
];

/** The PROTAC tour: the three parts of every molecule are already columns, so nothing is cut and the
 *  question is which of the three to change. */
export const PROTAC_HINTS: SarHint[] = [
  {
    before: () => sarTab('Summary'),
    anchor: at.scale,
    title: '2 792 degraders, split three ways',
    text: 'Warhead, linker and E3 ligand are columns of their own. The endpoint is predicted ' +
      'solubility in log units, higher is better.',
  },
  {
    anchor: at.coverage,
    title: 'What reached a matrix',
    text: 'How many degraders landed in a series, and how many had no analog to pair with.',
  },
  {
    before: () => openIfClosed(at.firstSeries()),
    anchor: () => document.querySelector('.chem-sar-sum-list'),
    position: ui.hints.POSITION.TOP,
    title: 'Which series to open first',
    text: 'The series to open first; the squares say why each is listed. A row expands for detail, as ' +
      'the first one here; its grid icon opens the matrix.',
  },
  {
    anchor: at.firstSeries,
    position: ui.hints.POSITION.RIGHT,
    title: 'Reading a row',
    text: 'One linker per series, with its level and badges for its fit and its best prediction. Below: ' +
      'degraders measured, their range and the typical value. On the right, the best one measured. Hover ' +
      'a row for the full counts.',
  },
  {
    before: () => at.components()?.scrollIntoView({block: 'nearest'}),
    anchor: at.components,
    position: ui.hints.POSITION.TOP,
    title: 'Which of the three to change',
    text: 'The three components ranked by how far each moves solubility, the other two held fixed.',
  },
  {
    before: () => showEvidence(
      document.querySelector('.chem-sar-sum-comp .chem-sar-sum-half:not(.chem-sar-sum-swap)')),
    anchor: () => document.querySelector('.chem-sar-sum-comp .chem-sar-sum-half') ?? at.components(),
    position: ui.hints.POSITION.TOP,
    title: 'Every linker, warhead and ligand ranked',
    text: 'A component\'s value opens its full ranking, drawn as structures, as the warheads show here; ' +
      'its swap opens the measured pairs. Both open right under the row.',
  },
  {
    before: () => showEvidence(at.swapRow()),
    anchor: at.swapRow,
    position: ui.hints.POSITION.TOP,
    title: 'And what was actually measured',
    text: 'Measured pairs alike in every part but one, per component. No model involved. Opened here: ' +
      'every warhead swap behind the best one.',
  },
  {
    anchor: at.making,
    position: ui.hints.POSITION.TOP,
    title: 'Combinations nobody has built',
    text: 'Warhead, linker and E3 ligand combinations never made, ranked on gain over their series\' best. ' +
      'Listed even when none clears the trust gate, as the line under the title says; the info icon says ' +
      'how to read them. The list scrolls on its own, and Add all puts every one on the Make list.',
  },
  {
    before: openMaking,
    anchor: at.analog,
    position: ui.hints.POSITION.LEFT,
    title: 'One combination',
    text: 'Its linker series and cell, n = the degraders behind the prediction, its gain and the predicted ' +
      'solubility. Open, as here, it shows the core and substituent it adds up from. The cart adds it to ' +
      'the Make list; the grid icon opens its matrix.',
  },
  {
    before: () => sarTab('SAR Matrix'),
    anchor: at.navCard,
    position: ui.hints.POSITION.RIGHT,
    title: 'One linker per series',
    text: 'Each card is one linker: warheads across the columns, the E3 ligand folded into the row. ' +
      'Best and Mean on the right are its solubility.',
  },
  {
    anchor: at.chips,
    title: 'What the open matrix holds',
    text: 'Degraders measured, cells the components could fill, predictions and the series\' R².',
  },
  {
    anchor: () => document.querySelector('.chem-sar-axis-input') ?? at.chips(),
    position: ui.hints.POSITION.LEFT,
    title: 'Pick what runs across',
    text: 'Choose which component fills the columns: warheads now, or E3 ligands. Every matrix is rebuilt; ' +
      'the Summary stays the same.',
  },
  {
    anchor: matrixCellAnchor,
    position: ui.hints.POSITION.RIGHT,
    title: 'Click a cell',
    text: 'One warhead against one E3 ligand on this linker, coloured by solubility.',
  },
  {
    anchor: at.contextPanel,
    position: ui.hints.POSITION.LEFT,
    title: 'A degrader that exists',
    text: 'The whole molecule assembled from its three parts, with its measured value. Nothing here ' +
      'is modelled.',
  },
  {
    anchor: () => document.querySelector('.chem-sar-control-bar .chem-sar-cart-icon') ??
      document.querySelector('.chem-sar-control-bar'),
    position: ui.hints.POSITION.LEFT,
    title: 'Collect what is worth making',
    text: 'A faint cell is a degrader nobody has built. The cart adds the selected cells to the Make ' +
      'list with their predicted values.',
  },
  {
    before: () => sarTab('SAR Transfer'),
    anchor: at.transferTab,
    title: 'Carry the SAR across linkers',
    text: 'Linkers tried with the same warheads whose solubility moves together inform each other.',
  },
  {
    anchor: at.transferSources,
    position: ui.hints.POSITION.RIGHT,
    title: 'Sources, nested by series',
    text: 'Each group is a linker lending its trend.',
  },
  {
    anchor: at.transferGrid,
    position: ui.hints.POSITION.LEFT,
    title: 'What the transfer argues for',
    text: 'Both matrices side by side. A source value over a target hole is an analog the transfer ' +
      'suggests.',
  },
  {
    before: () => sarTab('Make list'),
    anchor: at.makeList,
    title: 'The Make list',
    text: 'Everything collected, with its three parts, predicted value and origin. It exports like ' +
      'any table.',
  },
];
