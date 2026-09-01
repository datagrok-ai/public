import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {computeAllTransfers, Transfer, TransferSide, transferStats, TransferStats}
  from './sar-matrix-transfer';
import {SarMatrix} from './sar-matrix-types';
import {BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, CELL_H, COL_HEADER_H, CORE_BG_ARGB,
  GRID_SCROLLBAR_H, HEADER_ARGB, MatrixCellRef, PaneColumn, PaneGridSlot, PaneRow,
  renderMoleculeOnColor} from './sar-matrix-ui-common';

/** A transfer is identified by the core it starts from, so several targets collapse onto one card. */
function transferSourceKey(t: Transfer): string {
  return `${t.a.matrixIndex}:${t.a.rowIndex}`;
}

/** Ranks the transfer by the potency of the best analog it argues for — the matrix list's "Potent
 *  compounds" question asked of a transfer: which pairing points at the most promising unmade
 *  compound. Correlation is deliberately not a scheme: it takes three values across a whole run
 *  (0.87 / 0.95 / 1.00, most of them 1.00), so it cannot order a list. */
export const XFER_RANK_POTENCY = 'Potent compounds';
/** How faithfully the size of each step carries, not merely its direction. Two cores can agree on
 *  every substituent's rank and still disagree on what the change is worth. */
export const XFER_RANK_FOLD = 'Fold match';
/** How many compounds the transfer actually argues for making. Zero means the trend is real but both
 *  cores have already explored the same R-groups, so it proposes nothing. */
export const XFER_RANK_GAINED = 'Analogs gained';
export const XFER_RANK_SCHEMES = [XFER_RANK_POTENCY, XFER_RANK_FOLD, XFER_RANK_GAINED];

/** What the panel needs back from the viewer that owns it. Deliberately narrow: the panel reads the
 *  assembled matrices and the viewer's activity direction, and borrows its grid and score-box
 *  builders so a transfer pane is painted by exactly the same code as a matrix pane. */
export interface TransferPanelHost {
  readonly matrices: SarMatrix[];
  readonly higherIsBetter: boolean;
  readonly computing: boolean;
  readonly detached: boolean;
  readonly transferSimilarity: number;
  /** Whether this panel's tab is the one on screen — a scan that lands while it is hidden must not
   *  force it open. */
  readonly transferTabActive: boolean;
  releaseSlot(slot: PaneGridSlot): void;
  formatActivity(value: number): string;
  addAnalogsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void;
  cardScoreBox(lines: {value: string, label: string}[], tip: () => string): HTMLElement;
  buildPaneGrid(rows: PaneRow[], columns: PaneColumn[], slot: PaneGridSlot): HTMLElement;
}

/** The SAR Transfer tab: the list of transfers on the left and the selected pair laid out on the
 *  right. Owns its own grid slot, so switching transfers releases only this tab’s grid. */
export class TransferPanel {
  readonly root = ui.divH([], 'chem-sar-xfer-panel');
  /** Sorted by correlation, strongest first. */
  private transfers: Transfer[] = [];
  private transfersComputed = false;
  /** Bumped on every invalidation so a scan that finishes late can tell it is stale. */
  private transferGeneration = 0;
  private transfersComputing = false;
  private transferIndex = 0;
  private readonly transferSlot: PaneGridSlot = {state: null, subs: []};
  private rankScheme: string = XFER_RANK_POTENCY;
  /** Source series folded shut in the list, by matrix index. */
  private readonly collapsed = new Set<number>();
  private navCollapsed = false;
  /** Null when no filter is set; otherwise the floors a transfer must clear to be listed. */
  private filter: {r: number, fold: number, gained: number} | null = null;
  /** Per-render cache: transferStats is not free and ranking asks for every transfer's score. */
  private stats = new Map<Transfer, TransferStats>();

  constructor(private readonly host: TransferPanelHost) {}

  /** The transfers involving a given matrix, strongest first — the matrix pane’s `r` chip reads this. */
  transfersInvolving(matrixIndex: number): Transfer[] {
    return this.transfers.filter((t) => t.a.matrixIndex === matrixIndex || t.b.matrixIndex === matrixIndex);
  }

  /** Repaint the open pane in place, for a selection change that moves no rows. */
  invalidateGrid(): void {
    this.transferSlot.state?.grid.invalidate();
  }

  /** Drop the tab’s grid and its subscriptions; the panel stays usable and rebuilds on next open. */
  release(): void {
    this.host.releaseSlot(this.transferSlot);
  }

  showMessage(text: string): void {
    ui.empty(this.root);
    this.root.appendChild(ui.divText(text, 'chem-sar-empty-note'));
  }

  /** Drop the transfer list and tab content; the matrices it indexed are gone. Recomputes on the
   *  next transfer-tab visit. */
  invalidateTransfers(): void {
    this.transferGeneration++;
    this.host.releaseSlot(this.transferSlot);
    ui.empty(this.root);
    this.transfers = [];
    this.transferIndex = 0;
    this.transfersComputed = false;
  }

  /** Core label: "Series B · Core 1" with the series prefix, or just "Core 1" without it. */
  private coreLabel(side: TransferSide, withSeries: boolean): string {
    const matrix = this.host.matrices[side.matrixIndex];
    const core = matrix.rows[side.rowIndex].label;
    return withSeries ? `${matrix.label} · ${core}` : core;
  }

  /** Name a transfer side, prefixing the series whenever there is more than one matrix. */
  private sideLabel(side: TransferSide): string {
    return this.coreLabel(side, this.host.matrices.length > 1);
  }

  /** Both cores in one line. Two cores of the same series name it once — repeating the series on
   *  each side is the common case and would push the header's chips onto a second row. */
  private transferTitle(transfer: Transfer): string {
    if (transfer.a.matrixIndex !== transfer.b.matrixIndex)
      return `${this.sideLabel(transfer.a)} → ${this.sideLabel(transfer.b)}`;
    const pair = `${this.coreLabel(transfer.a, false)} → ${this.coreLabel(transfer.b, false)}`;
    return this.host.matrices.length > 1 ?
      `${this.host.matrices[transfer.a.matrixIndex].label} · ${pair}` : pair;
  }

  /** A transfer's target label with its R-position; series prefix dropped when the target shares the
   *  source's matrix. */
  private transferTargetLabel(transfer: Transfer): string {
    return `${this.coreLabel(transfer.b, transfer.b.matrixIndex !== transfer.a.matrixIndex)} · ${transfer.a.position}`;
  }

  /** Group transfers by source core (best-correlation order preserved), so a source reaching several
   *  targets collapses into one card whose pane dropdown switches between them. */
  private groupTransfersBySource(transfers: Transfer[]): Transfer[][] {
    const groups = new Map<string, Transfer[]>();
    for (const t of transfers) {
      const key = transferSourceKey(t);
      let group = groups.get(key);
      if (!group) {
        group = [];
        groups.set(key, group);
      }
      group.push(t);
    }
    return [...groups.values()];
  }

  /** Select a transfer and redraw the trend view; only the transfer tab is rebuilt. */
  private selectTransfer(transfer: Transfer): void {
    this.transferIndex = Math.max(0, this.transfers.indexOf(transfer));
    this.renderTransferPanel();
  }

  /** Bring the SAR Transfer tab up to date, computing the (quadratic) transfer list on first open
   *  behind a spinner. No-op when the list is already current. */
  activateTransferTab(): void {
    if (this.transfersComputing)
      return;
    if (this.host.computing) {
      // Matrices are still building; compute() ends in render(), which re-enters here.
      ui.empty(this.root);
      this.root.appendChild(ui.divText('Building SAR matrices...', 'chem-sar-empty-note'));
      return;
    }
    if (this.transfersComputed) {
      // List is current; rebuild the pane only if it's not on screen, else repaint the live grid.
      if (this.root.childElementCount === 0)
        this.renderTransferPanel();
      else
        this.transferSlot.state?.grid.invalidate();
      return;
    }
    this.transfersComputing = true;
    ui.empty(this.root);
    ui.setUpdateIndicator(this.root, true, 'Detecting SAR transfers...');
    // Timeout so the tab switch and indicator paint before the synchronous quadratic pass.
    setTimeout(async () => {
      // A matrix rebuild started meanwhile would make these transfers stale; its render() re-enters.
      if (this.host.detached || this.host.computing) {
        this.transfersComputing = false;
        ui.setUpdateIndicator(this.root, false);
        return;
      }
      // Captured before the await so a tab bounce / threshold change in the window can't publish a stale scan.
      const generation = this.transferGeneration;
      let detected: Transfer[] | null = null;
      try {
        detected = await computeAllTransfers(this.host.matrices, this.host.transferSimilarity,
          this.host.higherIsBetter);
      } catch (e) {
        _package.logger.error('SAR Matrix: transfer detection failed', {error: String(e)});
      }
      this.transfersComputing = false;
      ui.setUpdateIndicator(this.root, false);
      // A bumped generation means these indices no longer address the on-screen matrices. The scan is
      // restarted rather than abandoned: the host was emptied when the generation moved, so returning
      // here leaves the tab blank with nothing running and nothing to say until it is reopened.
      if (this.host.detached)
        return;
      if (generation !== this.transferGeneration) {
        if (this.host.transferTabActive)
          this.activateTransferTab();
        return;
      }
      if (detected === null) {
        ui.empty(this.root);
        this.root.appendChild(ui.divText('SAR transfer detection failed — see the browser console.',
          'chem-sar-empty-note'));
        return;
      }
      this.transfers = detected;
      this.transfersComputed = true;
      this.transferIndex = 0;
      this.renderTransferPanel();
    }, 100);
  }

  /** Cached per render: ranking scores every transfer, and transferStats walks both sides' cells. */
  private statsOf(t: Transfer): TransferStats {
    let s = this.stats.get(t);
    if (s === undefined) {
      s = transferStats(t, this.host.higherIsBetter);
      this.stats.set(t, s);
    }
    return s;
  }

  /** Higher is always better, whichever scheme is chosen, so one comparator orders them all; NaN
   *  marks an unscorable transfer, which sorts last. */
  private score(t: Transfer, scheme: string): number {
    if (scheme === XFER_RANK_FOLD)
      return this.statsOf(t).foldMatch ?? NaN;
    if (scheme === XFER_RANK_GAINED)
      return t.predictedSubstituents.length;
    const best = this.statsOf(t).benefiting;
    // The panel is direction-adjusted already: benefiting is picked as the best under the current
    // direction, so its raw value is comparable across transfers without re-applying it.
    return best === null ? NaN : (this.host.higherIsBetter ? best.value : -best.value);
  }

  private passesFilter(t: Transfer): boolean {
    if (this.filter === null)
      return true;
    const fold = this.statsOf(t).foldMatch;
    return t.correlation >= this.filter.r &&
      (fold === null ? this.filter.fold <= 0 : fold >= this.filter.fold) &&
      t.predictedSubstituents.length >= this.filter.gained;
  }

  /** The listed transfers: filtered, then ordered by the chosen scheme. Unscorable ones sort last
   *  rather than vanishing, so a scheme change never silently drops a transfer from the list. */
  private rankedGroups(): Transfer[][] {
    const groups = this.groupTransfersBySource(this.transfers.filter((t) => this.passesFilter(t)));
    const best = (g: Transfer[]): number =>
      g.reduce((m, t) => {const s = this.score(t, this.rankScheme); return Number.isNaN(s) ? m : Math.max(m, s);},
        Number.NEGATIVE_INFINITY);
    return groups.sort((a, b) => {
      const sa = best(a); const sb = best(b);
      const fa = Number.isFinite(sa); const fb = Number.isFinite(sb);
      if (fa !== fb)
        return fa ? -1 : 1;
      // Equal scores are routine here, so the source label breaks the tie and keeps the order stable.
      return sb !== sa ? sb - sa : this.sideLabel(a[0].a).localeCompare(this.sideLabel(b[0].a));
    });
  }

  /** Rebuild the SAR-transfer tab: the transfer list plus the trend view of the selected one. */
  renderTransferPanel(): void {
    // Release the tab's own grid (not the matrix's) so switching transfers doesn't stack up grids.
    this.host.releaseSlot(this.transferSlot);
    ui.empty(this.root);
    this.stats.clear();
    if (this.transfers.length === 0) {
      this.root.appendChild(ui.divText(this.host.matrices.length === 0 ?
        'No SAR matrices to compare — build matrices first on the SAR Matrix tab.' :
        'No SAR transfers detected: no two cores share enough R-groups, on scaffolds alike enough, ' +
        'with correlated potency trends. Lower "Transfer similarity" to accept less alike scaffolds.',
      'chem-sar-empty-note'));
      return;
    }
    if (this.transferIndex >= this.transfers.length)
      this.transferIndex = 0;
    this.root.appendChild(this.buildTransferNav());
    this.root.appendChild(this.buildTransferPane(this.transfers[this.transferIndex]));
  }

  /** The list panel: rank control, count chip and filter over a list that nests each source series'
   *  cores beneath it, mirroring the matrix navigator. */
  private buildTransferNav(): HTMLElement {
    const rankInput = ui.input.choice('Rank by', {
      value: this.rankScheme,
      items: XFER_RANK_SCHEMES,
      onValueChanged: (value) => {
        this.rankScheme = value!;
        this.renderTransferPanel();
      },
    });
    rankInput.caption = '';
    ui.tooltip.bind(rankInput.input, () => {
      const mark = (s: string): string => s === this.rankScheme ? '▸ ' : '   ';
      return ui.divV([
        ui.divText('Rank by — how the transfer list is ordered:'),
        ui.divText(`${mark(XFER_RANK_POTENCY)}Potent compounds — the best unmade analog the transfer points at.`),
        ui.divText(`${mark(XFER_RANK_FOLD)}Fold match — how far the size of each step carries, not its direction.`),
        ui.divText(`${mark(XFER_RANK_GAINED)}Analogs gained — how many compounds it argues for making.`),
      ], 'chem-sar-rank-tip');
    });

    const shown = this.rankedGroups();
    const listed = shown.reduce((n, g) => n + g.length, 0);
    const sources = new Set(shown.map((g) => g[0].a.matrixIndex)).size;
    const proposing = shown.filter((g) => g.some((t) => t.predictedSubstituents.length > 0)).length;
    const countBadge = ui.divText(`${shown.length} transfer${shown.length === 1 ? '' : 's'}`,
      'chem-sar-nav-badge');
    ui.tooltip.bind(countBadge, () =>
      `${shown.length} source cores across ${sources} series, ${listed} pairings in all. ` +
      `${proposing} propose at least one compound to make; ${shown.length - proposing} argue for nothing, ` +
      'their targets having already explored the same R-groups.');

    const matchCount = ui.divText('', 'chem-sar-nav-matches');
    const idleTip = 'Filter transfers by correlation, fold match and analogs gained';
    const filterIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.filterRoot(), 'chem-sar-struct-filters'), filterIcon, {vertical: true});
    }, idleTip);
    filterIcon.classList.add('chem-sar-struct-icon', 'chem-sar-filter-icon');
    if (this.filter !== null) {
      const total = this.groupTransfersBySource(this.transfers).length;
      matchCount.innerText = `${shown.length} of ${total} match`;
      filterIcon.classList.add('chem-sar-filter-on');
      ui.tooltip.bind(filterIcon, () => `${shown.length} of ${total} shown — filter active. Click to edit or clear.`);
    }

    const list = ui.div([], 'chem-sar-nav-list chem-sar-xfer-list');
    this.fillTransferList(list, shown);
    const header = ui.divV([
      ui.divH([rankInput.root, countBadge, filterIcon], 'chem-sar-nav-controls'),
      matchCount,
    ], 'chem-sar-nav-header');
    const nav = ui.divV([header, list], 'chem-sar-nav chem-sar-xfer-nav');
    const collapseBtn = ui.iconFA('chevron-left', () => {
      this.navCollapsed = !this.navCollapsed;
      nav.classList.toggle('chem-sar-collapsed', this.navCollapsed);
    }, 'Collapse / expand the transfer list');
    collapseBtn.classList.add('chem-sar-nav-collapse');
    nav.appendChild(collapseBtn);
    nav.classList.toggle('chem-sar-collapsed', this.navCollapsed);
    return nav;
  }

  /** Cards grouped under the series their source core belongs to. A series contributing one core gets
   *  no header: a heading over a single item is noise, and the card already names its series. */
  private fillTransferList(list: HTMLElement, groups: Transfer[][]): void {
    const byMatrix = new Map<number, Transfer[][]>();
    for (const g of groups) {
      const key = g[0].a.matrixIndex;
      const bucket = byMatrix.get(key);
      if (bucket === undefined)
        byMatrix.set(key, [g]);
      else
        bucket.push(g);
    }
    for (const [matrixIndex, cores] of byMatrix) {
      if (cores.length < 2) {
        list.appendChild(this.buildTransferCard(cores[0]));
        continue;
      }
      const matrix = this.host.matrices[matrixIndex];
      const open = !this.collapsed.has(matrixIndex);
      const twisty = ui.iconFA(open ? 'chevron-down' : 'chevron-right', (e: MouseEvent) => {
        e.stopPropagation();
        if (open)
          this.collapsed.add(matrixIndex);
        else
          this.collapsed.delete(matrixIndex);
        this.renderTransferPanel();
      });
      twisty.classList.add('chem-sar-card-twisty');
      const core = renderMoleculeOnColor(matrix.rows[0].keySmiles, CARD_CORE_W, CARD_CORE_H, CORE_BG_ARGB);
      core.classList.add('chem-sar-card-core');
      const header = ui.divH([twisty, core, ui.divV([
        ui.divText(matrix.label, 'chem-sar-card-name'),
        ui.divText(`${cores.length} cores transfer`, 'chem-sar-card-desc'),
      ], 'chem-sar-card-body')], 'chem-sar-card chem-sar-scaffold-card');
      list.appendChild(header);
      if (!open)
        continue;
      for (const g of cores) {
        const card = this.buildTransferCard(g);
        card.classList.add('chem-sar-card-nested');
        list.appendChild(card);
      }
    }
  }

  /** The filter popup: a floor per metric, since every one of them is "higher is better". */
  private filterRoot(): HTMLElement {
    const f = this.filter ?? {r: 0, fold: 0, gained: 0};
    const rIn = ui.input.float('Min correlation', {value: f.r});
    const foldIn = ui.input.float('Min fold match', {value: f.fold});
    const gainedIn = ui.input.int('Min analogs gained', {value: f.gained});
    ui.tooltip.bind(gainedIn.root, 'Set to 1 to hide the transfers that argue for nothing to make.');
    const apply = ui.button('Apply', () => {
      const next = {r: rIn.value ?? 0, fold: foldIn.value ?? 0, gained: gainedIn.value ?? 0};
      this.filter = (next.r <= 0 && next.fold <= 0 && next.gained <= 0) ? null : next;
      this.renderTransferPanel();
    });
    const clear = ui.button('Clear', () => {
      this.filter = null;
      this.renderTransferPanel();
    });
    return ui.divV([rIn.root, foldIn.root, gainedIn.root, ui.divH([apply, clear])]);
  }

  /** Every transfer sharing this one's source core — the pane's target-dropdown alternatives. */
  private transferSiblings(transfer: Transfer): Transfer[] {
    const key = transferSourceKey(transfer);
    return this.transfers.filter((t) => transferSourceKey(t) === key);
  }

  /** One selectable card for a group of transfers sharing a source core. */
  private buildTransferCard(group: Transfer[]): HTMLElement {
    const selected = this.transfers[this.transferIndex] ?? null;
    // Selected transfer when it belongs to this group, otherwise the strongest (first).
    const shown = (selected && group.includes(selected)) ? selected : group[0];

    const body = ui.divV([
      ui.divText(this.sideLabel(shown.a), 'chem-sar-card-title'),
    ], 'chem-sar-card-body');
    // Painted up front, not via the navigator's deferred observer, which belongs to the matrix list
    // and would leave these blank on refill.
    const sourceRow = this.host.matrices[shown.a.matrixIndex].rows[shown.a.rowIndex];
    const core = renderMoleculeOnColor(sourceRow.keySmiles, CARD_CORE_W, CARD_CORE_H, CORE_BG_ARGB);
    core.classList.add('chem-sar-card-core');
    // Correlation stays on every card whatever the ranking: it is what makes the pair a transfer at
    // all. The second line is the scheme the list is ordered by, so a card shows why it sits here.
    const stats = this.statsOf(shown);
    const lines = [{value: `r ${shown.correlation.toFixed(2)}`, label: 'transfer'}];
    if (this.rankScheme === XFER_RANK_FOLD)
      lines.push({value: stats.foldMatch === null ? '—' : stats.foldMatch.toFixed(2), label: 'fold'});
    else if (this.rankScheme === XFER_RANK_GAINED)
      lines.push({value: `${shown.predictedSubstituents.length}`, label: 'gained'});
    else {
      lines.push({value: stats.benefiting === null ? '—' : this.host.formatActivity(stats.benefiting.value),
        label: 'best'});
    }
    const scoreBox = this.host.cardScoreBox(lines,
      () => `${this.sideLabel(shown.a)} → ${this.transferTargetLabel(shown)}. Correlation of two ` +
        'different chemotypes’ potency trends across the R-groups they have both explored — the SAR ' +
        'learned on one should carry to the other. 1.00 = perfectly parallel. ' +
        (stats.foldMatch === null ? '' :
          `Fold match ${stats.foldMatch.toFixed(2)}: how far the size of each step carries, not just its ` +
          'direction. ') +
        `${shown.predictedSubstituents.length} analog(s) argued for.`);
    const card = ui.divH([core, body, scoreBox], 'chem-sar-card');
    if (selected && group.includes(selected))
      card.classList.add('selected');
    card.onclick = () => this.selectTransfer(shown);
    return card;
  }

  private transferPaneRow(side: TransferSide, colIdxs: number[]): PaneRow {
    return {matrix: this.host.matrices[side.matrixIndex], rowIndex: side.rowIndex, colIdxs,
      label: this.sideLabel(side), markBest: true};
  }

  /** Compact summary of the open transfer, in the matrix pane's chip style: the numbers that say
   *  whether the pairing is worth acting on, with the full wording on hover. */
  private buildTransferChips(transfer: Transfer): HTMLElement {
    const stats = this.statsOf(transfer);
    const chip = (text: string, tip: string): HTMLElement => {
      const el = ui.divText(text, 'chem-sar-chip-badge');
      ui.tooltip.bind(el, () => tip);
      return el;
    };
    const items = [
      chip(`r ${transfer.correlation.toFixed(2)}`,
        'Correlation of the two cores’ potency trends over the R-groups both have explored. ' +
        '1.00 means they rank the substituents identically.'),
      chip(stats.foldMatch === null ? 'fold —' : `fold ${stats.foldMatch.toFixed(2)}`,
        'How far the size of each step carries, not just its direction. 1.00 means a substituent ' +
        'swap is worth the same on both cores; a low value means the trend transfers but the ' +
        'magnitude does not.'),
      chip(`${transfer.substituents.length} shared`,
        `${transfer.substituents.length} R-groups at ${transfer.a.position} measured on both cores — ` +
        `the evidence the correlation rests on. Scaffold similarity ${transfer.similarity.toFixed(2)}.`),
    ];
    const gained = transfer.predictedSubstituents.length;
    items.push(chip(gained === 0 ? 'nothing to make' : `${gained} to make`, gained === 0 ?
      'Both cores have already explored the same R-groups, so this transfer argues for no new compound.' :
      `${gained} R-group(s) measured on one core and never made on the other — what this transfer ` +
      'argues for.'));
    return ui.divH(items, 'chem-sar-chips');
  }

  /** The cells this transfer argues for making. Each predicted pair has one measured side — the
   *  evidence — and one side where the compound does not exist yet; only the latter is a proposal. */
  private predictedCells(transfer: Transfer): MatrixCellRef[] {
    const out: MatrixCellRef[] = [];
    for (let i = 0; i < transfer.predictedSubstituents.length; i++) {
      for (const [side, ci] of [[transfer.a, transfer.predictedACols[i]],
        [transfer.b, transfer.predictedBCols[i]]] as const) {
        const matrix = this.host.matrices[side.matrixIndex];
        if (matrix.cells[side.rowIndex][ci].kind === 'virtual')
          out.push({matrix, ri: side.rowIndex, ci});
      }
    }
    return out;
  }

  /** The SAR transfer view: two cores whose potency trends run in parallel, as two rows of the same
   *  virtualized grid, each resolving against its own matrix. */
  private buildTransferPane(transfer: Transfer): HTMLElement {
    const from = this.sideLabel(transfer.a);
    // The target is named in the title even when the dropdown repeats it: a source with only one
    // target gets no dropdown, and the pane must still say which core it transfers to.
    const parts: HTMLElement[] = [
      ui.divText(this.transferTitle(transfer), 'chem-sar-main-title'),
      this.buildTransferChips(transfer),
    ];

    // Dropdown to switch targets when this source transfers to more than one. It sits in the title
    // row where the matrix pane keeps its Caption control, so both panes read the same way.
    const controls: HTMLElement[] = [];
    const siblings = this.transferSiblings(transfer);
    if (siblings.length > 1) {
      const targetInput = ui.input.choice('Transfer to', {
        value: this.transferTargetLabel(transfer),
        items: siblings.map((t) => this.transferTargetLabel(t)),
        onValueChanged: (value) => {
          const pick = siblings.find((t) => this.transferTargetLabel(t) === value) ?? transfer;
          this.selectTransfer(pick);
        },
      });
      // Named by its tooltip rather than a label, which keeps the header one row of aligned controls.
      targetInput.caption = '';
      targetInput.setTooltip(`Transfer to — the other cores ${from} carries its SAR to. ` +
        'Switching target keeps the source and re-reads the pairing against the core you pick.');
      controls.push(targetInput.root);
    }

    // Collects what the pane argues for in one action, which is the whole point of opening a
    // transfer; a single cell can still be added from its own Context Panel.
    const addIcon = ui.iconFA('cart-plus', () => this.host.addAnalogsToMakeList(this.predictedCells(transfer),
      'This transfer argues for no new compound — both cores have explored the same R-groups.'));
    addIcon.classList.add('chem-sar-struct-icon', 'chem-sar-cart-icon');
    ui.tooltip.bind(addIcon, () => {
      const count = this.predictedCells(transfer).length;
      return count === 0 ? 'Add to Make list — this transfer argues for no new compound.' :
        `Add the ${count} analog(s) this transfer argues for to the Make list`;
    });
    controls.push(addIcon);
    const controlBar = ui.divH(controls, 'chem-sar-control-bar');
    controlBar.classList.add('chem-sar-control-inline');
    parts.push(controlBar);
    const bar = ui.divH(parts, 'chem-sar-main-bar');

    // Each side keeps its own matrix. Measured pairs first, then the predicted analogs the transfer argues for.
    const rows = [
      this.transferPaneRow(transfer.a, [...transfer.aCols, ...transfer.predictedACols]),
      this.transferPaneRow(transfer.b, [...transfer.bCols, ...transfer.predictedBCols]),
    ];
    const matchedCount = transfer.substituents.length;
    const columns: PaneColumn[] = [...transfer.substituents, ...transfer.predictedSubstituents]
      .map((substSmiles, i) => ({substSmiles, position: transfer.a.position,
        caption: i < matchedCount ? '' : 'predicted'}));
    const gridHost = this.host.buildPaneGrid(rows, columns, this.transferSlot);
    // Two rows only: size to content, not stretch, so the statistics stay on screen.
    gridHost.style.flex = '0 0 auto';
    gridHost.style.height = `${COL_HEADER_H + rows.length * CELL_H + GRID_SCROLLBAR_H}px`;

    // Plain div, not ui.box (which would pin a fixed width).
    const scroll = ui.div([this.buildTransferStats(transfer)], 'chem-sar-main-scroll');
    return ui.divV([bar, gridHost, scroll], 'chem-sar-main');
  }

  /** What the transfer argues for: the single best analog it predicts, and where to make it. The
   *  correlation, fold match and shared-R-group counts are not repeated here — they are the chips. */
  private buildTransferStats(transfer: Transfer): HTMLElement {
    const stats = this.statsOf(transfer);
    const row = (label: string, value: HTMLElement): HTMLElement =>
      ui.divH([ui.divText(label, 'chem-sar-xfer-stat-l'), value], 'chem-sar-xfer-stat');
    const text = (s: string): HTMLElement => ui.divText(s, 'chem-sar-xfer-stat-v');

    // One row per fact rather than one packed line: the pane is often only a couple of hundred
    // pixels wide, where a core label and a depiction on the same line both get squeezed.
    const best = stats.benefiting;
    const rows = best === null ?
      [row('Make', text('nothing — both cores have explored the same R-groups'))] :
      [
        row('Make on', text(this.sideLabel(best.side === 'a' ? transfer.a : transfer.b))),
        row('Substituent',
          renderMoleculeOnColor(best.substSmiles, BENEFIT_MOL_W, BENEFIT_MOL_H, HEADER_ARGB)),
        row('Predicted', text(this.host.formatActivity(best.value))),
      ];
    return ui.divV([ui.divText('Best predicted analog', 'chem-sar-xfer-title'), ...rows],
      'chem-sar-xfer-stats');
  }
}
