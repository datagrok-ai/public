/* The building blocks every Summary segment draws with — rows, badges, depictions, number formats —
   and the rankings more than one segment reads. The panel extends this. */
import * as ui from 'datagrok-api/ui';
import {SarMatrix} from './sar-matrix-types';
import {MIN_SUPPORT, SeriesStat, SummaryData, SummaryHost, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {chipBadge, CORE_BG_ARGB, count, isStructure, MatrixCellRef, paintMoleculeOnColor, STRIP_MOL_H,
  STRIP_MOL_W, tipText} from './sar-matrix-ui-common';

export class SummaryKit {
  /** Depictions waiting for the tick after the pane lands. */
  pendingPaints: (() => void)[] = [];
  /**
   * The fold tier every ranking on this tab is read at, or null for all of them together.
   *
   * Not a cut of the data: the matrices are built once and all of them stay. This only decides which
   * of them the tab's rankings walk — a tier holds the same compounds as the tier below it over cores
   * cut one bond broader, so reading at one tier answers "what does the SAR look like at this breadth"
   * without the other breadths mixed in.
   */
  protected tierFilter: number | null = null;

  constructor(readonly host: SummaryHost) {}

  // ---- Elements ----------------------------------------------------------------------------------------

  /** Every clickable row on the tab. `.chem-sar-card` carries the hover, the left border and the row
   *  metric the two text lines and the value column are all sized against. */
  summaryRow(opts: {
    depiction: HTMLElement | null,
    name: string,
    badges: HTMLElement[],
    desc: string,
    /** A mark lane between the text column and the value, at a fixed width so it reads down a list. */
    mark?: HTMLElement,
    value: string,
    caption: string,
    valueTip?: string,
    faint?: boolean,
    cart?: MatrixCellRef,
    tip: string,
    onClick: () => void,
    /** Makes the chevron a toggle of its own rather than a repeat of the row's click. */
    onChevron?: () => void,
  }): HTMLElement {
    const value = ui.divText(opts.value,
      `chem-sar-sum-value${opts.faint ? ' chem-sar-sum-faint' : ''}`);
    const caption = ui.divText(opts.caption, 'chem-sar-sum-cap');
    const valueBox = ui.divV([value, caption], 'chem-sar-sum-valuebox');
    if (opts.valueTip)
      ui.tooltip.bind(valueBox, () => opts.valueTip!);
    const parts: HTMLElement[] = [];
    if (opts.depiction !== null)
      parts.push(opts.depiction);
    parts.push(ui.divV([
      ui.divH([ui.divText(opts.name, 'chem-sar-card-name'), ...opts.badges], 'chem-sar-card-title'),
      ui.divText(opts.desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body'));
    if (opts.mark !== undefined)
      parts.push(opts.mark);
    parts.push(valueBox);
    if (opts.cart !== undefined) {
      const ref = opts.cart;
      const cart = ui.iconFA('cart-plus', (e: MouseEvent) => {
        e.stopPropagation();
        this.host.addCellsToMakeList([ref], 'This cell has no structure to add.');
      }, 'Add this compound to the Make list');
      cart.classList.add('chem-sar-struct-icon', 'chem-sar-cart-icon');
      parts.push(cart);
    }
    const toggle = opts.onChevron;
    const go = toggle === undefined ? ui.iconFA('chevron-right') :
      ui.iconFA('chevron-down', (e: MouseEvent) => {
        e.stopPropagation();
        toggle();
      }, 'What is worth opening this series for');
    go.classList.add('chem-sar-sum-go');
    parts.push(go);
    const row = ui.divH(parts, 'chem-sar-card chem-sar-sum-row');
    ui.tooltip.bind(row, () => opts.tip);
    row.onclick = opts.onClick;
    return row;
  }

  /** Canvas now, RDKit later: a synchronous pass over every structure on the tab would stall the tab
   *  switch that asked for them. */
  depiction(smiles: string | null, w: number, h: number): HTMLElement {
    const canvas = ui.canvas(w, h);
    canvas.classList.add('chem-sar-card-core');
    if (smiles) {
      this.pendingPaints.push(() => {
        if (isStructure(smiles)) {
          paintMoleculeOnColor(canvas, smiles, w, h, CORE_BG_ARGB);
          return;
        }
        // A component given by name, such as "VHL", has nothing to draw and is written out instead.
        const name = ui.divText(smiles, 'chem-sar-cp-frag-name chem-sar-sum-name-art');
        name.style.width = `${w}px`;
        name.style.height = `${h}px`;
        canvas.replaceWith(name);
      });
    }
    return canvas;
  }

  flushPaints(): void {
    const paints = this.pendingPaints;
    this.pendingPaints = [];
    for (const paint of paints)
      paint();
  }

  badge(text: string, tip: string, partial = false): HTMLElement {
    return chipBadge(text, tip, partial ? 'chem-sar-chip-partial' : '');
  }

  /** Whether this series' fit was checked, and whether it held. No colour — `confidence` is null
   *  where too few cells could be cross-validated or every measured value was identical, and a traffic
   *  light would make "unchecked" read as a failing grade. */
  trustDot(matrix: SarMatrix): HTMLElement {
    const conf = matrix.confidence;
    const dot = ui.div([], 'chem-sar-sum-dot' + (!conf ? ' chem-sar-sum-dot-none' :
      conf.r2 < TRUST_R2 ? ' chem-sar-sum-dot-open' : ''));
    ui.tooltip.bind(dot, () => !conf ?
      'Unchecked: this series has fewer than four measured cells that could be held out, or every ' +
      'measured value in it is identical. Its predictions are unverified, not wrong.' :
      `R² ${conf.r2.toFixed(2)} ± ${this.host.formatActivity(conf.rmse)} predicting cells held out of ` +
      `the fit, over ${conf.n} of ${conf.total} measured — the fit ` +
      `${conf.r2 >= TRUST_R2 ? 'holds' : 'does not hold'} at R² ${TRUST_R2}. A property of the model, ` +
      'not of the data.');
    return dot;
  }

  /** The substituent a tile's conclusion is about, drawn. Undefined where there is no one structure —
   *  a tile that names two groups or declines to rank has nothing to draw. */
  answerArt(smiles: string | null): HTMLElement | undefined {
    return smiles ? this.depiction(smiles, STRIP_MOL_W, STRIP_MOL_H) : undefined;
  }

  /** A statement on the pane and its qualification on the hover. Method is read for the one fact it is
   *  opened for, and a paragraph per fact buries that fact in the others. */
  prose(text: string, tip: string): HTMLElement {
    return tipText(text, 'chem-sar-sum-prose', tip);
  }

  /** The same, in the faint style a row uses for its trailing clause. */
  hint(text: string, tip: string): HTMLElement {
    return tipText(text, 'chem-sar-cp-hint', tip);
  }

  /** A card with no qualifying rows still renders: a vanishing card reflows the grid and hides the
   *  fact that the analysis found nothing of that kind. */
  reason(text: string): HTMLElement[] {
    return [ui.divText(text, 'chem-sar-cp-hint')];
  }

  /** A chevron heading that shows and hides `body`; the caller owns whatever else sits on the row. */
  foldHead(title: string, body: HTMLElement, open: boolean, cls = '',
    onToggle?: (open: boolean) => void): HTMLElement {
    const icon = ui.iconFA('chevron-right');
    icon.classList.add('chem-sar-sum-fold-icon');
    icon.classList.toggle('chem-sar-sum-fold-open', open);
    body.style.display = open ? '' : 'none';
    const head = ui.divH([icon, ui.divText(title, 'chem-sar-sum-fold-title')],
      `chem-sar-sum-fold ${cls}`.trim());
    head.onclick = () => {
      const show = body.style.display === 'none';
      body.style.display = show ? '' : 'none';
      icon.classList.toggle('chem-sar-sum-fold-open', show);
      onToggle?.(show);
    };
    return head;
  }

  /** A body filled on its first opening, so an expand nobody opens draws nothing. */
  lazyBody(fill: () => HTMLElement[]): {body: HTMLElement, open: () => void, toggle: () => void} {
    const body = ui.divV([]);
    body.style.display = 'none';
    let filled = false;
    const open = (): void => {
      if (!filled) {
        filled = true;
        body.append(...fill());
        // Rows added after the paint tick has run, so these depictions need one of their own.
        this.flushPaints();
      }
      body.style.display = '';
    };
    const toggle = (): void => {
      if (body.style.display === 'none')
        open();
      else
        body.style.display = 'none';
    };
    return {body, open, toggle};
  }

  anchored(el: HTMLElement, anchor: string): HTMLElement {
    el.dataset.anchor = anchor;
    return el;
  }

  // ---- Text --------------------------------------------------------------------------------------------

  /** Signed, and finer than `formatActivity`: fitted effects and pooled differences live in tenths,
   *  where one decimal renders three differently-ranked R-groups as three rows all reading "+0.1". */
  formatEffect(value: number): string {
    const abs = Math.abs(value);
    return `${value < 0 ? '−' : '+'}${abs < 10 ? abs.toFixed(2) : abs.toFixed(1)}`;
  }

  /** A measured difference in the reader's own units: a signed log difference where the scale is a
   *  log, and the fold it corresponds to where it is not. A loss is written as its reciprocal fold:
   *  10^value rounds every loss past one log to "0.0×", which destroys the number it is reporting. */
  formatDelta(value: number): string {
    if (this.host.activityIsLog)
      return this.formatEffect(value);
    const fold = Math.pow(10, Math.abs(value));
    const shown = fold < 10 ? fold.toFixed(1) : fold.toFixed(0);
    return value < 0 ? `1/${shown}×` : `${shown}×`;
  }

  cellLocation(matrix: SarMatrix, ri: number, ci: number): string {
    return `${matrix.label} · ${matrix.rows[ri].label} × ${matrix.columns[ci].position}`;
  }

  openTip(matrix: SarMatrix, ri: number, ci: number): string {
    return `Open ${this.cellLocation(matrix, ri, ci)} in the SAR Matrix`;
  }

  // ---- Readings more than one segment shares -----------------------------------------------------------

  /** Whether the front-runner may be named alone: its margin over the runner-up has to clear the
   *  screen's own model error, or the two are one answer with two names. */
  leads(first: number | null, second: number | null, modelError: number | null): boolean {
    if (second === null)
      return true;
    if (first === null)
      return false;
    return first - second > (modelError ?? 0);
  }

  /** Why the whole fit cannot be read, or null when it can. */
  roleFitRefusal(data: SummaryData): string | null {
    const fit = data.roleFit;
    if (fit === null) {
      return `No component column carries two values with ${MIN_SUPPORT} or more compounds each, so ` +
        'there is nothing to compare.';
    }
    if (!fit.converged) {
      return `The fit did not settle within ${fit.sweeps} passes. Each component's offsets are ` +
        'still shrunk toward zero by a different amount and cannot be compared, so nothing is ranked.';
    }
    if (fit.cvR2 === null) {
      return 'Too few compounds to cross-validate a fit over ' +
        `${count(fit.roles.reduce((sum, role) => sum + role.fitted, 0))} component values. Unchecked, ` +
        'not wrong — nothing is ranked.';
    }
    if (fit.cvR2 < TRUST_R2) {
      return `Cross-validated R² ${fit.cvR2.toFixed(2)} — below ${TRUST_R2}, the additive reading does ` +
        'not hold on this table: knowing every component does not predict a compound the fit had not ' +
        'seen. The offsets are not ranked.';
    }
    return null;
  }

  /** Series a core leaderboard may rank at all: a matrix has to BE one core, since elsewhere a core is
   *  one row of one matrix and its fitted mean compares groupings, not scaffolds. */
  rankableCores(data: SummaryData): SeriesStat[] {
    return data.coresAreSeries ?
      data.series.filter((stat) => stat.realCells > 0 && stat.cpd >= MIN_SUPPORT) : [];
  }

  /**
   * One tier's cores, best first.
   *
   * Never two tiers at once: a folded tier holds the same compounds over a core cut one bond broader,
   * so a mixed ranking enters one compound twice under two cores and reads as two scaffolds agreeing.
   */
  topCores(data: SummaryData): SeriesStat[] {
    const tier = this.tierFilter ?? this.defaultCoreTier(data);
    if (tier === null)
      return [];
    const dir = this.host.higherIsBetter ? 1 : -1;
    return this.rankableCores(data).filter((stat) => stat.tier === tier)
      .sort((a, b) => dir * (b.typical - a.typical) || (a.matrix.id < b.matrix.id ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  /** The tier ranked when the reader has not picked one: the one holding the most cores, which is the
   *  cut depth this library actually recurs at. */
  private defaultCoreTier(data: SummaryData): number | null {
    const tiers = this.coreTiers(data);
    if (tiers.length === 0)
      return null;
    const counts = new Map(tiers.map((tier) =>
      [tier, this.rankableCores(data).filter((stat) => stat.tier === tier).length]));
    return tiers.reduce((best, tier) => counts.get(tier)! > counts.get(best)! ? tier : best, tiers[0]);
  }

  /** The fold tiers that carry rankable cores, shallowest first. */
  coreTiers(data: SummaryData): number[] {
    return [...new Set(this.rankableCores(data).map((stat) => stat.tier))].sort((a, b) => a - b);
  }
}
