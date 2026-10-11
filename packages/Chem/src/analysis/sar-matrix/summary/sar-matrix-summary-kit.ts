/* The building blocks every Summary segment draws with — rows, badges, depictions, number formats —
   and the rankings more than one segment reads. The panel extends this. */
import * as ui from 'datagrok-api/ui';
import {SarMatrix} from '../sar-matrix-types';
import {MIN_SUPPORT, SeriesStat, SummaryData, SummaryHost, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {chipBadge, CORE_BG_ARGB, isStructure, MatrixCellRef, paintMoleculeOnColor, STRIP_MOL_H,
  STRIP_MOL_W, tipText} from '../sar-matrix-ui-common';

export class SummaryKit {
  /** Depictions waiting for the tick after the pane lands. */
  pendingPaints: (() => void)[] = [];
  /** The fold tier every ranking on the tab reads, or null for all. It picks which matrices are walked;
   *  nothing is rebuilt. */
  protected tierFilter: number | null = null;

  constructor(readonly host: SummaryHost) {}

  // ---- Elements ----------------------------------------------------------------------------------------

  /** Every clickable row on the tab, styled as a navigator card. */
  summaryRow(opts: {
    depiction: HTMLElement | null,
    name: string,
    badges: HTMLElement[],
    /** Drawn under the description, where it lines up down the list and never squeezes the text. */
    mark?: HTMLElement,
    desc: string,
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
      `chem-sar-card-title${opts.faint ? ' chem-sar-sum-faint' : ''}`);
    const caption = ui.divText(opts.caption, 'chem-sar-card-desc');
    const valueBox = ui.divV([value, caption], 'chem-sar-sum-valuebox');
    if (opts.valueTip)
      ui.tooltip.bind(valueBox, () => opts.valueTip!);
    const parts: HTMLElement[] = [];
    if (opts.depiction !== null)
      parts.push(opts.depiction);
    const body = ui.divV([
      ui.divH([ui.divText(opts.name, 'chem-sar-card-name'), ...opts.badges], 'chem-sar-card-title'),
      ui.divText(opts.desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body');
    if (opts.mark !== undefined)
      body.appendChild(opts.mark);
    parts.push(body, valueBox);
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
    go.classList.add('chem-sar-card-twisty');
    parts.push(go);
    const row = ui.divH(parts, 'chem-sar-card');
    ui.tooltip.bind(row, () => opts.tip);
    row.onclick = opts.onClick;
    return row;
  }

  /** Canvas now, RDKit on the next tick, so a tab switch is not held up by drawing. */
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
        const name = ui.divText(smiles, 'chem-sar-cp-frag-name');
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

  /** Whether this series' fit was checked and held. No confidence means unchecked, not failing. */
  trustBadge(matrix: SarMatrix): HTMLElement {
    const conf = matrix.confidence;
    return this.badge(!conf ? 'unchecked' : `R² ${conf.r2.toFixed(2)}`, !conf ?
      'Too few cells to cross-validate, or all values identical: unverified, not wrong.' :
      `± ${this.host.formatActivity(conf.rmse)} on held-out cells (${conf.n} of ${conf.total}); the fit ` +
      `${conf.r2 >= TRUST_R2 ? 'holds' : 'does not hold'} at R² ${TRUST_R2}.`, !!conf && conf.r2 < TRUST_R2);
  }

  /** Two structures and the glyph between them: the shape every swap is drawn in. */
  pairArt(from: string, to: string, cls: string, glyph = '→', w = STRIP_MOL_W, h = STRIP_MOL_H): HTMLElement {
    return ui.divH([this.depiction(from, w, h), ui.divText(glyph, 'chem-sar-sum-arrow'), this.depiction(to, w, h)],
      cls);
  }

  /** The substituent a tile's conclusion is about, drawn; undefined where there is no single one. */
  answerArt(smiles: string | null): HTMLElement | undefined {
    return smiles ? this.depiction(smiles, STRIP_MOL_W, STRIP_MOL_H) : undefined;
  }

  /** A statement on the pane, with its qualification on the hover. */
  prose(text: string, tip: string): HTMLElement {
    return tipText(text, 'chem-sar-sum-prose', tip);
  }

  /** The same, in the faint style a row uses for its trailing clause. */
  hint(text: string, tip: string): HTMLElement {
    return tipText(text, 'chem-sar-cp-hint', tip);
  }

  /** Why a card has no rows; the card still renders, so finding nothing is visible. */
  reason(text: string): HTMLElement[] {
    return [ui.divText(text, 'chem-sar-cp-hint')];
  }

  /** A chevron heading that shows and hides `body`; the caller owns whatever else sits on the row. */
  foldHead(title: string, body: HTMLElement, open: boolean, cls = '',
    onToggle?: (open: boolean) => void): HTMLElement {
    const icon = ui.iconFA(open ? 'chevron-down' : 'chevron-right');
    icon.classList.add('chem-sar-card-twisty');
    body.style.display = open ? '' : 'none';
    const head = ui.divH([icon, ui.divText(title, 'chem-sar-cp-section-title')],
      `chem-sar-sum-fold ${cls}`.trim());
    head.onclick = () => {
      const show = body.style.display === 'none';
      body.style.display = show ? '' : 'none';
      icon.classList.toggle('fa-chevron-down', show);
      icon.classList.toggle('fa-chevron-right', !show);
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
        // Added after the paint tick, so these depictions need their own.
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

  /** Signed, at two decimals: effects live in tenths, where one decimal makes distinct rows read alike. */
  formatEffect(value: number): string {
    const abs = Math.abs(value);
    return `${value < 0 ? '−' : '+'}${abs < 10 ? abs.toFixed(2) : abs.toFixed(1)}`;
  }

  /** A measured difference: a signed log difference on a log scale, otherwise a fold. A loss is written
   *  as a reciprocal fold, since 10^value rounds every large loss to "0.0×". */
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

  /** Whether the front-runner leads by more than the model error, so it may be named alone. */
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
    if (fit.cvR2 === null)
      return 'Too few compounds to cross-validate the fit: unchecked, not wrong, so nothing is ranked.';

    if (fit.cvR2 < TRUST_R2) {
      return `Cross-validated R² ${fit.cvR2.toFixed(2)}, below ${TRUST_R2}: the components do not predict ` +
        'unseen compounds, so the offsets are not ranked.';
    }
    return null;
  }

  /** Series a core ranking may include: only where a matrix is one core. */
  rankableCores(data: SummaryData): SeriesStat[] {
    return data.coresAreSeries ?
      data.series.filter((stat) => stat.realCells > 0 && stat.cpd >= MIN_SUPPORT) : [];
  }

  /** One tier's cores, best first; mixing tiers would count one compound under two cores. */
  topCores(data: SummaryData): SeriesStat[] {
    const tier = this.tierFilter ?? this.defaultCoreTier(data);
    if (tier === null)
      return [];
    const dir = this.host.higherIsBetter ? 1 : -1;
    return this.rankableCores(data).filter((stat) => stat.tier === tier)
      .sort((a, b) => dir * (b.typical - a.typical) || (a.matrix.id < b.matrix.id ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  /** The tier ranked by default: the one with the most cores. */
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
