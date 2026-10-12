import * as ui from 'datagrok-api/ui';
import {SarMatrix} from '../sar-matrix-types';
import {MIN_SUPPORT, SeriesStat, SummaryData, SummaryHost, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {chipBadge, CORE_BG_ARGB, isStructure, MatrixCellRef, paintMoleculeOnColor, SUM_MOL_H,
  SUM_MOL_W, tipText} from '../sar-matrix-ui-common';

export class SummaryKit {
  pendingPaints: (() => void)[] = [];
  /** Tier the summary is read at; null for all. */
  protected tierFilter: number | null = null;

  constructor(readonly host: SummaryHost) {}

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
    /** Makes a click on the row expand it, returning whether it is now open; without it the row is still and
     *  only the open button acts. */
    onChevron?: () => boolean,
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
    const open = ui.iconFA('th', (e: MouseEvent) => {
      e.stopPropagation();
      opts.onClick();
    }, toggle === undefined ? opts.tip : 'Open in the SAR Matrix');
    open.classList.add('chem-sar-struct-icon');
    parts.push(open);
    if (toggle === undefined)
      return ui.divH(parts, 'chem-sar-card chem-sar-sum-static');
    const [twisty, point] = this.twisty();
    parts.push(twisty);
    const row = ui.divH(parts, 'chem-sar-card');
    ui.tooltip.bind(row, () => opts.tip);
    row.onclick = () => point(toggle());
    return row;
  }

  /** A chevron that points down while its section is open. */
  private twisty(): [HTMLElement, (open: boolean) => void] {
    const icon = ui.iconFA('chevron-right');
    icon.classList.add('chem-sar-card-twisty');
    return [icon, (open) => {
      icon.classList.toggle('fa-chevron-down', open);
      icon.classList.toggle('fa-chevron-right', !open);
    }];
  }

  /** Canvas now, RDKit on the next tick, so a tab switch is not held up by drawing. */
  depiction(smiles: string | null): HTMLElement {
    const canvas = ui.canvas(SUM_MOL_W, SUM_MOL_H);
    canvas.classList.add('chem-sar-card-core');
    if (smiles) {
      this.pendingPaints.push(() => {
        if (isStructure(smiles)) {
          paintMoleculeOnColor(canvas, smiles, SUM_MOL_W, SUM_MOL_H, CORE_BG_ARGB);
          return;
        }
        // A component given by name, such as "VHL", has nothing to draw and is written out instead.
        const name = ui.divText(smiles, 'chem-sar-cp-frag-name');
        name.style.width = `${SUM_MOL_W}px`;
        name.style.height = `${SUM_MOL_H}px`;
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

  /** Only a fit that needs a caveat gets a badge; no confidence means unchecked, not failing. */
  trustBadges(matrix: SarMatrix): HTMLElement[] {
    const conf = matrix.confidence;
    if (!conf) {
      return [this.badge('fit unchecked', 'Too few measured cells to cross-validate this series\' fit, or all ' +
        'values identical: its predictions are unverified, not wrong.')];
    }
    if (conf.r2 >= TRUST_R2)
      return [];
    return [this.badge('fit weak', `The additive fit does not predict held-out compounds well (R² ` +
      `${conf.r2.toFixed(2)}, below ${TRUST_R2}; ± ${this.host.formatActivity(conf.rmse)} over ${conf.n} of ` +
      `${conf.total} cells): treat this series' predictions with caution.`, true)];
  }

  pairArt(from: string, to: string, glyph = '→'): HTMLElement {
    return ui.divH([this.depiction(from), ui.divText(glyph, 'chem-sar-sum-arrow'), this.depiction(to)],
      'chem-sar-sum-pair');
  }

  answerArt(smiles: string | null): HTMLElement | undefined {
    return smiles ? this.depiction(smiles) : undefined;
  }

  prose(text: string, tip: string): HTMLElement {
    return tipText(text, 'chem-sar-sum-prose', tip);
  }

  hint(text: string, tip: string): HTMLElement {
    return tipText(text, 'chem-sar-cp-hint', tip);
  }

  reason(text: string): HTMLElement[] {
    return [ui.divText(text, 'chem-sar-cp-hint')];
  }

  /** A closed fold whose body is built on its first opening. */
  lazyFold(title: string, fill: () => HTMLElement[]): [HTMLElement, HTMLElement] {
    const lazy = this.lazyBody(fill);
    const [icon, point] = this.twisty();
    const head = ui.divH([icon, ui.divText(title, 'chem-sar-cp-section-title')], 'chem-sar-sum-fold');
    head.onclick = () => point(lazy.toggle());
    return [head, lazy.body];
  }

  /** A body filled on its first opening, so an expand nobody opens draws nothing. */
  lazyBody(fill: () => HTMLElement[]): {body: HTMLElement, toggle: () => boolean} {
    const body = ui.divV([]);
    body.style.display = 'none';
    let filled = false;
    const toggle = (): boolean => {
      if (body.style.display !== 'none') {
        body.style.display = 'none';
        return false;
      }
      if (!filled) {
        filled = true;
        body.append(...fill());
        // Added after the paint tick, so these depictions need their own.
        this.flushPaints();
      }
      body.style.display = '';
      return true;
    };
    return {body, toggle};
  }

  /** Which way a better value points: 1 when higher is better. */
  get dir(): number {
    return this.host.higherIsBetter ? 1 : -1;
  }

  /** Signed, at two decimals: effects live in tenths, where one decimal makes distinct rows read alike. */
  formatEffect(value: number): string {
    const abs = Math.abs(value);
    return `${value < 0 ? '−' : '+'}${abs < 10 ? abs.toFixed(2) : abs.toFixed(1)}`;
  }

  /** A loss is a reciprocal fold, since 10^value rounds every large loss to "0.0×". */
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

  openMatrixTip(matrix: SarMatrix): string {
    return `Open ${matrix.label} in the SAR Matrix`;
  }

  /** A front-runner is named alone only when it clears the model error. */
  leads(first: number | null, second: number | null, modelError: number | null): boolean {
    if (second === null)
      return true;
    if (first === null)
      return false;
    return first - second > (modelError ?? 0);
  }

  roleFitRefusal(data: SummaryData): string | null {
    const fit = data.roleFit;
    if (fit === null) {
      return `No component column carries two values with ${MIN_SUPPORT} or more compounds each, so ` +
        'there is nothing to compare.';
    }
    if (!fit.converged) {
      return `The fit did not settle within ${fit.sweeps} passes. Each component's effects are ` +
        'still shrunk toward zero by a different amount and cannot be compared, so nothing is ranked.';
    }
    if (fit.cvR2 === null)
      return 'Too few compounds to cross-validate the fit: unchecked, not wrong, so nothing is ranked.';

    if (fit.cvR2 < TRUST_R2) {
      return `Cross-validated R² ${fit.cvR2.toFixed(2)}, below ${TRUST_R2}: the components do not predict ` +
        'unseen compounds, so the effects are not ranked.';
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
    return this.rankableCores(data).filter((stat) => stat.tier === tier)
      .sort((a, b) => this.dir * (b.typical - a.typical) || (a.matrix.id < b.matrix.id ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  private defaultCoreTier(data: SummaryData): number | null {
    const tiers = this.coreTiers(data);
    if (tiers.length === 0)
      return null;
    const counts = new Map(tiers.map((tier) =>
      [tier, this.rankableCores(data).filter((stat) => stat.tier === tier).length]));
    return tiers.reduce((best, tier) => counts.get(tier)! > counts.get(best)! ? tier : best, tiers[0]);
  }

  coreTiers(data: SummaryData): number[] {
    return [...new Set(this.rankableCores(data).map((stat) => stat.tier))].sort((a, b) => a - b);
  }
}
