/* The Summary's Overview segment, the one the analysis opens on: the answers, the component band
   and its measured swaps, and Start here. */
import * as ui from 'datagrok-api/ui';
import {RoleSummary} from './sar-matrix-role-fit';
import {MIN_SUPPORT, RGROUP_MIN_SERIES, SummaryData, SWAP_MIN_PAIRS,
  SWAP_MIN_SERIES} from './sar-matrix-summary-data';
import {count, STRIP_MOL_H, STRIP_MOL_W, tipText} from '../sar-matrix-ui-common';
import {FINDING_ROWS, PANE_EFFECTS, orient, Answer, shortSmiles, swapAnchor,
  roleRanks} from './sar-matrix-summary-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {SummaryStartHere} from './sar-matrix-summary-start-here';

export class SummaryOverview {
  private readonly startHere: SummaryStartHere;

  constructor(private readonly kit: SummaryPanel) {
    this.startHere = new SummaryStartHere(kit);
  }

  build(data: SummaryData): HTMLElement {
    return ui.divV([
      this.coverage(data),
      this.buildAnswers(data),
      ui.divText('Start here', 'chem-sar-cp-section-title'),
      ui.divV(this.startHere.build(data), 'chem-sar-sum-list'),
      this.buildTotals(data),
    ], 'chem-sar-sum-overview');
  }

  /** The verdict on the fits, and a warning where values could not be scaled. */
  private buildTotals(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const trust = tipText(`${count(data.fitHolds)} fits hold · ${count(data.unchecked)} unchecked · ` +
      `${count(data.lowR2.length)} do not`, 'chem-sar-cp-hint', 'Per series, not per compound: one compound ' +
      'sits in several series. Unchecked means unverified, not wrong.');

    const parts = [trust];
    if (host.unscalableCount > 0) {
      parts.push(tipText(`${count(host.unscalableCount)} values cannot be scaled by ${host.scalingLabel}`,
        'chem-sar-cp-hint chem-sar-chip-partial', `${host.scalingLabel} needs positive values, so these are ` +
        'in no matrix. Set Scaling to "none" to use them.'));
    }
    return ui.divH(parts, 'chem-sar-sum-chips');
  }

  /** What fraction of the table the analysis sees: a bar over a label naming its segments in the same
   *  order, so the label is the legend. */
  private coverage(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const total = Math.max(1, host.hostRowCount);
    // Kept apart: never assayed calls for a plate, no analog means a library too diverse for any core.
    const unplaced = Math.max(0, host.assayedCount - data.compounds - host.unscalableCount);
    const never = Math.max(0, total - host.assayedCount - data.untested);
    const segs = [
      {n: data.compounds, word: 'measured', kind: 'real'},
      {n: data.untested, word: 'held, never assayed', kind: 'held'},
      {n: host.unscalableCount, word: 'assayed, cannot be scaled', kind: 'bad'},
      {n: unplaced, word: 'assayed, no analog to pair with', kind: 'out'},
      {n: never, word: 'never assayed', kind: 'none'},
    ].filter((s) => s.n > 0);
    const bar = ui.divH(segs.map((s) => {
      const seg = ui.div([], `chem-sar-sum-cov-seg chem-sar-sum-cov-${s.kind}`);
      seg.style.flexGrow = `${s.n}`;
      return seg;
    }), 'chem-sar-sum-cov');
    const label = this.kit.hint(segs.map((s) => `${count(s.n)} ${s.word}`).join(' · '),
      `${count(total)} table rows. Rows, not compounds: the denominator is the table's raw row count, so ` +
      'this ignores any filter on the source grid.');
    return ui.divV([bar, ui.divH([label, ui.divText('·', 'chem-sar-cp-hint'), this.setupLine(data)],
      'chem-sar-sum-chips')], 'chem-sar-sum-coverage');
  }

  /** Which R-group, which core and which swap buys the most. Each tile carries its own limit, or it reads
   *  as a ranking of the chemistry rather than of what was made. */
  private buildAnswers(data: SummaryData): HTMLElement {
    const findings = this.buildFindings(data);
    if (findings !== null)
      return findings;
    return ui.divH([
      this.answerTile(`Best ${data.axisRole ?? 'R-group'}`, this.rgroupAnswer(data), 'rgroup',
        () => this.kit.expandTopRGroup = true),
      this.answerTile(`Best ${data.coreRole ?? 'core'}`, this.coreAnswer(data), 'cores'),
      this.answerTile('What improves potency', this.potencyAnswer(data), 'swaps'),
    ], 'chem-sar-sum-answers');
  }

  /** One row per component column and one per best measured swap, full width so structures can be
   *  drawn. Null when there is no component fit; the three tiles show instead. */
  private buildFindings(data: SummaryData): HTMLElement | null {
    const fit = data.roleFit;
    if (fit === null || this.kit.roleFitRefusal(data) !== null)
      return null;
    const rows = fit.roles.map((role, index) => this.componentRow(role, index));
    const body: HTMLElement[] = [
      ui.divText(`spans = the range of ${this.kit.host.activityColumnName} across that component's ` +
        `values · offsets are against the library mean of ${this.kit.host.formatActivity(fit.mean)}`,
      'chem-sar-cp-hint'),
      ...rows.slice(0, FINDING_ROWS),
    ];
    // Fold the components that move the endpoint least, so the swaps stay on the first screen.
    if (rows.length > FINDING_ROWS)
      body.push(this.moreComponents(rows.slice(FINDING_ROWS)));

    const swapRows = fit.roles.map((role) => this.swapRow(data, role.name));
    const withPairs = fit.roles.filter((role) => (data.swapsByRole.get(role.name) ?? []).length > 0);
    // A folded half still shows its headline.
    const lead = fit.roles[0];
    return ui.divV([
      this.bandGroup('comp', `What to change · ${count(fit.roles.length)} components`,
        `${lead.name} leads, spans ${lead.spread.toFixed(2)}`, ui.divV(body)),
      this.bandGroup('swap', `Best measured swap · ${count(fit.roles.length)} components`,
        withPairs.length === 0 ? 'none pooled' :
          `${withPairs[0].name} ${this.swapHeadline(data, withPairs[0].name)}`, ui.divV(swapRows)),
    ], 'chem-sar-sum-comp chem-sar-sum-card');
  }

  /** Which column is the core, which one the matrix columns enumerate, and which fold into the row. */
  private setupLine(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    if (data.axisRole === null) {
      return this.kit.hint('Cores and substituents found by cutting the molecules.',
        'Substituent labels from cutting are local to their series and mean nothing across series.');
    }
    const folded = host.roleColumns.filter((name) => name !== data.axisRole);
    const parts = [`core: ${data.coreRole}`, `across the matrix columns: ${data.axisRole}`];
    if (folded.length > 0)
      parts.push(`folded into the row: ${folded.join(', ')}`);
    return this.kit.hint(parts.join(' · '),
      'Which component runs across the matrix columns only changes the grid\'s layout; every component ' +
      'is ranked and has its swaps pooled either way.');
  }

  /** One foldable half of the landing band, whose heading states its answer when shut. */
  private bandGroup(key: string, title: string, headline: string, body: HTMLElement): HTMLElement {
    const open = !this.kit.foldedBands.has(key);
    const lead = ui.divText(headline, 'chem-sar-cp-hint');
    lead.style.display = open ? 'none' : '';
    const head = this.kit.foldHead(title, body, open, '', (show) => {
      lead.style.display = show ? 'none' : '';
      if (show)
        this.kit.foldedBands.delete(key);
      else
        this.kit.foldedBands.add(key);
    });
    head.appendChild(lead);
    return ui.divV([head, body], `chem-sar-sum-band-${key}`);
  }

  /** The swap's conclusion in one phrase, for the group heading when it is shut. */
  private swapHeadline(data: SummaryData, role: string): string {
    const top = data.swapsByRole.get(role)![0];
    const {worst, widest} = orient(top);
    return worst > 0 ? `≥ ${this.kit.formatDelta(worst)} in ${top.n} pairs` :
      `${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)}`;
  }

  /** One component: how much it moves the endpoint, and its best value. */
  private componentRow(role: RoleSummary, index: number): HTMLElement {
    const top = roleRanks(role) ? role.levels[0] : undefined;
    const parts: HTMLElement[] = [
      ui.divText(role.name, 'chem-sar-card-title chem-sar-sum-col'),
      ui.divText(`spans ${role.spread.toFixed(2)}`, 'chem-sar-card-desc chem-sar-sum-col'),
    ];
    if (top === undefined)
      parts.push(ui.divH([ui.divText('not ranked', 'chem-sar-cp-hint')], 'chem-sar-sum-slot'));
    else {
      parts.push(ui.divH([top.value === '' ? this.kit.hint('nothing at this position',
        `Compounds with ${role.name} left blank score best: hydrogen, or nothing recorded.`) :
        this.kit.depiction(top.value, STRIP_MOL_W, STRIP_MOL_H)], 'chem-sar-sum-slot'));
      parts.push(ui.divText(this.kit.formatEffect(top.coef), 'chem-sar-card-title chem-sar-sum-col'));
      parts.push(ui.divText(`over ${count(top.n)} compounds`, 'chem-sar-cp-hint'));
    }
    const row = ui.divH(parts, 'chem-sar-card');
    ui.tooltip.bind(row, () => `${role.name} moves ${this.kit.host.activityColumnName} by ` +
      `${role.spread.toFixed(2)} across its values. Click for its full ranking.`);
    row.onclick = () => this.kit.showPane(PANE_EFFECTS, `role-${index}`);
    return row;
  }

  /** The components that move the endpoint least, folded away behind their own count. */
  private moreComponents(rows: HTMLElement[]): HTMLElement {
    const body = ui.divV(rows);
    const head = this.kit.foldHead(`${rows.length} more`, body, false);
    return ui.divV([head, body]);
  }

  /** The best measured swap of one component, laid out like the component rows above it. */
  private swapRow(data: SummaryData, role: string): HTMLElement {
    const parts: HTMLElement[] = [
      ui.divText(role, 'chem-sar-card-title chem-sar-sum-col'),
      ui.divText('swapped', 'chem-sar-card-desc chem-sar-sum-col'),
    ];
    const top = data.swapsByRole.get(role)?.[0];
    const only = ` Pairs alike in every component but ${role}, pooled across the rest of the table.`;
    let tip: string;
    if (top !== undefined) {
      const {from, to, worst, widest, mean} = orient(top);
      parts.push(this.pairArt(from, to, 'chem-sar-sum-slot'));
      parts.push(ui.divText(worst > 0 ? `≥ ${this.kit.formatDelta(worst)}` :
        `mean ${this.kit.formatDelta(mean)}`, 'chem-sar-card-title chem-sar-sum-col'));
      parts.push(ui.divText(`${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)} · ` +
        `${count(top.n)} pairs · ${count(top.roots.size)} contexts`, 'chem-sar-cp-hint'));
      tip = (worst > 0 ? 'Improved potency in every measured pair.' : 'No swap improved every pair; this ' +
        'one has the best floor.') + ' Measured, not fitted.' + only + ' Click for the full pool.';
    } else {
      parts.push(ui.divH([ui.divText(`nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} contexts`,
        'chem-sar-cp-hint')], 'chem-sar-sum-slot'));
      tip = 'A swap needs the same substitution measured on several backgrounds.' + only + ' Click for details.';
    }
    const row = ui.divH(parts, 'chem-sar-card');
    ui.tooltip.bind(row, () => tip);
    row.onclick = () => this.kit.showPane(PANE_EFFECTS, swapAnchor(role));
    return row;
  }

  /** One conclusion: the name, the number and one clause; the evidence is on Effects. */
  private answerTile(title: string, text: Answer, anchor: string, onOpen?: () => void): HTMLElement {
    const go = ui.iconFA('chevron-right');
    go.classList.add('chem-sar-card-twisty');
    const body: HTMLElement[] = [
      ui.divH([ui.divText(title, 'chem-sar-cp-section-title'), go], 'chem-sar-sum-answer-head'),
    ];
    if (text.art !== undefined)
      body.push(text.art);
    body.push(ui.divText(text.answer, 'chem-sar-sum-prose'));
    // Not `chem-sar-card-desc`: this is a sentence, and that class clips to one line.
    body.push(ui.divText(text.negative, 'chem-sar-cp-hint'));
    const tile = ui.divV(body, 'chem-sar-sum-answer');
    // The value line clips, so the full answer is on the hover.
    ui.tooltip.bind(tile, () => `${text.answer} — ${text.negative}. The evidence behind this answer is ` +
      'on the Effects segment; click to go to it.');
    tile.onclick = () => {
      onOpen?.();
      this.kit.showPane(PANE_EFFECTS, anchor);
    };
    return tile;
  }

  private rgroupAnswer(data: SummaryData): Answer {
    if (!this.kit.host.activityIsLog) {
      return {answer: 'Not answerable on a raw scale',
        negative: 'set Scaling to lg or −lg, or mark the column higher-is-better if it is a pIC50'};
    }
    const loser = data.rgroupLosers[0];
    const negative = loser === undefined ?
      'over each series\' own most-common substituent, not over your whole library' :
      `stop making ${shortSmiles(loser.subst)} — last in ${loser.k} of ${loser.tried} tried`;
    const rows = data.rgroups;
    if (rows.length === 0) {
      return {answer: `None in ${RGROUP_MIN_SERIES} series`,
        negative: 'no substituent comes first at its position in that many series whose additive fit holds'};
    }
    const top = rows[0];
    const partner = data.axisRole === null ? 'series' : 'core';
    if (!this.kit.leads(top.magnitude, rows[1]?.magnitude ?? null, data.modelError)) {
      return {answer: `Two lead equally`,
        negative: `the ${partner} decides — ${top.k} of ${top.tried} against ${rows[1].k} of ` +
          `${rows[1].tried}`,
        art: this.pairArt(top.subst, rows[1].subst, 'chem-sar-sum-pair', 'or')};
    }
    return {answer: (top.magnitude === null ? '' : `${this.kit.formatEffect(top.magnitude)} · `) +
      `first on ${top.k} of ${top.tried}`, negative, art: this.kit.answerArt(top.subst)};
  }

  private coreAnswer(data: SummaryData): Answer {
    const host = this.kit.host;
    if (!data.coresAreSeries) {
      return {answer: 'Not comparable across series',
        negative: 'a core recurs only inside its own lineage — nothing to pool across series'};
    }
    const cores = this.kit.topCores(data);
    if (cores.length === 0)
      return {answer: `None with ${MIN_SUPPORT} compounds`, negative: 'mean of what was made'};
    // Only cores cut at the same depth compare, so the depth is named.
    const depth = this.kit.coreTiers(data).length > 1 ? `among the L${cores[0].tier} cores; ` : '';
    const negative = `${depth}mean of what was made — a core only ever paired with good caps scores high`;
    const dir = host.higherIsBetter ? 1 : -1;
    const top = cores[0];
    if (!this.kit.leads(dir * top.typical, cores.length > 1 ? dir * cores[1].typical : null, data.modelError)) {
      return {answer: `No core runs clear`,
        negative: `${depth}the top ${cores.length} sit inside ` +
          `${host.formatActivity(data.modelError ?? 0)} model error of each other`};
    }
    return {answer: `${top.matrix.label} · typical ${host.formatActivity(top.typical)} over ` +
      `${count(top.cpd)} cpd`, negative,
    art: this.kit.answerArt(top.matrix.rows[0]?.keySmiles ?? null)};
  }

  private potencyAnswer(data: SummaryData): Answer {
    const negative = 'both compounds measured, nothing fitted — one R-group swapped, the rest identical';
    const top = data.swaps[0];
    if (top !== undefined) {
      const {from, to, worst, widest, mean} = orient(top);
      const art = this.pairArt(from, to, 'chem-sar-sum-pair');
      if (worst > 0)
        return {answer: `≥ ${this.kit.formatDelta(worst)} in ${top.n} pairs`, negative, art};
      // A negative floor means the swap lost at least once, so it is not stated as an improvement.
      return {answer: `${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)} over ${top.n} pairs`,
        negative: `no swap improves potency in every pair it was measured in — this one is the best ` +
          `floor; mean ${this.kit.formatDelta(mean)}`, art};
    }
    // Nothing pooled: what the gate rejected is the finding, and there is nothing to draw.
    const gate = data.swapCandidates === 0 ? 'No row carries two measured substituents' :
      `Nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} series`;
    const best = data.swapBest;
    if (best === null || best.best === null)
      return {answer: gate, negative};
    return {answer: gate, negative: `the largest single measured move is ` +
      `${this.kit.formatDelta(Math.abs(best.best.delta))} in ${best.best.matrix.label}`};
  }

  /** Two structures and the glyph between them: the shape every swap is drawn in. */
  private pairArt(from: string, to: string, cls: string, glyph = '→'): HTMLElement {
    return ui.divH([this.kit.depiction(from, STRIP_MOL_W, STRIP_MOL_H),
      ui.divText(glyph, 'chem-sar-sum-arrow'),
      this.kit.depiction(to, STRIP_MOL_W, STRIP_MOL_H)], cls);
  }
}
