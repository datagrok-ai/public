import * as ui from 'datagrok-api/ui';
import {RoleSummary} from './sar-matrix-role-fit';
import {MIN_SUPPORT, RGROUP_MIN_SERIES, SummaryData, SWAP_MIN_PAIRS,
  SWAP_MIN_SERIES} from './sar-matrix-summary-data';
import {count, tipText} from '../sar-matrix-ui-common';
import {FINDING_ROWS, NOTHING_AT_POSITION, orient, Answer, shortSmiles, roleRanks} from './sar-matrix-summary-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {SummaryStartHere} from './sar-matrix-summary-start-here';

export class SummaryOverview {
  private readonly startHere: SummaryStartHere;

  constructor(private readonly kit: SummaryPanel) {
    this.startHere = new SummaryStartHere(kit);
  }

  head(data: SummaryData): HTMLElement {
    return this.kit.scaleRow(data, this.coverage(data));
  }

  keySeries(data: SummaryData): HTMLElement {
    const title = tipText('Key series', 'chem-sar-cp-section-title', `Series fits: ${count(data.fitHolds)} ` +
      `hold, ${count(data.unchecked)} unchecked, ${count(data.lowR2)} weak. Unchecked means unverified, ` +
      'not wrong.');
    return ui.divV([ui.divH([title, this.kit.levelInput()], 'chem-sar-sum-title-row'),
      ui.divV(this.startHere.build(data), 'chem-sar-sum-list')], 'chem-sar-sum-card chem-sar-sum-overview');
  }

  private coverage(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const total = Math.max(1, host.hostRowCount);
    // Kept apart: never assayed calls for a plate, no analog means a library too diverse for any core.
    const unplaced = Math.max(0, host.assayedCount - data.compounds - host.unscalableCount);
    const never = Math.max(0, total - host.assayedCount - data.untested);
    const segs = [
      {n: data.compounds, word: 'measured'},
      {n: data.untested, word: 'held, never assayed'},
      {n: host.unscalableCount, word: 'assayed, cannot be scaled'},
      {n: unplaced, word: 'assayed, no analog'},
      {n: never, word: 'never assayed'},
    ].filter((s) => s.n > 0);
    const label = this.kit.hint(segs.map((s) => `${count(s.n)} ${s.word}`).join(' · '),
      `${count(total)} table rows. Rows, not compounds: the denominator is the table's raw row count, so ` +
      'this ignores any filter on the source grid.' + (host.unscalableCount === 0 ? '' : ` ${host.scalingLabel} ` +
      'needs positive values; set Scaling to "none" to use the values it cannot scale.') +
      (data.axisRole !== null ? '' : ' Cores and substituents were found by cutting the molecules, so ' +
      'substituent labels are local to their series.'));
    return ui.divH([ui.divText('·', 'chem-sar-cp-hint'), label], 'chem-sar-sum-chips chem-sar-sum-coverage');
  }

  /** Each tile carries its own limit, or it reads as a ranking of the chemistry, not of what was made. */
  answers(data: SummaryData): HTMLElement {
    const findings = this.buildFindings(data);
    if (findings !== null)
      return findings;
    const detail = this.detail();
    const effects = this.kit.effects;
    const tiles = ui.divH([
      this.answerTile(`Best ${data.axisRole ?? 'R-group'}`, this.rgroupAnswer(data), detail,
        () => effects.rgroupCardFor(data)),
      this.answerTile(`Best ${data.coreRole ?? 'core'}`, this.coreAnswer(data), detail,
        () => effects.coreCardFor(data)),
      this.answerTile('What improves potency', this.potencyAnswer(data), detail, () => effects.swapCard(data)),
    ], 'chem-sar-sum-answers');
    if (data.axisRole === null)
      return ui.divV([tiles, detail.el]);
    const refusal = ui.divText(this.kit.roleFitRefusal(data) ?? '', 'chem-sar-cp-hint');
    return ui.divV([ui.divH([refusal, this.setupLine(data)], 'chem-sar-sum-title-row'), tiles, detail.el]);
  }

  private buildFindings(data: SummaryData): HTMLElement | null {
    const fit = data.roleFit;
    if (fit === null || this.kit.roleFitRefusal(data) !== null)
      return null;
    const head = ui.divH([ui.div([], 'chem-sar-sum-col'), ui.divText('Best value', 'chem-sar-sum-half-head'),
      ui.divText('Best measured swap', 'chem-sar-sum-half-head')], 'chem-sar-sum-comp-row');
    const body: HTMLElement[] = [
      ui.divH([ui.divText(`What to change · ${count(fit.roles.length)} components`, 'chem-sar-cp-section-title'),
        this.setupLine(data)], 'chem-sar-sum-title-row'),
      head,
      ...fit.roles.slice(0, FINDING_ROWS).map((role) => this.componentRow(data, role)),
    ];
    if (fit.roles.length > FINDING_ROWS) {
      body.push(...this.kit.lazyFold(`${fit.roles.length - FINDING_ROWS} more`, () => fit.roles.slice(FINDING_ROWS)
        .map((role) => this.componentRow(data, role))));
    }
    return ui.divV(body, 'chem-sar-sum-comp chem-sar-sum-card');
  }

  /** One shared slot; a second click on the same trigger closes it. */
  private detail(): {el: HTMLElement, toggle: (trigger: HTMLElement, build: () => HTMLElement) => void} {
    const el = ui.div([], 'chem-sar-sum-detail');
    let open: HTMLElement | null = null;
    const toggle = (trigger: HTMLElement, build: () => HTMLElement): void => {
      open?.classList.remove('chem-sar-sum-open');
      ui.empty(el);
      open = open === trigger ? null : trigger;
      if (open === null)
        return;
      open.classList.add('chem-sar-sum-open');
      el.appendChild(build());
      this.kit.flushPaints();
    };
    return {el, toggle};
  }

  private setupLine(data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const folded = host.roleColumns.filter((name) => name !== data.axisRole);
    const parts = [`core: ${data.coreRole}`, `columns: ${data.axisRole}`];
    if (folded.length > 0)
      parts.push(`folded into rows: ${folded.join(', ')}`);
    return this.kit.hint(parts.join(' · '),
      'Which component runs across the matrix columns only changes the grid\'s layout; every component ' +
      'is ranked and has its swaps pooled either way.');
  }

  private componentRow(data: SummaryData, role: RoleSummary): HTMLElement {
    const host = this.kit.host;
    const detail = this.detail();
    const top = roleRanks(role) ? role.levels[0] : undefined;
    const best = ui.divH(top === undefined ? [ui.divText('not ranked', 'chem-sar-cp-hint')] : [
      top.value === '' ? ui.divText(NOTHING_AT_POSITION, 'chem-sar-cp-hint') :
        this.kit.depiction(top.value),
      this.figure(this.kit.formatEffect(top.coef), `vs average · ${count(top.n)} cpd`),
    ], 'chem-sar-sum-half');
    ui.tooltip.bind(best, () => `${role.name} typically moves ${host.activityColumnName} by ` +
      `${role.spread.toFixed(2)} (a standard deviation, not a range). Click for its full ranking.`);
    best.onclick = () => detail.toggle(best, () => this.kit.effects.roleCard(data, role));

    const pool = data.swapsByRole.get(role.name)?.[0];
    const only = ` Pairs alike in every component but ${role.name}, pooled across the rest of the table.`;
    let tip: string;
    let parts: HTMLElement[];
    if (pool !== undefined) {
      const {from, to, worst, widest, mean} = orient(pool);
      parts = [this.kit.pairArt(from, to),
        this.figure(worst > 0 ? `≥ ${this.kit.formatDelta(worst)}` : `mean ${this.kit.formatDelta(mean)}`,
          `${count(pool.n)} pairs`)];
      tip = `${this.kit.formatDelta(worst)} to ${this.kit.formatDelta(widest)}. ` + (worst > 0 ?
        'Improved potency in every measured pair.' : 'No swap improved every pair; this one has the best floor.') +
        ' Measured, not fitted.' + only + ' Click for the full pool.';
    } else {
      parts = [ui.divText(`nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} contexts`,
        'chem-sar-cp-hint')];
      tip = 'A swap needs the same substitution measured on several backgrounds.' + only + ' Click for details.';
    }
    const swap = ui.divH(parts, 'chem-sar-sum-half chem-sar-sum-swap');
    ui.tooltip.bind(swap, () => tip);
    swap.onclick = () => detail.toggle(swap, () => this.kit.effects.swapCard(data, role.name));

    const name = ui.divV([ui.divText(role.name, 'chem-sar-card-title'),
      ui.divText(`spread ${role.spread.toFixed(2)}`, 'chem-sar-card-desc')], 'chem-sar-sum-col');
    ui.tooltip.bind(name, () => `Spread: how far ${host.activityColumnName} typically moves with the ${role.name} ` +
      'value, noise removed. A standard deviation, not a range. Click a value or a swap for the evidence.');
    return ui.divV([ui.divH([name, best, swap], 'chem-sar-sum-comp-row'), detail.el]);
  }

  private figure(value: string, caption: string): HTMLElement {
    return ui.divV([ui.divText(value, 'chem-sar-card-title'), ui.divText(caption, 'chem-sar-cp-hint')],
      'chem-sar-sum-figure');
  }

  private answerTile(title: string, text: Answer, detail: ReturnType<SummaryOverview['detail']>,
    build: () => HTMLElement): HTMLElement {
    const go = ui.iconFA('chevron-down');
    go.classList.add('chem-sar-card-twisty');
    const body: HTMLElement[] = [
      ui.divH([ui.divText(title, 'chem-sar-cp-section-title'), go], 'chem-sar-sum-answer-head'),
    ];
    if (text.art !== undefined)
      body.push(text.art);
    if (text.answer !== undefined)
      body.push(ui.divText(text.answer, 'chem-sar-sum-prose'));
    // Not `chem-sar-card-desc`: this is a sentence, and that class clips to one line.
    body.push(ui.divText(text.negative, 'chem-sar-cp-hint'));
    const tile = ui.divV(body, 'chem-sar-sum-answer');
    ui.tooltip.bind(tile, () => `${text.answer === undefined ? '' : `${text.answer} — `}${text.negative}. Click for ` +
      'the evidence behind it.');
    tile.onclick = () => detail.toggle(tile, build);
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
      // Name the series floor that the empty list actually failed.
      return {answer: `None in ${data.rgroupsThin.length ? RGROUP_MIN_SERIES : SWAP_MIN_SERIES} series`,
        negative: 'no substituent comes first at its position in that many series whose additive fit holds'};
    }
    const top = rows[0];
    const partner = data.axisRole === null ? 'series' : 'core';
    if (!this.kit.leads(top.magnitude, rows[1]?.magnitude ?? null, data.modelError)) {
      return {answer: 'Two lead equally',
        negative: `the ${partner} decides — ${top.k} of ${top.tried} against ${rows[1].k} of ` +
          `${rows[1].tried}`,
        art: this.kit.pairArt(top.subst, rows[1].subst, 'or')};
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
    const dir = this.kit.dir;
    const top = cores[0];
    if (!this.kit.leads(dir * top.typical, cores.length > 1 ? dir * cores[1].typical : null, data.modelError)) {
      return {answer: 'No core runs clear',
        negative: `${depth}the top ${cores.length} sit inside ` +
          `${host.formatActivity(data.modelError ?? 0)} model error of each other`};
    }
    return {answer: `${top.matrix.label} · typical ${host.formatActivity(top.typical)} over ` +
      `${count(top.cpd)} cpd`, negative,
    art: this.kit.answerArt(top.matrix.rows[0].keySmiles)};
  }

  private potencyAnswer(data: SummaryData): Answer {
    const negative = 'both compounds measured, nothing fitted — one R-group swapped, the rest identical';
    const top = data.swaps[0];
    if (top !== undefined) {
      const {from, to, worst, mean} = orient(top);
      // A negative floor means the swap lost at least once, so it is not stated as an improvement.
      return {art: this.kit.pairArt(from, to),
        negative: worst > 0 ? negative : 'no swap improves potency in every pair it was measured in — this ' +
          `one is the best floor; mean ${this.kit.formatDelta(mean)}`};
    }
    return {answer: data.swapCandidates === 0 ? 'No row carries two measured substituents' :
      `Nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} series`, negative};
  }
}
