import * as ui from 'datagrok-api/ui';
import {MIN_SUPPORT, SeriesStat, SummaryData, SWAP_MIN_PAIRS, SWAP_MIN_SERIES,
  SWAP_ROW_CAP} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import type {SummaryPanel} from './sar-matrix-summary-panel';
import {orient, roleRanks} from './sar-matrix-summary-common';
import type {RoleSummary} from './sar-matrix-role-fit';
import {SummaryRoleCards} from './sar-matrix-summary-role-cards';
import {SummaryRGroupCard} from './sar-matrix-summary-rgroup-card';

export class SummaryEffects {
  private readonly roleCards: SummaryRoleCards;
  private readonly rgroupCard: SummaryRGroupCard;

  constructor(private readonly kit: SummaryPanel) {
    this.roleCards = new SummaryRoleCards(kit, this);
    this.rgroupCard = new SummaryRGroupCard(kit, this);
  }

  rgroupCardFor(data: SummaryData): HTMLElement {
    return this.rgroupCard.build(data);
  }

  coreCardFor(data: SummaryData): HTMLElement {
    return data.coresAreSeries ? this.buildCoreCard(data, this.kit.topCores(data)) :
      ui.divText('Cores do not compare across series here; each series\' best core is in its Key series row.',
        'chem-sar-cp-hint');
  }

  roleCard(data: SummaryData, role: RoleSummary): HTMLElement {
    return this.roleCards.roleCard(data, role);
  }

  /** Whether the component fit ranks this component, so its own card covers it. */
  roleRanked(data: SummaryData, name: string | null): boolean {
    const role = data.roleFit?.roles.find((other) => other.name === name);
    return role !== undefined && this.kit.roleFitRefusal(data) === null && roleRanks(role);
  }

  card(title: string, subtitle: string, body: HTMLElement[], footer?: HTMLElement): HTMLElement {
    const parts = [ui.divText(title, 'chem-sar-cp-section-title'),
      ui.divText(subtitle, 'chem-sar-cp-hint'), ...body];
    if (footer !== undefined)
      parts.push(footer);
    return ui.divV(parts, 'chem-sar-sum-card');
  }

  swapCard(data: SummaryData, role?: string): HTMLElement {
    const host = this.kit.host;
    const pools = role === undefined ? data.swaps : data.swapsByRole.get(role) ?? [];
    const title = role === undefined ? 'What buys potency inside these series' :
      `Swapping ${role} — what it bought`;
    const subtitle = 'Measured pairs alike in everything but one R-group.';
    // With component columns every pair has its own background, so a count of those only repeats the pairs.
    const across = (pool: {roots: Set<string>}): string =>
      role === undefined ? ` · ${pool.roots.size} series` : '';
    const show = (value: number): string => this.kit.formatDelta(value);
    const rows = pools.map((pool) => {
      const {forward, from, to, worst, widest: bestOf, mean, up} = orient(pool);
      const allUp = up === pool.n;
      // Null for a pool grouped off the component columns: no single cell to open.
      const best = pool.best;
      const ci = best === null ? -1 : forward ? best.ciTo : best.ciFrom;
      const position = best === null ? '' : best.matrix.columns[ci].position;
      // Badges, not descriptor text, which clips to one line.
      const badges: HTMLElement[] = [];
      if (allUp) {
        badges.push(this.kit.badge('every pair improved', 'Every measured pair improved, none got worse. With only ' +
          '3 pairs that still happens one time in four by chance.'));
      }
      if (pool.sampled) {
        badges.push(this.kit.badge('sampled row', 'One contributing row carried more than ' +
          `${SWAP_ROW_CAP} measured cells; its most and least potent halves were kept and mid-range ` +
          'pairs dropped.', true));
      }
      // Only on a log scale, where the delta and the model error share units.
      if (host.activityIsLog && data.modelError !== null && Math.abs(worst) < data.modelError) {
        badges.push(this.kit.badge('under resolution', 'The weakest pair is smaller than this analysis\' own ' +
          'typical prediction error, so the whole row may be noise.', true));
      }
      return this.kit.summaryRow({
        depiction: this.kit.pairArt(from, to),
        name: best === null ? `${role} swap` : `${best.matrix.label} · ${position}`,
        badges,
        desc: `mean ${show(mean)} (${show(worst)} → ${show(bestOf)}) · ${pool.n} pairs${across(pool)}`,
        value: show(worst),
        caption: 'weakest pair',
        valueTip: `The smallest change among its ${pool.n} measured pairs${across(pool)}: every pair moved at ` +
          `least ${show(worst)}. Measured, not fitted or predicted.` +
          (pool.sampled ? ` A group over ${SWAP_ROW_CAP} compounds kept only its extremes.` : ''),
        tip: best === null ? `Select every compound whose ${role} is the one on the right` :
          this.kit.openTip(best.matrix, best.ri, ci),
        onClick: best === null ? () => host.selectRoleValue(data.roleColumn.get(role!) ?? role!, to) :
          () => host.revealCell(best.matrix, best.ri, ci, position),
      });
    });
    if (rows.length > 0)
      return this.card(title, subtitle, rows);
    return this.card(title, subtitle, this.kit.reason(data.swapCandidates === 0 && role === undefined ?
      'No row of any series carries two measured substituents, so there is no swap to measure.' :
      `No ${role === undefined ? '' : `${role} `}swap clears ${SWAP_MIN_PAIRS} measured pairs in ` +
      `${SWAP_MIN_SERIES} ${role === undefined ? 'series' : 'contexts'}.`));
  }

  private buildCoreCard(data: SummaryData, cores: SeriesStat[]): HTMLElement {
    const host = this.kit.host;
    const title = `${data.coreRole ?? 'Cores'} — the average of what was made on each`;
    const subtitle = 'Mean activity of each core\'s measured compounds, not adjusted for its substituents: a ' +
      'core only tried with good groups scores high. Unlike the fitted effects, these compare between cores.';
    if (cores.length === 0)
      return this.card(title, subtitle, this.kit.reason(`No core holds ${MIN_SUPPORT} measured compounds.`));
    const body = cores.map((stat) => this.kit.summaryRow({
      depiction: this.kit.depiction(stat.matrix.rows[0].coreSmiles),
      name: stat.matrix.label,
      badges: this.kit.trustBadges(stat.matrix),
      desc: `${count(stat.cpd)} cpd · best ${host.formatActivity(stat.best!.value)} · ` +
        `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}`,
      value: host.formatActivity(stat.typical),
      caption: 'typical',
      valueTip: 'Count-weighted fitted mean of what was made on this core; good caps lift it.',
      tip: this.kit.openMatrixTip(stat.matrix),
      onClick: () => host.revealMatrix(stat.matrix),
    }));
    const dir = this.kit.dir;
    let holder: SeriesStat | null = null;
    // Inside the ranked tier only: other depths do not compare.
    for (const stat of this.kit.rankableCores(data).filter((s) => s.tier === cores[0].tier)) {
      if (stat.best !== null && (holder === null || dir * stat.best.value > dir * holder.best!.value))
        holder = stat;
    }
    if (holder !== null && holder.matrix.id !== cores[0].matrix.id) {
      body.push(ui.divText(`${cores[0].matrix.label} runs highest on average. Your single best compound ` +
        `is on ${holder.matrix.label}, whose typical is ${host.formatActivity(holder.typical)}: one ` +
        'lucky substituent, not a better scaffold.', 'chem-sar-cp-hint'));
    }
    return this.card(title, subtitle, body);
  }
}
