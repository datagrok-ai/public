import * as ui from 'datagrok-api/ui';
import {RoleLevel, RoleSummary} from './sar-matrix-role-fit';
import {LOSER_ROWS, MIN_SUPPORT, SummaryData, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import {NOTHING_AT_POSITION} from './sar-matrix-summary-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRoleCards {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  /** So the card's compound count is not mistaken for the table's. */
  private roleDropped(dropped: number): HTMLElement[] {
    if (dropped === 0)
      return [];
    return [ui.divText(`${count(dropped)} compounds share no component value with the rest and are not ` +
      'fitted.', 'chem-sar-sum-prose')];
  }

  /** More than one loser: the single worst is often the least-supported value. */
  private roleRows(role: RoleSummary): {level: RoleLevel, faint: boolean}[] {
    const top = role.levels.slice(0, SUM_ROWS);
    return [...top.map((level) => ({level, faint: false})),
      ...role.levels.slice(Math.max(top.length, role.levels.length - LOSER_ROWS))
        .map((level) => ({level, faint: true}))];
  }

  roleCard(data: SummaryData, role: RoleSummary): HTMLElement {
    const host = this.kit.host;
    const fit = data.roleFit!;
    const title = `${role.name} — effect vs the average compound`;
    const subtitle = `How far each value moves ${host.activityColumnName} from the library average of ` +
      `${host.formatActivity(fit.mean)}, the other components held fixed.`;
    const footer = this.kit.hint('vs average: distance from the library average, the other components held fixed.',
      'Fitted over the whole table. Observational, not causal: neither a potency nor a prediction. The ' +
      'measured swap beside it counts pairs made side by side instead.');
    const dropped = this.roleDropped(fit.dropped);
    if (role.levels.length < 2) {
      return this.effects.card(title, subtitle, [...this.kit.reason(`${role.levels.length === 1 ? 'Only one' : 'No'} ` +
        `${role.name} value carries ${MIN_SUPPORT} or more compounds, so there is nothing to compare.`),
      ...dropped], footer);
    }
    // The halves check is the only test that this component's share of the credit means anything.
    if (role.repeat !== null && role.repeat < TRUST_R2) {
      return this.effects.card(title, subtitle, [...this.kit.reason(`${role.name} effects do not repeat across ` +
        `halves of the table (r ${role.repeat.toFixed(2)}, below ${TRUST_R2}): too entangled with the other ` +
        'components to rank.'), ...dropped], footer);
    }

    const chips = [
      this.kit.badge(`spread ${role.spread.toFixed(2)}`, `How far ${role.name} typically moves the activity, noise ` +
        'removed. Not a range.'),
      role.thin === 0 ?
        this.kit.badge(`${role.levels.length} values`, 'Every value has enough compounds to be ranked.') :
        this.kit.badge(`${role.levels.length} of ${role.fitted} ranked`, `${role.thin} values have fewer than ` +
          `${MIN_SUPPORT} compounds: kept in the fit, not ranked.`),
    ];
    const body: HTMLElement[] = [ui.divH(chips, 'chem-sar-sum-chips')];
    const split = role.split;
    if (split !== null) {
      body.push(ui.divText(`Bimodal: ${split.hiCount} values at ${this.kit.formatEffect(split.hiMean)} ` +
        `(n = ${count(split.hiN)}) and ${split.loCount} at ${this.kit.formatEffect(split.loMean)} ` +
        `(n = ${count(split.loN)}), a gap of ${split.gap.toFixed(2)} against the ` +
        `± ${host.formatActivity(fit.cvRmse!)} this fit resolves.`, 'chem-sar-sum-prose'));
    }

    const best = data.roleBest.get(role.name)!;
    const others = fit.roles.filter((other) => other !== role).map((other) => other.name).join(', ');
    for (const {level, faint} of this.roleRows(role)) {
      const ref = best.get(level.value)!;
      body.push(this.kit.summaryRow({
        depiction: level.value === '' ? null : this.kit.depiction(level.value),
        name: ref.matrix.label,
        badges: [this.roleCountBadge(role.name, level)],
        desc: `best measured ${host.formatActivity(ref.value)}` +
          (level.value === '' ? ` · ${NOTHING_AT_POSITION}` : ''),
        value: this.kit.formatEffect(level.coef),
        caption: 'vs average',
        valueTip: `Over ${count(level.n)} compounds, adjusted for ${others}. A difference from the average, ` +
          'not a potency.',
        faint,
        tip: role.name === host.coreRole ? this.kit.openMatrixTip(ref.matrix) :
          this.kit.openTip(ref.matrix, ref.ri, ref.ci),
        onClick: () => role.name === host.coreRole ? host.revealMatrix(ref.matrix) :
          host.revealCell(ref.matrix, ref.ri, ref.ci),
      }));
    }
    body.push(...dropped);
    return this.effects.card(title, subtitle, body, footer);
  }

  private roleCountBadge(role: string, level: RoleLevel): HTMLElement {
    const badge = this.kit.badge(`${count(level.n)} cpd`, 'Measured compounds with this value. Click to select them.');
    badge.classList.add('chem-sar-sum-click');
    badge.onclick = (e: MouseEvent) => {
      e.stopPropagation();
      this.kit.host.selectRows(level.molIdx);
    };
    return badge;
  }
}
