import * as ui from 'datagrok-api/ui';
import {RoleLevel, RoleSummary} from './sar-matrix-role-fit';
import {LOSER_ROWS, MIN_SUPPORT, SummaryData, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {count} from '../sar-matrix-ui-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRoleCards {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  roleCard(data: SummaryData, index: number): HTMLElement {
    return this.buildRoleCard(data, data.roleFit!.roles[index]);
  }

  /** So the card's compound count is not mistaken for the table's. */
  private roleDropped(data: SummaryData): HTMLElement[] {
    const dropped = data.roleFit?.dropped ?? 0;
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

  private buildRoleCard(data: SummaryData, role: RoleSummary): HTMLElement {
    const host = this.kit.host;
    const fit = data.roleFit!;
    const title = `${role.name} — effect vs the average compound`;
    const subtitle = `How far each value moves ${host.activityColumnName} from the library average of ` +
      `${host.formatActivity(fit.mean)}, the other components held fixed.`;
    const footer = this.roleFooter();
    const dropped = this.roleDropped(data);
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

    const chips: HTMLElement[] = [
      this.kit.badge(`spread ${role.spread.toFixed(2)}`, 'Count-weighted sd of this component\'s effects, ' +
        'with estimation noise subtracted so columns with different numbers of values compare. ' +
        'Not a range.'),
      this.kit.badge(`${role.levels.length} of ${role.fitted} values ranked`,
        `${role.thin} values seen fewer than ${MIN_SUPPORT} times stay in the fit but carry no readable ` +
        'effect of their own.'),
    ];
    if (role.repeat !== null) {
      chips.push(this.kit.badge(`repeats at r ${role.repeat.toFixed(2)}`,
        'Correlation of the effects fitted on each half of the table; below ' + TRUST_R2 + ' they are not ranked.'));
    }
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
      // The axis is the one role pinned in the Vary filter, so its row can land on the position too.
      const position = role.name === data.axisRole ? ref.matrix.columns[ref.ci].position : undefined;
      body.push(this.kit.summaryRow({
        depiction: level.value === '' ? null : this.kit.depiction(level.value),
        name: ref.matrix.label,
        badges: [this.roleCountBadge(role.name, level)],
        desc: `best measured ${host.formatActivity(ref.value)}` +
          (level.value === '' ? ' · nothing at this position' : ''),
        value: this.kit.formatEffect(level.coef),
        caption: 'vs average',
        valueTip: `Over ${count(level.n)} compounds, adjusted for ${others}. A difference from the average, ` +
          'not a potency.',
        faint,
        tip: role.name === data.coreRole ? `Open ${ref.matrix.label} in the SAR Matrix` :
          this.kit.openTip(ref.matrix, ref.ri, ref.ci),
        onClick: () => role.name === data.coreRole ? host.revealMatrix(ref.matrix) :
          host.revealCell(ref.matrix, ref.ri, ref.ci, position),
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
      this.kit.host.selectRoleValue(role, level.value);
    };
    return badge;
  }

  private roleFooter(): HTMLElement {
    return this.kit.hint('vs average: distance from the library average, the other components held fixed.',
      'Fitted over the whole table. Observational, not causal: neither a potency nor a prediction. The ' +
      'measured swap beside it counts pairs made side by side instead.');
  }
}
