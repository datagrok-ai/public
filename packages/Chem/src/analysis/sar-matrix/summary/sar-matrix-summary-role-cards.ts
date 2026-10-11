/* The Effects segment's role cards: one fit over every component column, a card per component. */
import * as ui from 'datagrok-api/ui';
import {RoleLevel, RoleSummary} from './sar-matrix-role-fit';
import {LOSER_ROWS, MIN_SUPPORT, SummaryData, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {CARD_CORE_H, CARD_CORE_W, count} from '../sar-matrix-ui-common';
import {ROLE_SPREAD_TIE, PANE_SERIES, roleValueName, roleRanks} from './sar-matrix-summary-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRoleCards {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  /** One card per component, ordered by spread, or one card saying why nothing is ranked. Empty when the
   *  cores came from fragmentation. */
  build(data: SummaryData): {label: string, note: string | null, el: HTMLElement}[] {
    if (data.axisRole === null)
      return [];
    const refusal = this.kit.roleFitRefusal(data);
    if (refusal !== null) {
      return [{label: 'Components', note: null,
        el: this.kit.anchored(this.effects.card('Component contributions', refusal, this.roleDropped(data)),
          'role-0')}];
    }
    const roles = data.roleFit!.roles;
    return roles.map((role, i) => ({
      label: role.name,
      note: roleRanks(role) ? role.spread.toFixed(2) : null,
      el: this.kit.anchored(this.buildRoleCard(data, role), `role-${i}`),
    }));
  }

  /** Which component moves the endpoint most, above the cards. Stated as a rank, since the gaps are not
   *  stable; omitted where any component refuses to rank. */
  roleOrdering(data: SummaryData): HTMLElement | null {
    const roles = data.roleFit?.roles;
    if (roles === undefined || this.kit.roleFitRefusal(data) !== null || !roles.every(roleRanks))
      return null;
    // Components closer than the tie threshold are joined by ≈, not >.
    const chain = roles.map((role, i) => i === 0 ? role.name :
      `${roles[i - 1].spread - role.spread < ROLE_SPREAD_TIE ? '≈' : '>'} ${role.name}`).join(' ');
    const el = this.kit.hint(`Changing ${roles[0].name} moves ${this.kit.host.activityColumnName} most: ${chain}.`,
      'Each component\'s spread, corrected for its number of values so components compare. "≈" marks a ' +
      'gap too small to order.');
    return el;
  }

  /** The compounds the connected-block prune left out, so the card's count is not mistaken for the table. */
  private roleDropped(data: SummaryData): HTMLElement[] {
    const dropped = data.roleFit?.dropped ?? 0;
    if (dropped === 0)
      return [];
    return [ui.divText(`${count(dropped)} compounds share no component value with the rest and are not ` +
      'fitted.', 'chem-sar-sum-prose')];
  }

  /** The best three rows, then the worst two: one loser alone is often the least-supported value. */
  private roleRows(role: RoleSummary): {level: RoleLevel, faint: boolean}[] {
    const top = role.levels.slice(0, SUM_ROWS);
    return [...top.map((level) => ({level, faint: false})),
      ...role.levels.slice(Math.max(top.length, role.levels.length - LOSER_ROWS))
        .map((level) => ({level, faint: true}))];
  }

  /** Which component the matrix columns enumerate, and the control that changes it, so a card about R4
   *  beside a grid of R1 does not read as a disagreement. */
  private axisSwitchNote(name: string, data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const note = this.kit.hint(`Fitted over the whole table; the SAR Matrix columns show ${data.axisRole}.`,
      `Putting ${name} across the columns rebuilds every matrix.`);
    // The hover says the ranking itself will not change: it is one fit either way.
    const button = ui.button(`Put ${name} across the columns`, () => {
      if (!host.computing)
        host.setColumnAxis(name);
    }, `Rebuilds every matrix with ${name} across the columns. The offsets on this card stay the same.`);
    return ui.divH([note, button], 'chem-sar-sum-chips');
  }

  private buildRoleCard(data: SummaryData, role: RoleSummary): HTMLElement {
    const host = this.kit.host;
    const fit = data.roleFit!;
    const title = `${role.name} — offsets from the additive fit`;
    const subtitle = `Offsets vs the library mean of ${host.formatActivity(fit.mean)}, adjusted for the other ` +
      `components. R² ${fit.cvR2!.toFixed(2)} ± ${host.formatActivity(fit.cvRmse!)} on unseen compounds, ` +
      `n = ${count(fit.compounds)}.`;
    const footer = this.roleFooter();
    // On every card, since each tab shows one: a subtitle counting fewer compounds than the tab does
    // needs its explanation on the same card.
    const dropped = this.roleDropped(data);
    if (role.levels.length < 2) {
      return this.effects.card(title, subtitle, [...this.kit.reason(`${role.levels.length === 1 ? 'Only one' : 'No'} ` +
        `${role.name} value carries ${MIN_SUPPORT} or more compounds, so there is nothing to compare.`),
      ...dropped], footer);
    }
    // The halves check is the only test that this component's share of the credit means anything.
    if (role.repeat !== null && role.repeat < TRUST_R2) {
      return this.effects.card(title, subtitle, [...this.kit.reason(`${role.name} offsets do not repeat across ` +
        `halves of the table (r ${role.repeat.toFixed(2)}, below ${TRUST_R2}): too entangled with the other ` +
        'components to rank.'), ...dropped], footer);
    }

    const chips: HTMLElement[] = [
      this.kit.badge(`spread ${role.spread.toFixed(2)}`, 'Count-weighted sd of this component\'s offsets, ' +
        'with estimation noise subtracted so columns with different numbers of values compare. ' +
        'Not a range.'),
      this.kit.badge(`${role.levels.length} of ${role.fitted} values ranked`,
        `${role.thin} values seen fewer than ${MIN_SUPPORT} times stay in the fit but carry no readable ` +
        'offset of their own.'),
    ];
    if (role.repeat !== null) {
      chips.push(this.kit.badge(`repeats at r ${role.repeat.toFixed(2)}`,
        'Correlation of the offsets fitted on each half of the table; below ' + TRUST_R2 + ' they are not ranked.'));
    }
    const body: HTMLElement[] = [ui.divH(chips, 'chem-sar-sum-chips')];
    if (role.name !== data.axisRole && host.roleColumns.includes(role.name))
      body.push(this.axisSwitchNote(role.name, data));
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
        depiction: level.value === '' ? null : this.kit.depiction(level.value, CARD_CORE_W, CARD_CORE_H),
        name: roleValueName(level.value),
        badges: [this.roleCountBadge(role.name, level)],
        desc: `best measured ${host.formatActivity(ref.value)} · ${ref.matrix.label}`,
        value: this.kit.formatEffect(level.coef),
        caption: 'offset',
        valueTip: `Over ${count(level.n)} compounds, adjusted for ${others}. An offset, not a potency.`,
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

  /** The compound count, which selects those compounds so an offset can be checked against them. */
  private roleCountBadge(role: string, level: RoleLevel): HTMLElement {
    const badge = this.kit.badge(`${count(level.n)} cpd`, 'Measured compounds with this value. Click to select them.');
    badge.classList.add('chem-sar-sum-click');
    badge.onclick = (e: MouseEvent) => {
      e.stopPropagation();
      this.kit.host.selectRoleValue(role, level.value);
    };
    return badge;
  }

  /** The negative clause and the reconciliation, on every role card: each tab shows one card, so
   *  neither can be carried by a neighbour. */
  private roleFooter(): HTMLElement {
    return ui.divV([
      this.kit.prose('Offset = the mean difference of the compounds carrying this value, adjusted for the ' +
        'other components.',
      'Observational, not causal: neither a potency nor a prediction.'),
      this.kit.prose(`Fitted over the whole table. "${PANE_SERIES}" counts instead.`,
        `"${PANE_SERIES}" only sees pairs made side by side; this fit sees all of them.`),
    ]);
  }
}
