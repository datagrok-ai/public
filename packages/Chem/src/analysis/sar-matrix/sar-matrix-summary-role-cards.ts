/* The Effects segment's role cards: one fit over every component column, a card per component. */
import * as ui from 'datagrok-api/ui';
import {RoleLevel, RoleSummary} from './sar-matrix-role-fit';
import {LOSER_ROWS, MIN_SUPPORT, SummaryData, SUM_ROWS, TRUST_R2} from './sar-matrix-summary-data';
import {CARD_CORE_H, CARD_CORE_W, chipBadge, count} from './sar-matrix-ui-common';
import {ROLE_SPREAD_TIE, PANE_SERIES, roleValueName, roleRanks} from './sar-matrix-summary-common';
import type {SummaryEffects} from './sar-matrix-summary-effects';
import type {SummaryPanel} from './sar-matrix-summary-panel';

export class SummaryRoleCards {
  constructor(private readonly kit: SummaryPanel, private readonly effects: SummaryEffects) {}

  /**
   * One card per role column, ordered by that role's noise-corrected spread, or one card naming why
   * nothing is ranked.
   *
   * Absent rather than empty outside fragment-columns mode: a substituent label discovered by
   * fragmentation is local to its own series, so there is no scale on which one pooled offset could be
   * read and the question does not apply at all.
   */
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
    // One scale across every role card, against the tab's usual rule: each of these offsets comes from
    // one fit against one library average, so they are on one scale and a shared bar is the honest
    // rendering. Per-matrix effects elsewhere are centred inside their own matrix and are not. Over
    // the rows that will be drawn only — a coefficient on no card would shorten every bar on every
    // other one for a reason the reader cannot see.
    let scale = 0;
    for (const role of roles.filter(roleRanks)) {
      for (const {level} of this.roleRows(role))
        scale = Math.max(scale, Math.abs(level.coef));
    }
    return roles.map((role, i) => ({
      label: role.name,
      note: roleRanks(role) ? role.spread.toFixed(2) : null,
      el: this.kit.anchored(this.buildRoleCard(data, role, scale), `role-${i}`),
    }));
  }

  /** The one finding that is about every component at once, so it sits above the tabs rather than on
   *  any one of them. The order of the roles is stable; the size of the gap between two of them is not,
   *  so it is named as a rank and never quoted as a ratio. Omitted where a role is refusing to rank its
   *  own values, since that refusal is what invalidates the spread ordering. */
  roleOrdering(data: SummaryData): HTMLElement | null {
    const roles = data.roleFit?.roles;
    if (roles === undefined || this.kit.roleFitRefusal(data) !== null || !roles.every(roleRanks))
      return null;
    // Consecutive roles closer than the tie threshold are joined by an approximation sign rather than
    // an inequality: the order is stable, a gap that small is not, and one chain says so however many
    // component columns there are.
    const chain = roles.map((role, i) => i === 0 ? role.name :
      `${roles[i - 1].spread - role.spread < ROLE_SPREAD_TIE ? '≈' : '>'} ${role.name}`).join(' ');
    const el = this.kit.hint(`Changing ${roles[0].name} moves ${this.kit.host.activityColumnName} most: ${chain}.`,
      'Each component\'s range between its best and worst value, corrected for the fact that a column ' +
      'with more values gets a wider range by chance alone. That correction is what makes columns with ' +
      'different numbers of values comparable. "≈" marks a gap too small to order.');
    el.classList.add('chem-sar-sum-sub-lead');
    return el;
  }

  /** Stated wherever the fit is reported: the connected-block prune drops measured compounds by design,
   *  so the compound count the card prints is not the table. */
  private roleDropped(data: SummaryData): HTMLElement[] {
    const dropped = data.roleFit?.dropped ?? 0;
    if (dropped === 0)
      return [];
    return [ui.divText(`${count(dropped)} compounds share no component value with the rest of the table; ` +
      'only the largest connected block is fitted, because offsets from separate blocks rest on ' +
      'unrelated baselines.', 'chem-sar-sum-prose')];
  }

  /** The rows one role card shows: its best three, then its worst two. Two losers rather than one —
   *  where a role falls into two groups the single lowest value is routinely its least-supported
   *  member, and one row would then stand for a whole cluster by its smallest. */
  private roleRows(role: RoleSummary): {level: RoleLevel, faint: boolean}[] {
    const top = role.levels.slice(0, SUM_ROWS);
    return [...top.map((level) => ({level, faint: false})),
      ...role.levels.slice(Math.max(top.length, role.levels.length - LOSER_ROWS))
        .map((level) => ({level, faint: true}))];
  }

  /**
   * Which component the matrix columns actually enumerate, and the one control that changes it.
   *
   * A card ranking R4 sits beside a grid whose columns are R1, which reads as the two disagreeing. They
   * are answering different questions: the card is the fit over the whole table, the grid enumerates one
   * component at a time.
   */
  private axisSwitchNote(name: string, data: SummaryData): HTMLElement {
    const host = this.kit.host;
    const note = this.kit.hint('Ranked from the fit over the whole table. The SAR Matrix columns enumerate ' +
      `${data.axisRole}, so these values are not the ones across the top there.`,
    `Putting ${name} across the columns means rebuilding: every matrix is reassembled and every fit ` +
      'recomputed, which takes as long as the first run did.');
    // What will and will not move, on the control: this ranking is the same fit either way, so a reader
    // who clicks expecting these rows to change waits out a rebuild for a screen that looks identical.
    const pill = chipBadge(`Put ${name} across the columns — rebuilds`,
      `Reassembles every matrix with ${name} across the columns, and pools its measured pairs instead of ` +
      `${data.axisRole}'s. The offsets on this card do not change — they come from one fit over the whole ` +
      'table, whichever component the columns enumerate.', 'chem-sar-sum-role');
    pill.onclick = () => {
      if (!host.computing)
        host.setColumnAxis(name);
    };
    return ui.divH([note, pill], 'chem-sar-sum-chips');
  }

  private buildRoleCard(data: SummaryData, role: RoleSummary, scale: number): HTMLElement {
    const host = this.kit.host;
    const fit = data.roleFit!;
    const title = `${role.name} — offsets from the additive fit`;
    const subtitle = `Additive fit over ${fit.roles.length} component columns. Each offset is against ` +
      `the library mean of ${host.formatActivity(fit.mean)}, adjusted for the other components. ` +
      `R² ${fit.cvR2!.toFixed(2)} ± ${host.formatActivity(fit.cvRmse!)} predicting unseen compounds, ` +
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
    // Prediction is nearly unique even where the split of credit between the roles is not, so this is
    // the only check that this role's share of it means anything.
    if (role.repeat !== null && role.repeat < TRUST_R2) {
      return this.effects.card(title, subtitle, [...this.kit.reason(`${role.name} offsets do not repeat: fitted on ` +
        `each half of the table separately they correlate at r ${role.repeat.toFixed(2)}, below ` +
        `${TRUST_R2}. The other components in this table predict which ${role.name} a compound carries ` +
        'closely enough that the fit cannot separate their contributions on a subset, so these values ' +
        'are not ranked.'), ...dropped], footer);
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
        'Fitted on each half of the table separately, these offsets correlate this closely. Below ' +
        TRUST_R2 + ' they would not be ranked.'));
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
        mark: this.effects.effectBar(level.coef, scale, fit.cvRmse),
        value: this.kit.formatEffect(level.coef),
        caption: 'offset',
        valueTip: `${this.kit.formatEffect(level.coef)} against the library mean of ` +
          `${host.formatActivity(fit.mean)}, over ${count(level.n)} compounds, adjusted for ${others}. ` +
          'An offset, not a potency.',
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

  /** The compound count, which also selects those compounds. An offset is a statement about a subset
   *  of the table, and the only way to check one is to look at the subset it is about. */
  private roleCountBadge(role: string, level: RoleLevel): HTMLElement {
    const badge = this.kit.badge(`${count(level.n)} cpd`, 'Measured compounds carrying this value — click ' +
      'to select them in the table. No ± is printed: one from this count alone would ignore the ' +
      'degrees of freedom the other components cost.');
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
      'Observational, not causal: it is what the compounds that happen to carry this value did, so it ' +
        'is neither a potency nor a prediction.'),
      this.kit.prose(`Fitted over the whole table. The "${PANE_SERIES}" tab counts instead.`,
        `A count on "${PANE_SERIES}" only sees comparisons that were actually made inside one row — a ` +
        'substitution nobody tried side by side is invisible there, and visible here.'),
    ]);
  }
}
