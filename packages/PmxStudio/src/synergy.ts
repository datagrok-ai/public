/* Whether a combination does more than its two drugs would together: Bliss independence on tumour growth
   (Demidenko & Miller 2019, Crux stats.tumor_growth_synergy): the control, drug A, drug B and A+B, each group's
   growth rate of ln volume by a random-intercept model; synergy when A's and B's effects add up to less than the
   combination's (β_A + β_B − β_AB − β_C > 0), one-sided. Analyze's menu (Analyze > Combination synergy…). */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {crux, rows} from './crux';
import type {StudyInputs} from './explore';
import {pValue} from './ui';

export interface Synergy {slopes: {control: number, a: number, b: number, ab: number}, dif: number, se: number,
  z: number, p: number, pTwoSided: number}

/** The test over `groups`' measurements up to `lastDay`. */
export async function synergyOf(s: StudyInputs, groups: {control: string, a: string, b: string, ab: string},
  lastDay = Infinity): Promise<Synergy> {
  const [an, gr, dy, vol] = [s.animal, s.group, s.day, s.volume].map((c) => s.study.col(c)!);
  const wanted = new Set(Object.values(groups));
  const keep = Array.from({length: s.study.rowCount}, (_, i) => i).filter((i) => wanted.has(String(gr.get(i))) &&
    Number.isFinite(dy.get(i)) && dy.get(i) <= lastDay && vol.get(i) > 0);
  const r = await crux('stats.tumor_growth_synergy', {group: keep.map((i) => String(gr.get(i))),
    animal: keep.map((i) => String(an.get(i))), time: Float64Array.from(keep, (i) => dy.get(i)),
    volume: Float64Array.from(keep, (i) => vol.get(i))}, groups);
  const t = new Map(rows(r.main as any).map((x) => [x.term as string, x.estimate as number]));
  return {slopes: {control: t.get('slope_control')!, a: t.get('slope_a')!, b: t.get('slope_b')!,
    ab: t.get('slope_ab')!}, dif: t.get('dif')!, se: t.get('se')!, z: t.get('z')!, p: t.get('p_synergy')!,
  pTwoSided: t.get('p_two_sided')!};
}

/** What the test says, in a sentence. */
export function synergyWords(x: Synergy, ab: string, alpha = 0.05): string {
  return x.dif > 0 && x.p < alpha ? `Synergy: ${ab} grows slower than its two drugs' effects add up to (Bliss ` +
    `independence, one-sided p ${pValue(x.p)})` : x.dif < 0 && x.pTwoSided < alpha ? `Antagonism: ${ab} grows ` +
    `faster than its two drugs' effects add up to (p ${pValue(x.pTwoSided)})` : `No evidence of synergy: ${ab} is ` +
    `what its two drugs' effects add up to (Bliss independence, one-sided p ${pValue(x.p)})`;
}

/** Picks the groups (the combination guessed from its name: +, combo), runs the test and shows it. */
export function synergyDialog(s: StudyInputs, control: string, treated: string[], lastDay: number): void {
  if (treated.length < 3) {
    grok.shell.warning('A combination needs three treated groups: drug A, drug B and A+B');
    return;
  }
  const ab0 = treated.find((t) => /\+|combo|comb/i.test(t)) ?? treated[treated.length - 1];
  const [a0, b0] = treated.filter((t) => t !== ab0);
  const a = ui.input.choice('Drug A', {items: treated, value: a0});
  const b = ui.input.choice('Drug B', {items: treated, value: b0});
  const ab = ui.input.choice('A + B', {items: treated, value: ab0});
  ui.dialog('Combination synergy')
    .add(ui.divText(`Against ${control}, on the measurements to day ${lastDay}.`))
    .add(a).add(b).add(ab)
    .onOK(() => {
      const g = {control, a: a.value!, b: b.value!, ab: ab.value!};
      if (new Set([g.a, g.b, g.ab]).size < 3) {
        grok.shell.warning('Pick three different groups');
        return;
      }
      void synergyOf(s, g, lastDay).then((x) => {
        const rate = (v: number) => `${v.toFixed(3)} ln V per day`;
        ui.dialog('Combination synergy').add(ui.divText(synergyWords(x, g.ab)))
          .add(ui.tableFromMap({[`${control} (control)`]: rate(x.slopes.control), [g.a]: rate(x.slopes.a),
            [g.b]: rate(x.slopes.b), [g.ab]: rate(x.slopes.ab),
            'β_A + β_B − β_AB − β_C': `${x.dif.toFixed(3)} (SE ${x.se.toFixed(3)}, z ${x.z.toFixed(2)})`,
            'Method': 'Each group\'s growth rate of ln volume by a random-intercept linear mixed model (REML); ' +
              'Bliss independence on growth (Demidenko & Miller 2019)'})).show();
      }).catch((e) => grok.shell.error(String(e?.message ?? e)));
    }).show();
}
