/* Whether a combination does more than its two drugs would together: Bliss independence on tumour growth
   (Demidenko & Miller 2019, Crux stats.tumor_growth_synergy): the control, drug A, drug B and A+B, each group's
   growth rate of ln volume by a random-intercept model; synergy when A's and B's effects add up to less than the
   combination's (β_A + β_B − β_AB − β_C > 0), one-sided. Opened from the combination's line in Analyze's results. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {findCombinations} from './arm-name';
import {crux, named, rows} from './crux';
import {numbersOf, type StudyInputs} from './explore';
import {declarationOf, type PlanInputs} from './frames';
import {groupColors} from './colors';
import {metaOf} from './study-meta';
import {parsePlan} from './plan-file';
import {dialogLine, okNamed, pValue, showCentred, swatch} from './ui';

export interface Synergy {slopes: {control: number, a: number, b: number, ab: number}, dif: number, se: number,
  z: number, p: number, pTwoSided: number}

/** A combination's groups: the control, drug A, drug B and A+B. */
export type Combination = {control: string, a: string, b: string, ab: string};

/** The study's columns and the rows of `groups` (`keep`) from randomization to `lastDay` whose volume passes `ok`. */
function readRows(s: StudyInputs, groups: Combination, lastDay: number, ok: (volume: number) => boolean) {
  const [an, gr] = [s.animal, s.group].map((c) => s.study.col(c)!);
  const [dy, vol] = [s.day, s.volume].map((c) => numbersOf(s.study.col(c)!));
  const wanted = new Set(Object.values(groups));
  const keep = Array.from({length: s.study.rowCount}, (_, i) => i).filter((i) => wanted.has(String(gr.get(i))) &&
    Number.isFinite(dy[i]) && dy[i] >= 0 && dy[i] <= lastDay && ok(vol[i]));
  return {an, gr, dy, vol, keep};
}

/** The test over `groups`' measurements up to `lastDay`. */
export async function synergyOf(s: StudyInputs, groups: Combination, lastDay = Infinity): Promise<Synergy> {
  // from randomization on, as the analysis: growth before any dose would pull every rate together
  const {an, gr, dy, vol, keep} = readRows(s, groups, lastDay, (v) => v > 0);
  const r = await crux('stats.tumor_growth_synergy', {group: keep.map((i) => String(gr.get(i))),
    animal: keep.map((i) => String(an.get(i))), time: Float64Array.from(keep, (i) => dy[i]),
    volume: Float64Array.from(keep, (i) => vol[i])}, groups);
  const t = new Map(rows(r.main as any).map((x) => [x.term as string, x.estimate as number]));
  return {slopes: {control: t.get('slope_control')!, a: t.get('slope_a')!, b: t.get('slope_b')!,
    ab: t.get('slope_ab')!}, dif: t.get('dif')!, se: t.get('se')!, z: t.get('z')!, p: t.get('p_synergy')!,
  pTwoSided: t.get('p_two_sided')!};
}

/** The smallest volume Crux measures (its default floor, mm³): a tumor at or under it counts as gone. */
const FLOOR = 3.4;

/** Crux stats.tumor_combination_effects on `groups`' measurements from day 0 to `lastDay`, and the rows it read. */
async function effectsOf(s: StudyInputs, groups: Combination, beta: number, lastDay: number) {
  // a 0 kept: a tumor gone (Crux floors it, and it counts as gone)
  const {an, gr, dy, vol, keep} = readRows(s, groups, lastDay, (v) => v >= 0);
  // an animal is its group and its ID: studies number the animals in each group from 1, and Crux refuses an ID in
  // two groups
  const id = (i: number) => `${gr.get(i)}\u0000${an.get(i)}`;
  // a fixed seed: the bootstrap gives the same step each time the study is shown
  const r = await crux('stats.tumor_combination_effects', {animal: keep.map(id),
    group: keep.map((i) => String(gr.get(i))), day: Float64Array.from(keep, (i) => dy[i]),
    volume: Float64Array.from(keep, (i) => vol[i])}, {...groups, beta, start_day: 0, level: 0.975, seed: 20261124});
  return {r, rows: keep.map((i) => ({animal: id(i), group: String(gr.get(i)), day: dy[i], volume: vol[i]}))};
}

/** Where A+B's mean course would run if its two drugs acted independently (each animal the better of an A and a B
 *  animal, Palmer & Sorger 2017), day by day, as validated (analysis/out_combo_ladder/SEAL_band.md, REPORT_band.md,
 *  AUDIT_band.md): a 90 % band (87–90 % in the validation) for the mean of log V − log V(first day). `base`: A+B's
 *  mean log V on its first day (the band in mm³ is e^(base + value)). `verdict` on the last shared day: A+B's mean
 *  below the band kills more than independent action, above it less. Null with fewer than 3 shared days after the
 *  first or fewer than 3 animals in A, B or A+B. */
export interface IndependenceBand {days: number[], lower: number[], mid: number[], upper: number[], base: number,
  verdict: 'more' | 'less' | 'inside', day: number}

export async function independenceBand(s: StudyInputs, groups: Combination, beta: number, lastDay = Infinity):
  Promise<IndependenceBand | null> {
  const {r, rows: read} = await effectsOf(s, groups, beta, lastDay);
  // κ̂ by animal (effectsOf's IDs: group and animal)
  const kappa = new Map(named(r, 'animals').filter((x) => Number.isFinite(x.kappa))
    .map((x) => [x.animal as string, x.kappa as number]));
  // each animal's log V by day (the floor's at most), and its first day's
  const byAnimal = new Map<string, {group: string, at: Map<number, number>}>();
  for (const x of read) {
    const a = byAnimal.get(x.animal) ?? {group: x.group, at: new Map<number, number>()};
    a.at.set(x.day, Math.log(Math.max(x.volume, FLOOR)));
    byAnimal.set(x.animal, a);
  }
  const curves = (g: string) => [...byAnimal].filter(([, a]) => a.group === g).map(([id, a]) => {
    const first = Math.min(...a.at.keys());
    return {k: kappa.get(id) ?? -Infinity, base: a.at.get(first)!,
      at: new Map([...a.at].filter(([d]) => d > first).map(([d, v]) => [d, v - a.at.get(first)!]))};
  });
  const [ca, cb, cab] = [curves(groups.a), curves(groups.b), curves(groups.ab)];
  const n = Math.min(ca.length, cb.length, cab.length);
  const days = [...new Set(cab.flatMap((c) => [...c.at.keys()]))].sort((x, y) => x - y)
    .filter((d) => [ca, cb, cab].every((g) => g.some((c) => c.at.has(d))));
  if (n < 3 || days.length < 3)
    return null;
  const mean = (v: number[]) => v.reduce((t, x) => t + x, 0) / v.length;
  // seeded (mulberry32): the same band each time the study is shown
  let seed = 20261125;
  const u = () => {
    seed = (seed + 0x6D2B79F5) >>> 0;
    let t = Math.imul(seed ^ (seed >>> 15), seed | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  const pick = <T>(v: T[]) => v[Math.floor(u() * v.length)];
  const draws = days.map(() => [] as number[]);
  for (let b = 0; b < 1000; b++) {
    const A = ca.map(() => pick(ca));
    const B = cb.map(() => pick(cb));
    const pseudo = cab.map(() => {
      const [x, y] = [pick(A), pick(B)];
      return x.k >= y.k ? x : y;
    });
    days.forEach((d, j) => {
      const v = pseudo.map((c) => c.at.get(d)).filter((x): x is number => x !== undefined);
      if (v.length)
        draws[j].push(mean(v));
    });
  }
  // the percentile bootstrap widened for small groups: t(n − 1) / z and √(n / (n − 1))
  const z = 1.6449;
  const df = n - 1;
  const t = z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df * df);
  const f = (t / z) * Math.sqrt(n / df);
  const q = (v: number[], p: number) => {
    const sorted = [...v].sort((x, y) => x - y);
    const h = (sorted.length - 1) * p;
    const lo = Math.floor(h);
    return sorted[lo] + (h - lo) * (sorted[Math.min(lo + 1, sorted.length - 1)] - sorted[lo]);
  };
  const mid = draws.map((v) => q(v, 0.5));
  const lower = draws.map((v, j) => mid[j] - f * (mid[j] - q(v, 0.05)));
  const upper = draws.map((v, j) => mid[j] + f * (q(v, 0.95) - mid[j]));
  const last = days.length - 1;
  const seen = mean(cab.map((c) => c.at.get(days[last])).filter((x): x is number => x !== undefined));
  return {days, lower, mid, upper, base: mean(cab.map((c) => c.base)), day: days[last],
    verdict: seen < lower[last] ? 'more' : seen > upper[last] ? 'less' : 'inside'};
}

/** The band's verdict in words, for under Monitor's growth chart. */
export function bandWords(ab: string, x: IndependenceBand): string {
  return `Day ${x.day}: ${ab} ` + {more: 'slows growth more than its two drugs acting independently would',
    less: 'slows growth less than its two drugs acting independently would',
    inside: 'is where its two drugs acting independently would put it'}[x.verdict];
}

/** How the band was read, for its tooltip: what it is and what it cannot say. */
export function bandWhy(a: string, b: string): string {
  return `The dashed band: where the combination's mean would be if each animal responded as the better of an ` +
    `animal on ${a} and one on ${b}, from this study's animals (90 %, read day by day). Inside it is not proof that ` +
    'the drugs act independently, and with 10 animals or fewer per group a combination that slows growth ' +
    'more is often still inside.';
}

/** Where A+B stands against its two drugs, from worse than the better one alone to more than the two together. */
export type CombinationStep = 'antagonism' | 'no_added_benefit' | 'better_than_each' | 'beyond_independent_action';

/** A+B on the ladder (Crux stats.tumor_combination_effects: each animal's kill rate against the control's mean, at the
 *  line's `beta`), as validated (analysis/out_combo_ladder/SEAL_v2.md, REPORT_v2.md): intervals at 0.975, and "more
 *  than the two together" not said when an A+B tumor ends at the volume floor (beyond it the study cannot see).
 *  `floored`: an A+B tumor ended at the floor. Null: undecided (too few animals). */
export async function combinationStep(s: StudyInputs, groups: Combination, beta: number, lastDay = Infinity):
  Promise<{step: CombinationStep, floored: boolean} | null> {
  const {r} = await effectsOf(s, groups, beta, lastDay);
  const step = named(r, 'verdict')[0]?.step as string | undefined;
  if (!step || step === 'undecided')
    return null;
  const floored = named(r, 'animals').some((x) => x.role === 'ab' && x.cured === true);
  return floored && step === 'beyond_independent_action' ? {step: 'better_than_each', floored} :
    {step: step as CombinationStep, floored};
}

/** The step in words, for under A+B's result. */
export function stepWords(step: CombinationStep): string {
  return {antagonism: 'Worse than the better drug alone', no_added_benefit: 'No better than the better drug alone',
    better_than_each: 'Better than each drug alone',
    beyond_independent_action: 'More than its two drugs acting independently'}[step];
}

/** How the step was read, for its tooltip. */
export function stepWhy(a: string, b: string, floored: boolean): string {
  return `How fast each animal's tumor shrinks under the combination, against under ${a} and under ${b}.` +
    (floored ? ' Some combination tumors shrank to the smallest measurable volume: how much more than the two ' +
      'drugs together it does cannot be seen.' : '');
}

/** The α of the test: the study's plan's if one is attached (the plan's text), else `analysisAlpha`, Analyze's. */
export function synergyAlpha(planText: string | null | undefined, analysisAlpha: number): number {
  try {
    const a = planText?.trim() ? declarationOf(parsePlan(planText).inputs as unknown as PlanInputs).alpha : null;
    return a !== null && Number.isFinite(a) && a > 0 ? a : analysisAlpha;
  } catch {
    return analysisAlpha;
  }
}

/** The arms the plan's author unticked as combinations (the plan's text; none without a plan or one that reads). */
export function notCombined(planText: string | null | undefined): string[] {
  try {
    return planText?.trim() ? parsePlan(planText).inputs.arms.filter((a) => a.noCombination).map((a) => a.name) : [];
  } catch {
    return [];
  }
}

/** The test's verdict: more than the two drugs acting independently (synergy), less (antagonism), or about the
 *  same. The plan's convention for α, as `synergyWords`. */
export function synergyVerdict(x: Synergy, alpha: number): 'more' | 'about' | 'less' {
  return x.dif > 0 && x.p < alpha / 2 ? 'more' : x.dif < 0 && x.pTwoSided < alpha ? 'less' : 'about';
}

/** The test's result: the verdict in a line, each group's growth a day as a bar (the combination's with a mark where
 *  the two drugs acting independently would put it), the statistics folded under Details. */
function synergyResult(x: Synergy, g: Combination, alpha: number, colors: {[g: string]: string}): HTMLElement {
  const v = synergyVerdict(x, alpha);
  const [words, color] = {more: ['More than the two drugs acting independently', 'var(--green-2, #2e7d32)'],
    about: ['About what the two drugs give independently', 'var(--grey-5, #68727f)'],
    less: ['Less than the two drugs acting independently', 'var(--orange-2, #b26a00)']}[v];
  const p = v === 'less' ? x.pTwoSided : x.p;
  const verdict = ui.divH([ui.span(['●'], {style: {color}}), ui.span([words], {style: {fontWeight: '500'}}),
    ui.span([`p ${pValue(p)}`], {style: {color: 'var(--grey-5, #68727f)'}})], {style: {gap: '8px',
    alignItems: 'baseline', marginBottom: '12px'}});
  ui.tooltip.bind(verdict, `${g.ab} against ${g.a} and ${g.b} alone, each against ${g.control}. ` +
    (v === 'less' ? 'Two-sided p.' : 'One-sided p: the chance of a result this far beyond the two drugs if they ' +
    'only acted independently.'));
  // growth a day in %, the bars from 0 (a shrinking tumor's to the left)
  const pct = (slope: number) => 100 * (Math.exp(slope) - 1);
  const expected = x.slopes.a + x.slopes.b - x.slopes.control;
  const bars: [string, number][] = [[g.control, x.slopes.control], [g.a, x.slopes.a], [g.b, x.slopes.b],
    [g.ab, x.slopes.ab]];
  const all = [...bars.map(([, s]) => pct(s)), pct(expected)];
  const lo = Math.min(0, ...all);
  const span = Math.max(0, ...all) - lo || 1;
  const at = (q: number) => `${(100 * (q - lo) / span).toFixed(1)}%`;
  const grid = ui.div([], {style: {display: 'grid', gridTemplateColumns: 'max-content 220px max-content',
    gap: '6px 10px', alignItems: 'center'}});
  for (const [name, slope] of bars) {
    const q = pct(slope);
    const track = ui.div([], {style: {position: 'relative', height: '12px'}});
    // the zero axis, one line down the rows (the grid's gap is 6px): a shrinking tumor's bar lies left of it
    const zero = ui.div([], {style: {position: 'absolute', top: '-3px', bottom: '-3px', left: at(0), width: '0',
      borderLeft: '1px solid var(--grey-4, #9aa5b1)'}});
    ui.tooltip.bind(zero, 'No growth: right of it the tumors grow, left of it they shrink');
    track.append(zero, ui.div([], {style: {position: 'absolute', top: '0', bottom: '0', left: at(Math.min(0, q)),
      width: `${(100 * Math.abs(q) / span).toFixed(1)}%`, background: colors[name] ?? '#8a959f',
      borderRadius: '2px'}}));
    if (name === g.ab) {
      const mark = ui.div([], {style: {position: 'absolute', top: '-3px', bottom: '-3px', left: at(pct(expected)),
        width: '0', borderLeft: '2px dashed var(--grey-6, #2d3d4f)'}});
      ui.tooltip.bind(mark, `Expected if ${g.a} and ${g.b} acted independently: ${pct(expected).toFixed(1)} %/day`);
      track.append(mark);
    }
    grid.append(swatch(name === g.control ? `${name} (control)` : name, colors[name] ?? '#8a959f', name), track,
      ui.span([`${q.toFixed(1)} %/day`]));
  }
  const key = ui.divText(`┆ expected if independent: ${pct(expected).toFixed(1)} %/day`,
    {style: {color: 'var(--grey-5, #68727f)', marginTop: '6px'}});
  ui.tooltip.bind(key, 'Where the combination would grow if each drug slowed growth as it does alone, no more');
  const details = document.createElement('details');
  details.style.marginTop = '12px';
  const summary = document.createElement('summary');
  summary.textContent = 'Details';
  summary.style.cursor = 'pointer';
  details.append(summary, ui.tableFromMap({
    'Combination beyond its two drugs': `${x.dif.toFixed(3)} (SE ${x.se.toFixed(3)}, z ${x.z.toFixed(2)})`,
    'p (one-sided, two-sided)': `${pValue(x.p)}, ${pValue(x.pTwoSided)}`,
    'Method': 'Bliss independence on tumor growth (Demidenko & Miller 2019)'}));
  return ui.divV([verdict, grid, key, details]);
}

/** What the test says, in a sentence. The plan's convention for α: synergy when the one-sided p is under α / 2,
 *  antagonism when the two-sided p is under α. */
export function synergyWords(x: Synergy, ab: string, alpha: number): string {
  const v = synergyVerdict(x, alpha);
  return v === 'more' ? `Synergy: ${ab} grows slower than its two drugs' effects add up to (Bliss ` +
    `independence, one-sided p ${pValue(x.p)})` : v === 'less' ? `Antagonism: ${ab} grows ` +
    `faster than its two drugs' effects add up to (p ${pValue(x.pTwoSided)})` : `No clear synergy: ${ab} does ` +
    `not do clearly better than its two drugs together would (p ${pValue(x.p)})`;
}

/** Picks the groups (the combination found in the names of the groups, unless the plan unticked it: `notCombined`),
 *  runs the test and shows it, at `alpha`. */
export function synergyDialog(s: StudyInputs, control: string, treated: string[], lastDay: number,
  alpha: number, notCombined: string[] = []): void {
  if (treated.length < 3) {
    grok.shell.warning('A combination needs three treated groups: drug A, drug B and A+B');
    return;
  }
  const c = findCombinations(treated).find((x) => !notCombined.includes(treated[x.ab]));
  const [a0, b0, ab0] = c ? [treated[c.a], treated[c.b], treated[c.ab]] :
    [treated[0], treated[1], treated[treated.length - 1]];
  // the combination first (what is tested), then its two drugs alone
  const ab = ui.input.choice('Combination', {items: treated, value: ab0});
  const a = ui.input.choice('Drug A alone', {items: treated, value: a0});
  const b = ui.input.choice('Drug B alone', {items: treated, value: b0});
  const same = ui.divText('Pick three different groups', {style: {color: 'var(--red-3, #c62828)'}});
  const d = ui.dialog('Combination synergy (exploratory)')
    .add(dialogLine(`Against ${control}, to day ${lastDay}`, 'Exploratory: the plan did not set this test. ' +
      'Bliss independence on tumor growth (Demidenko & Miller 2019).'))
    .add(ui.inputs([ab, a, b])).add(same)
    .onOK(() => {
      const g = {control, a: a.value!, b: b.value!, ab: ab.value!};
      void synergyOf(s, g, lastDay).then((x) => {
        const r = ui.dialog('Combination synergy (exploratory)')
          .add(synergyResult(x, g, alpha, groupColors(treated, control, metaOf(s.study)?.colors,
            metaOf(s.study)?.slots))).onOK(() => {});
        okNamed(r, 'Close');
        // one way out: Close
        const cancel = r.getButton('CANCEL');
        if (cancel)
          cancel.style.display = 'none';
        showCentred(r);
      }).catch((e) => grok.shell.error(String(e?.message ?? e)));
    });
  // three different groups or no test: said under the fields, the button waits
  const ok = okNamed(d, 'Test');
  const check = () => {
    const bad = new Set([a.value, b.value, ab.value]).size < 3;
    same.style.display = bad ? '' : 'none';
    ok.disabled = bad;
  };
  [a, b, ab].forEach((x) => x.onChanged.subscribe(check));
  check();
  showCentred(d);
}
