/** The text of Plan's study protocol for an ethics application: the facts of the plan, already formatted where the
 *  protocol formats them (percentages, alpha, diameter), in the committee's wording. Pure: no platform imports. */
export interface EthicsFacts {
  /** one per treated arm; `effect` as the protocol writes it ('TGI 60%'), `power` a formatted percentage */
  arms: {name: string, therapy: string, effect: string, n: number | null, power: string}[];
  line: string;
  papers: number;
  /** the target power, a formatted percentage */
  target: string;
  alpha: string;
  alphaPer: string;
  correction: 'holm' | 'bonferroni' | 'none';
  nFixed: boolean;
  /** the minimum group size (MIN_GROUP) */
  minGroup: number;
  /** formatted percentages */
  preLoss: string;
  postLoss: string;
  confidence: string;
  lossRule: string;
  implant: number;
  randomise: number;
  /** null when the plan is not sized */
  nTotal: number | null;
  humane: number;
  diameter: string;
  lastDayChosen: boolean;
  end: number;
  endpoint: number;
}

export function ethicsText(f: EthicsFacts): string {
  const k = f.arms.length;
  const uniq = (xs: string[]) => [...new Set(xs)].join(' and ');
  const unreachable = f.arms.filter((a) => a.n === null);
  const same = f.arms.every((a) => a.n === f.arms[0].n && a.power === f.arms[0].power && a.effect === f.arms[0].effect);
  const nMax = Math.max(...f.arms.map((a) => a.n ?? 0));
  const control = k === 1 ? ` and ${nMax} controls` : ', with a control group as large as the largest treated group';
  const request = unreachable.length > 0 ?
    `No group size up to 40 reaches ${f.target} power for ${uniq(unreachable.map((a) => a.effect))}; the study is ` +
    'not adequately powered to detect it.' :
    same ? `We request ${nMax} analyzable animals per treated group${control}, which gives ${f.arms[0].power} ` +
      `power to detect ${f.arms[0].effect} at a two-sided significance level of ${f.alpha}.` :
      `We request ${f.arms.map((a, i) => `${a.n}${i === 0 ? ' analyzable animals' : ''} for ${a.name} ` +
        `(${a.power} power to detect ${a.effect})`).join(', ')}${control}, at a two-sided significance level of ` +
        `${f.alpha}.`;
  const size = unreachable.length > 0 ? '' : f.nFixed ?
    ' The group size was chosen by us; the power stated above is what this size provides.' :
    ` This is the smallest group size that reaches ${f.target} power, estimated by simulating tumor growth from the ` +
    `controls of ${f.papers} published ${f.line} studies. ` + (nMax > f.minGroup ?
      `Fewer animals per group would fall below ${f.target} power and risk an inconclusive study that may need to ` +
      'be repeated.' : 'We use five per group as a minimum, so that losing one animal still leaves a usable ' +
      'comparison.');
  const multiple = k < 2 ? '' : f.correction === 'none' ?
    ` Each of the ${k} comparisons is tested at α = ${f.alpha}, without correction for multiple comparisons.` :
    ` With ${k} comparisons, ${f.correction === 'holm' ? 'Holm' : 'Bonferroni'} correction is applied; each group ` +
    `is sized at α = ${f.alphaPer}.`;
  const arms = f.arms.map((a) => `${a.name} (${a.therapy}, expected ${a.effect})`).join(', ');
  const effects = [...new Set(f.arms.map((a) => a.effect))];
  const p1 = `${k + 1} groups: a control and ${arms}; each treated group is compared with the control. ${request} ` +
    (effects.length === 1 ? `[Explain why an effect of ${effects[0]} is expected or is the smallest effect of ` +
      'interest.]' : `[Explain why effects of ${effects.join(' and ')} are expected or are the smallest effects of ` +
      'interest.]') + `${size}${multiple}`;
  const lost = `About ${f.preLoss} of implanted animals are expected to be lost before randomization (no tumor take, ` +
    `tumor outside the volume window) and ${f.postLoss} after it (ulceration, deaths, technical failures).`;
  const animals = f.nTotal === null ? lost : `${lost} We therefore implant ${f.implant}, randomize ` +
    `${f.randomise} and expect ${f.nTotal} analyzable animals; ${f.lossRule === 'expected' ?
      'the counts cover the average loss' : `each count is enough with at least ${f.confidence} probability`}.`;
  const reduction = `Reduction: ${k === 2 ? 'one control group serves both treated groups; ' :
    k > 2 ? `one control group serves all ${k} treated groups; ` : ''}the size is ` +
    'based on the published growth of this tumor model, which avoids using animals in a pilot study for that ' +
    'purpose, and the extra animals only cover the expected losses. An animal reaches its humane endpoint when the ' +
    `tumor exceeds ${f.humane} mm³ (equivalent to a mean diameter of about ${f.diameter} cm); the study ends ` +
    (f.lastDayChosen ? `on day ${f.end}, chosen.` : `when the median control reaches ${f.endpoint} mm³, expected ` +
      `around day ${f.end}.`);
  return `${p1}\n\n${animals} ${reduction}`;
}
