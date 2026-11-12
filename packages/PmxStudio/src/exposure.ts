// exposure per group from the NCA package, when it is installed; nothing here is saved or enters the analysis
import * as DG from 'datagrok-api/dg';

import {drugKey, parseArmName} from './arm-name';

/** A group's exposure as NCA computed it: AUClast and Cmax (mean and CV%, as NCA labels and units them). */
export interface Exposure {
  auc: number, aucCv: number | null, aucUnit: string, cmax: number | null, cmaxCv: number | null, cmaxUnit: string,
  /** the dose and route NCA read, and the arm's own schedule: for the tooltip, never to scale */
  dose: string, route: string, schedule: string | null, profiles: number,
  /** the drug NCA names (the group's name without a drug column); a pooled profile of animals sampled sparsely */
  drug: string, pooled: boolean,
}
/** An arm: its exposure, or why it has none, in a plain sentence. */
export type ArmExposure = {exposure: Exposure, why?: undefined} | {exposure?: undefined, why: string};

const nca = (name: string) => DG.Func.find({package: 'Nca', name})[0] ?? null;
/** Whether the NCA package is installed (its functions are looked up, never imported). */
export const ncaInstalled = () => nca('ncaCompute') !== null;

// a column of the concentration table that names the group as the study does
const GROUP = /^(group|groups|arm|treatment|trt|grp|cohort|dosegroup)$/;
const norm = (s: string) => s.toLowerCase().replace(/[\s_-]+/g, '');
// NCA's routes (IV-bolus, IV-infusion, PO, SC, IM, other) and what an arm's name says
const ROUTE_OF: {[arm: string]: string} = {'i.v.': 'IV', 'p.o.': 'PO', 's.c.': 'SC', 'i.m.': 'IM', 'i.p.': 'other',
  'i.t.': 'other'};
const sameRoute = (arm: string | null, theirs: string) => arm === null ||
  (ROUTE_OF[arm] ?? arm) === (theirs.startsWith('IV') ? 'IV' : theirs);
const unitOf = (u: string) => u.replace(/[µμ]/g, 'u').toLowerCase();

/** NCA's mapping of `table`, or null when it does not read as concentrations over time; `issues`, what NCA missed. */
export async function concentrationMapping(table: DG.DataFrame, issues: string[] = []): Promise<any | null> {
  const detect = nca('ncaDetectMapping');
  if (!detect)
    return null;
  try {
    const d = JSON.parse(await detect.apply({table}));
    issues.push(...Array.isArray(d?.issues) ? d.issues : []);
    return d?.mapping && d.confidence === 'high' ? d.mapping : null;
  } catch {
    return null;
  }
}

/** NCA's per-profile results and per-group summary of `table`, its group column (if any) kept as a key. */
export async function runNca(table: DG.DataFrame, mapping: any):
Promise<{results: DG.DataFrame, summary: DG.DataFrame, groupCol: string | null}> {
  const used = new Set([mapping.subject, mapping.time, mapping.concentration, mapping.dose, mapping.route,
    ...mapping.groupBy ?? []]);
  const groupCol = table.columns.toList().find((c) => c.type === DG.COLUMN_TYPE.STRING && !used.has(c.name) &&
    GROUP.test(norm(c.name)))?.name ?? null;
  const m = {...mapping, groupBy: [...mapping.groupBy ?? [], ...groupCol ? [groupCol] : []]};
  const results: DG.DataFrame = await nca('ncaCompute')!.apply({table, mapping: JSON.stringify(m), rules: ''});
  const summary: DG.DataFrame = await nca('ncaSummarize')!.apply({results});
  return {results, summary, groupCol};
}

/** NCA's reasons for a group it could not compute, in plain words. */
function plain(reasons: string[]): string {
  const r = reasons.join('; ');
  if (/fewer than \d+ animals/.test(r))
    return 'too few animals sampled (a pooled profile needs at least 3)';
  if (/route column/.test(r))
    return 'NCA computes intravenous, oral, subcutaneous and intramuscular doses only, not this route';
  if (/mixed sampling/.test(r))
    return 'some animals sampled fully and others only a few times';
  if (/time points?|λz|lambda/i.test(r))
    return 'too few time points per animal';
  return r || 'no reason given';
}

/** Each arm's exposure from NCA's results and summary: the group of the same name, else (one drug only) the drug
 *  with the same dose and route. Never a guess: several candidates, or none, is said so. */
export function exposureByArm(arms: string[], results: DG.DataFrame, summary: DG.DataFrame,
  groupCol: string | null): Map<string, ArmExposure> {
  const keys: string[] = JSON.parse(summary.getTag('nca.keyColumns') || '[]');
  const units = JSON.parse(results.getTag('nca.units') || '{}');
  const doseUnit = typeof units.dose === 'string' ? unitOf(units.dose) : null;
  const drugCol = keys.find((k) => k !== groupCol && !/^(study|period|visit|occasion)/i.test(k)) ?? null;
  const str = (df: DG.DataFrame, c: string, i: number) => df.col(c)?.getString(i) ?? '';
  const num = (df: DG.DataFrame, c: string, i: number): number | null => {
    const col = df.col(c);
    return col && !col.isNone(i) && Number.isFinite(col.get(i)) ? col.get(i) : null;
  };
  // the rows of `df` that are this arm's, or why it cannot be matched; a row NCA blocked has no dose or route, so by
  // drug it is the drug's alone
  const rowsOf = (df: DG.DataFrame, arm: string, blocked = false): number[] | string => {
    const all = Array.from({length: df.rowCount}, (_, i) => i);
    if (groupCol)
      return all.filter((i) => str(df, groupCol, i) === arm);
    const parts = parseArmName(arm);
    if (parts.length !== 1)
      return 'a combination is matched by its exact name only (add a group column to the concentrations)';
    const [p] = parts;
    if (!drugCol)
      return 'the concentrations name no drug or group';
    if (!p.dose || doseUnit === null || unitOf(p.dose.unit) !== doseUnit)
      return 'its dose is not given in the same unit as NCA’s';
    return all.filter((i) => drugKey(str(df, drugCol, i)) === p.key && (blocked ?
      str(df, 'VERDICT', i) === 'Blocked' : sameRoute(p.route, str(df, 'ROUTE', i)) &&
      Math.abs((num(df, 'DOSE', i) ?? NaN) - p.dose!.value) <= 1e-9 * Math.max(1, p.dose!.value)));
  };
  const out = new Map<string, ArmExposure>();
  for (const arm of arms) {
    const no = (why: string) => out.set(arm, {why: `NCA could not compute exposure for ${arm}: ${why}`});
    const inSummary = rowsOf(summary, arm);
    if (typeof inSummary === 'string') {
      no(inSummary);
      continue;
    }
    if (inSummary.length > 1) {
      const sexes = new Set(inSummary.map((i) => str(summary, 'SEX', i)));
      no(sexes.size > 1 ? 'NCA gives its males and females apart' :
        'more than one NCA group fits it (by period or study)');
      continue;
    }
    if (inSummary.length === 0) {
      const inResults = rowsOf(results, arm) as number[];
      const blocked = groupCol ? inResults : rowsOf(results, arm, true) as number[];
      const reasons = blocked.map((i) => str(results, 'REASONS', i)).filter((r) => r);
      no(reasons.length ? `${groupCol ? '' : 'some doses of this drug were not computed, '}${plain(reasons)}` :
        inResults.length === 0 ? 'no concentrations found for it' : 'no reason given');
      continue;
    }
    const i = inSummary[0];
    const auc = num(summary, 'MEAN_AUC_LAST', i);
    if (auc === null) {
      no('no AUClast in NCA’s results');
      continue;
    }
    const unitTag = (c: string) => summary.col(c)?.getTag('units') ?? '';
    const parts = parseArmName(arm);
    out.set(arm, {exposure: {auc, aucCv: num(summary, 'CV_PCT_AUC_LAST', i), aucUnit: unitTag('MEAN_AUC_LAST'),
      cmax: num(summary, 'MEAN_CMAX', i), cmaxCv: num(summary, 'CV_PCT_CMAX', i), cmaxUnit: unitTag('MEAN_CMAX'),
      dose: `${num(summary, 'DOSE', i) ?? '?'}${doseUnit ? ` ${units.dose}` : ''}`, route: str(summary, 'ROUTE', i),
      schedule: parts.length === 1 ? parts[0].schedule : null, profiles: num(summary, 'N', i) ?? 0,
      drug: drugCol ? str(summary, drugCol, i) : arm, pooled: (rowsOf(results, arm) as number[])
        .some((j) => str(results, 'SUBJECT', j) === '(composite)')}});
  }
  return out;
}

/** "1,230 (CV 18%)": three significant figures; no CV for a single profile. */
export const meanCv = (mean: number, cv: number | null) =>
  `${Number(mean.toPrecision(3)).toLocaleString('en-US')}${cv === null ? '' : ` (CV ${Math.round(cv)}%)`}`;
