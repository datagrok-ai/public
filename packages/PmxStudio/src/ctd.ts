/* The regulatory table for nonclinical efficacy: "Primary pharmacodynamics: in vivo tumour models", one study a
   block of rows and each group a row, in the column style of ICH M4S(R2) Appendix B 2.6.3.4 (M4S gives no layout for
   primary pharmacodynamics), with the 2.6.3.1 overview line above it. Analyze gives it for one study, the Overview
   for all the studies of a compound. It takes each study's facts, not where they are stored. */
import type * as DG from 'datagrok-api/dg';

import {ArriveItem} from './arrive';
import {findCombinations} from './arm-name';
import {AnalysisFrames, BASES, intervalNote} from './analyze';
import {Block, Cell, fill} from './docx';
import {isFiniteNumber as has} from './engine';
import {GroupDetails, StudyMeta, rowsOf} from './study-meta';
import type {SavedStudy} from './study-file';
import {aucUsed, sig as sig3} from './exposure';
import {entered} from './sheet';
import type {NcaByArm} from './exposure';
import {modelName} from './priors';

/** A group's results: `tgi`, `tc` and their intervals in %, `p` against the control (adjusted as the analysis says);
 *  `nadir` the lowest mean body-weight change from day 0 (%), `lost` the animals dead or removed for their weight
 *  (null: no body weight). A control has no effect. */
export interface CtdGroup {
  name: string;
  control: boolean;
  /** compound, dose and schedule ('Cpd-X 10 mg/kg QD'); empty: not recorded */
  given: string;
  route: string;
  n: number;
  tgi?: number | null;
  tgiLower?: number | null;
  tgiUpper?: number | null;
  tc?: number | null;
  tcLower?: number | null;
  tcUpper?: number | null;
  /** the result of the test on the growth rates (analyze.ts Activity) */
  activity?: string | null;
  p?: number | null;
  /** complete and partial regressions (mRECIST), of `n` */
  regressions?: {complete: number, partial: number} | null;
  nadir: number | null;
  lost: number | null;
  /** NCA's exposure at the group's dose (AUC as the multiple uses it, Cmax), or why there is none; absent: no
   *  concentrations analyzed with the study */
  pk?: CtdPk | {why: string};
}

export interface CtdPk {auc: number, aucUnit: string, aucKind: string, cmax: number | null, cmaxUnit: string,
  profiles: number, pooled: boolean}

export interface CtdStudy {
  /** the study number (PMX-2026-014); empty: not given yet */
  id: string;
  name: string;
  compounds: string[];
  /** species, strain, sex: what is recorded */
  animals: string;
  line: string;
  endpointDay: number;
  /** the multiplicity correction, in words ('Holm over 2 comparisons'); empty: one comparison */
  correction: string;
  groups: CtdGroup[];
  /** each combination as indices into `groups` (A+B, A, B; arm-name.ts findCombinations), with where it stands
   *  against its drugs in words (Analyze's ladder; absent: not decided); absent: none, or kept before there were */
  combinations?: {ab: number, a: number, b: number, step?: string}[];
}

export const CTD_TITLE = 'Primary pharmacodynamics: in vivo tumour models';
export const M4S_URL = 'https://database.ich.org/sites/default/files/M4S_R2_Guideline.pdf';
export const CTD_CREDIT = 'Layout after ICH M4S(R2) Appendix B (2.6.3.1; columns in the style of 2.6.3.4, as M4S ' +
  'gives no layout for primary pharmacodynamics). GLP: No; primary pharmacodynamic studies need not be conducted ' +
  'under GLP (ICH S7A, 2.11). Prepared with PMX Studio.';

const r0 = (x: number) => String(Math.round(x));
const signed = (x: number) => `${x > 0 ? '+' : x < 0 ? '−' : ''}${Math.abs(x).toFixed(1)}%`;
const pText = (p: number) => p < 0.001 ? '<0.001' : p.toFixed(3);
const routes = (s: CtdStudy) => [...new Set(s.groups.filter((g) => !g.control).map((g) => g.route).filter((r) => r))]
  .join(', ');

/** A group's cells, after its study's two (said on its first row only). */
function groupCells(g: CtdGroup): Cell[] {
  const tgi = !has(g.tgi) ? '—' : r0(g.tgi) + (has(g.tgiLower) && has(g.tgiUpper) ?
    ` [${r0(g.tgiLower)}, ${r0(g.tgiUpper)}]` : '');
  const reg = g.regressions ? `${g.regressions.complete} CR, ${g.regressions.partial} PR` : '—';
  const bw = !has(g.nadir) ? '—' : `${signed(g.nadir)}; ${g.lost ?? 0}`;
  return [g.given || (g.control ? `${g.name} (control)` : fill(`${g.name}: dose and schedule`)),
    g.route || (g.control ? '' : fill('route')), String(g.n), g.control ? 'Control' : tgi,
    g.control || !has(g.tc) ? '—' : r0(g.tc), g.control || !has(g.p) ? '—' : pText(g.p),
    g.control ? '—' : reg, bw];
}

/** The overview line (2.6.3.1) of each study. */
export function overviewRows(studies: CtdStudy[]): Cell[][] {
  return [['Type of study', 'Test system', 'Method of administration', 'Testing facility', 'Study number',
    'Location (Vol./Section)'], ...studies.map((s) => [`Primary pharmacodynamics: antitumour activity, ${s.line}`,
    s.animals || fill('species, strain'), routes(s) || fill('route'), fill('Testing facility'),
    s.id || fill('Study number'), fill('Vol./Section')])];
}

/** The table: a header, then each study's groups, the study and its test system on its first row. */
export function ctdRows(studies: CtdStudy[]): Cell[][] {
  const pk = withPk(studies);
  const cells = (g: CtdGroup) => [...groupCells(g), ...pk ? [pkCell(g)] : []];
  return [['Study No.', 'Test system (species, strain, sex; model)', 'Dose and schedule', 'Route', 'n',
    'TGI % [95% CI]', 'T/C %', 'p', 'Regressions', 'BW nadir; deaths', ...pk ? ['AUC; Cmax'] : [], 'GLP',
    'Location'],
  ...studies.flatMap((s) => s.groups.map((g, i): Cell[] => i > 0 ? ['', '', ...cells(g), '', ''] :
    [s.id || fill('Study number'), `${s.animals ? `${s.animals}; ` : ''}${s.line} tumour model; endpoint day ` +
      `${s.endpointDay}`, ...cells(g), 'No', fill('Vol./Section')]))];
}

const withPk = (studies: CtdStudy[]) => studies.some((s) => s.groups.some((g) => g.pk));
const pkOf = (g: CtdGroup): CtdPk | null => g.pk && !('why' in g.pk) ? g.pk : null;
/** A group's AUC and Cmax; never one read from another dose: — with the reason in the notes. */
function pkCell(g: CtdGroup): string {
  const x = pkOf(g);
  return !x ? '—' : `${sig3(x.auc)} ${x.aucUnit}${x.cmax === null ? '' : `; ${sig3(x.cmax)} ${x.cmaxUnit}`}`;
}

/** The lowest exposure among the single-drug groups with a significant effect (p < 0.05): the minimal efficacious
 *  exposure, with its group and study; null: no such group with NCA's exposure. */
export function minimalEfficacious(studies: CtdStudy[]): {study: CtdStudy, group: CtdGroup, pk: CtdPk} | null {
  let best: {study: CtdStudy, group: CtdGroup, pk: CtdPk} | null = null;
  for (const study of studies) {
    for (const group of study.groups) {
      const pk = pkOf(group);
      if (pk && !group.control && has(group.p) && group.p < 0.05 && !group.given.includes(' + ') &&
        (!best || (pk.aucUnit === best.pk.aucUnit && pk.auc < best.pk.auc)))
        best = {study, group, pk};
    }
  }
  return best;
}

/** The PK table in the layout of M4S 2.6.5.3 (one column per group with exposure), what NCA could not give, and the
 *  minimal efficacious exposure; none without NCA. */
function pkBlocks(studies: CtdStudy[]): Block[] {
  if (!withPk(studies))
    return [];
  const cols = studies.flatMap((s) => s.groups.filter((g) => pkOf(g)).map((g) => ({s, g, x: pkOf(g)!})));
  const missing = studies.flatMap((s) => s.groups.filter((g) => g.pk && 'why' in g.pk)
    .map((g) => `${s.id || s.name}, ${g.name}: ${(g.pk as {why: string}).why}`));
  const row = (label: string, f: (c: typeof cols[0]) => Cell): Cell[] => [label, ...cols.map(f)];
  const m = minimalEfficacious(studies);
  const table: Block[] = !cols.length ? [] : [{head: true, rows: [
    row('Study No.', (c) => c.s.id || fill('Study number')),
    row('Species', (c) => c.s.animals || fill('species, strain')),
    row('Gender (M/F) / number of animals', (c) => [fill('sex'), `; ${c.x.profiles}` +
      (c.x.pooled ? ' (pooled, sparse sampling)' : '')]),
    row('Vehicle/formulation', () => fill('vehicle')),
    row('Method of administration', (c) => c.g.route || fill('route')),
    row('Dose and schedule', (c) => c.g.given || fill('dose')), row('Sample', () => fill('plasma, serum')),
    row('Analyte', (c) => c.s.compounds.join(', ') || fill('analyte')), row('Assay', () => fill('e.g. LC-MS/MS')),
    row('AUC', (c) => `${sig3(c.x.auc)} ${c.x.aucUnit} (${c.x.aucKind})`),
    row('Cmax', (c) => c.x.cmax === null ? '—' : `${sig3(c.x.cmax)} ${c.x.cmaxUnit}`)]}];
  return [{heading: 'Pharmacokinetics in the efficacy studies (layout of 2.6.5.3)', level: 2}, ...table,
    ...missing.length ? [{p: `No exposure given, none extrapolated: ${missing.join('; ')}.`, small: true}] : [],
    {p: m ? `Minimal efficacious exposure: ${m.group.given || m.group.name} (${m.study.id || m.study.name}, ` +
      `${m.study.line}), the lowest exposure with a significant effect (TGI ${has(m.group.tgi) ?
        `${r0(m.group.tgi)}%` : '—'}, p ${has(m.group.p) ? pText(m.group.p) : '—'}): AUC ${sig3(m.pk.auc)} ` +
      `${m.pk.aucUnit}${m.pk.cmax === null ? '' : `, Cmax ${sig3(m.pk.cmax)} ${m.pk.cmaxUnit}`}. A starting point ` +
      'for the pharmacologically active dose.' :
      'Minimal efficacious exposure: no group with NCA\'s exposure has a significant effect.'}];
}

/** The document: the test article, the overview line, the table and its notes. `drug`: a summary of one compound
 *  (the Overview's, over studies of it): that is the test article, the others only given with it in combination. */
export function ctdBlocks(studies: CtdStudy[], o: {sub: string, drug?: string}): Block[] {
  const compounds = [...new Set(studies.flatMap((s) => s.compounds))];
  const others = o.drug ? compounds.filter((c) => c.toLowerCase() !== o.drug!.toLowerCase()) : [];
  const article = o.drug ? `${o.drug}${others.length ? ` (given with ${list(others)} in combination)` : ''}` :
    compounds.join(', ');
  const corrections = [...new Set(studies.map((s) => s.correction).filter((c) => c))];
  return [{heading: CTD_TITLE}, {p: o.sub, small: true},
    {p: [{text: 'Test article: ', bold: true}, article || fill('INN or code')]},
    {heading: 'Pharmacology: overview (2.6.3.1)', level: 2},
    {rows: overviewRows(studies), head: true, widths: [3, 2.2, 1.6, 1.6, 1.4, 1.4]},
    {heading: 'Primary pharmacodynamics', level: 2},
    {rows: ctdRows(studies), head: true, widths: [1.3, 2.3, 2.2, 0.9, 0.5, 1.4, 0.7, 0.8, 1.2, 1.2,
      ...withPk(studies) ? [1.6] : [], 0.6, 1.4]},
    {p: 'TGI: tumour growth inhibition against the control on the endpoint day, with its 95% confidence interval. ' +
      'T/C: treated over control tumour volume on that day. p: Welch\'s t test of each group\'s tumour growth ' +
      `rate against the control${corrections.length ? `, ${corrections.join('; ')}` : ''}. Regressions: animals ` +
      'with a complete or partial response (mRECIST). BW nadir: the lowest mean body-weight change from day 0; ' +
      'deaths: animals dead or removed for their weight. —: not applicable or not recorded.', small: true},
    ...studies.flatMap((s) => {
      const note = intervalNote(s.groups.filter((g) => !g.control && g.activity).map((g) => ({activity: g.activity!,
        tcLower: g.tcLower ?? null, tcUpper: g.tcUpper ?? null})), s.endpointDay);
      return note ? [{p: studies.length > 1 ? `${s.name || s.id}: ${note}` : note, small: true}] : [];
    }),
    ...pkBlocks(studies),
    ...combinationBlocks(studies),
    {heading: 'For the Investigator\'s Brochure: nonclinical pharmacology', level: 2},
    ...ibParagraphs(studies, o.drug).map((p): Block => ({p})),
    {p: 'A draft from the tables above, to adapt.', small: true},
    {p: CTD_CREDIT, small: true}];
}

const tgiText = (g: CtdGroup) => !has(g.tgi) ? 'TGI not estimated' : `TGI ${r0(g.tgi)}%` +
  (has(g.tgiLower) && has(g.tgiUpper) ? ` (95% CI ${r0(g.tgiLower)} to ${r0(g.tgiUpper)})` : '') +
  (has(g.p) ? `, p ${pText(g.p)}` : '');
const label = (g: CtdGroup) => g.given || g.name;

/** Each combination against its two drugs alone, a row each: the study, the model, A+B, A and B with their effect
 *  against the control, and where A+B stands against them. */
export function combinationRows(studies: CtdStudy[]): Cell[][] {
  const cell = (g: CtdGroup) => `${label(g)}: ${tgiText(g)}`;
  return [['Study No.', 'Model', 'Combination', 'First drug alone', 'Second drug alone', 'Against its drugs'],
    ...studies.flatMap((s) => (s.combinations ?? []).map((c) => [s.id || fill('Study number'), s.line,
      cell(s.groups[c.ab]), cell(s.groups[c.a]), cell(s.groups[c.b]), c.step ?? '—']))];
}

/** The combinations' table (ICH S9: the rationale of a combination), none without one. */
function combinationBlocks(studies: CtdStudy[]): Block[] {
  const rows = combinationRows(studies);
  return rows.length < 2 ? [] : [{heading: 'Combinations: each against its drugs alone (ICH S9)', level: 2},
    {rows, head: true, widths: [1.2, 1, 2.4, 2.4, 2.4, 2]},
    {p: 'TGI and p against the control, as in the table above. Against its drugs: each animal\'s tumour growth ' +
      'under the combination against under each drug alone, from worse than the better drug alone to more than ' +
      'the two acting independently (Bliss); —: not decided, too few animals.', small: true}];
}

const list = (xs: string[]) => xs.length < 2 ? xs.join('') : `${xs.slice(0, -1).join(', ')} and ${xs.at(-1)}`;

/** The Investigator's Brochure's paragraphs on nonclinical pharmacology, from the studies' facts: the models, each
 *  single-drug group's effect, the combinations, tolerability, and the exposure when NCA gives it. Nothing beyond
 *  the tables: what is not recorded is not said. `drug`: the one compound summarized (ctdBlocks). */
export function ibParagraphs(studies: CtdStudy[], drug?: string): string[] {
  const compounds = drug ? [drug] : [...new Set(studies.flatMap((s) => s.compounds))];
  const species = [...new Set(studies.map((s) => s.animals.split(',')[0].trim()).filter((x) => x))];
  const lines = [...new Set(studies.map((s) => s.line))];
  const effects = studies.flatMap((s) => s.groups.filter((g) => !g.control && !g.given.includes(' + ') &&
    !(s.combinations ?? []).some((c) => s.groups[c.ab] === g))
    .map((g) => `${label(g)}${g.route ? ` ${g.route}` : ''} in ${s.line}: ${tgiText(g)} on day ${s.endpointDay}`));
  const out = [`The in vivo antitumour activity of ${list(compounds) || '[test article]'} was studied in ` +
    `${studies.length} ${studies.length === 1 ? 'study' : 'studies'} in ${lines.length} tumour ` +
    `model${lines.length === 1 ? '' : 's'} (${list(lines)})${species.length ? ` in the ${list(species)}` : ''}. ` +
    (effects.length ? `Against the control: ${effects.join('; ')}.` : '')];
  const combos = studies.flatMap((s) => (s.combinations ?? []).map((c) => `${label(s.groups[c.ab])} in ${s.line}: ` +
    `${tgiText(s.groups[c.ab])}${c.step ? `; ${c.step[0].toLowerCase()}${c.step.slice(1)}` : ''}`));
  if (combos.length)
    out.push(`In combination: ${combos.join('; ')}.`);
  const weighed = studies.flatMap((s) => s.groups.filter((g) => !g.control && has(g.nadir)).map((g) => ({s, g})));
  if (weighed.length) {
    const worst = weighed.reduce((a, x) => x.g.nadir! < a.g.nadir! ? x : a);
    const lost = weighed.reduce((a, x) => a + (x.g.lost ?? 0), 0);
    out.push(`The lowest mean body-weight change from day 0 in a treated group was ${signed(worst.g.nadir!)} ` +
      `(${label(worst.g)}, ${worst.s.line}); ${lost} treated animal${lost === 1 ? '' : 's'} died or ` +
      `${lost === 1 ? 'was' : 'were'} removed for weight loss.`);
  }
  const pk = studies.flatMap((s) => s.groups.filter((g) => pkOf(g)).map((g) => `${label(g)}: ${pkCell(g)}`));
  const m = minimalEfficacious(studies);
  if (pk.length) {
    out.push(`Exposure in the efficacy studies (AUC; Cmax): ${pk.join('; ')}.` + (m ? ` The lowest exposure with a ` +
      `significant effect was AUC ${sig3(m.pk.auc)} ${m.pk.aucUnit} (${label(m.group)}, ${m.study.line}).` : ''));
  }
  return out;
}

/** The dialog's rows for one study: what the table can say, and what the study's details or the user add. */
export function ctdItems(s: CtdStudy): ArriveItem[] {
  const treated = s.groups.filter((g) => !g.control);
  const item = (n: number, name: string, done: boolean, says: string, asks: string,
    fix?: 'details'): ArriveItem => ({n, item: name, done, status: done ? 'met' : 'todo', says, asks, fix});
  const given = treated.every((g) => g.given && g.route);
  const bw = s.groups.some((g) => has(g.nadir));
  return [
    item(1, 'Test article', s.compounds.length > 0, s.compounds.join(', ') || 'Add the compound of each group',
      'The INN or code of the compound', s.compounds.length ? undefined : 'details'),
    item(2, 'Test system', !!s.animals, s.animals ? `${s.animals}; ${s.line}` : 'Add the species, strain and sex',
      'Species, strain and sex; the tumour model', s.animals ? undefined : 'details'),
    item(3, 'Dose, schedule and route', given, given ? `${treated.length} treated groups` :
      'Add the dose, schedule and route of each group', 'Each group\'s dose (mg/kg), schedule and route',
    given ? undefined : 'details'),
    item(4, 'Efficacy', true, `TGI, T/C and p on day ${s.endpointDay}`, 'The noteworthy findings: TGI with its ' +
      'interval, T/C and p against the control; regressions'),
    item(5, 'Body weight', bw, bw ? 'Nadir and deaths per group' : 'No body weight in the study',
      'The lowest mean body-weight change and the deaths, the limited safety endpoints ICH S9 names'),
    item(6, 'Study number', !!s.id, s.id || 'Given when the study is saved; or write the report number',
      'The study or report number'),
    item(7, 'GLP', true, 'No: primary pharmacodynamics need not be GLP', 'The GLP status of each study'),
    item(8, 'Testing facility and location', false, 'Write them in the document',
      'Where the study was done, and its volume and section in the dossier'),
    ...!withPk([s]) ? [] : [pkItem(s)],
    ...!s.combinations?.length ? [] : [combinationItem(s)],
    {n: 11, item: 'Investigator\'s Brochure', done: true, status: 'met' as const, says: 'A paragraph from these facts',
      detail: ibParagraphs([s]).join('\n'), asks: 'The nonclinical pharmacology, in a few sentences'},
  // the app's own list, not a template's numbering: a row a study has no use for (exposure, combinations) leaves no
  // gap ("1–8, 10, 11")
  ].map((x, i): ArriveItem => ({...x, n: i + 1}));
}

/** The combinations' row: how many, and where they stand against their drugs. */
function combinationItem(s: CtdStudy): ArriveItem {
  const c = s.combinations!;
  const decided = c.filter((x) => x.step);
  return {n: 10, item: 'Combinations (ICH S9)', done: true, status: 'met',
    says: `${c.length} against ${c.length === 1 ? 'its' : 'their'} drugs alone` +
      (decided.length < c.length ? `; ${c.length - decided.length} not decided` : ''),
    detail: c.map((x) => `${s.groups[x.ab].name}: ${x.step ?? 'not decided'}`).join('\n'),
    asks: 'Each combination\'s effect against each of its drugs alone, the rationale ICH S9 asks for'};
}

/** The exposure row: how many groups NCA gives, the lowest efficacious one, what it could not give. */
function pkItem(s: CtdStudy): ArriveItem {
  const treated = s.groups.filter((g) => !g.control);
  const given = treated.filter((g) => pkOf(g)).length;
  const m = minimalEfficacious([s]);
  // met with any group's exposure: a group without one (a combination, a drug or dose NCA has no profile for) is
  // said in the detail, and the document never extrapolates it
  const done = given > 0;
  const why = treated.filter((g) => g.pk && 'why' in g.pk).map((g) => `${g.name}: ${(g.pk as {why: string}).why}`);
  return {n: 9, item: 'Exposure (NCA)', done, status: done ? 'met' : 'todo',
    says: `${given} of ${treated.length} groups with AUC and Cmax${m ? `; lowest efficacious: ${m.group.name}` : ''}`,
    asks: 'AUC and Cmax at each group\'s dose, from NCA; never extrapolated', detail: why.join('\n') || undefined};
}

/** A study's facts from its analysis: the groups' effects (comparisons), regressions (mRECIST per animal), body
 *  weight as weighed (tolerability), and what each group receives from the study's details. */
export function ctdStudyOf(r: AnalysisFrames, meta: StudyMeta | null, name: string, nca: NcaByArm | null = null,
  steps: Map<string, string> = new Map()): CtdStudy {
  const v = (c: string) => r.analysis.get(c, 0);
  const control = String(v('control'));
  const c = r.comparisons;
  const arms = Array.from({length: c.rowCount}, (_, i) => String(c.get('arm', i)));
  const weighed = new Map(frameRows(r.tolerability).filter((g) => g.basis === BASES.measured)
    .map((g) => [String(g.group), g]));
  const animals = frameRows(r.animals);
  const scored = animals.some((a) => typeof a.mrecist === 'string' && a.mrecist);
  const count = (group: string, re: RegExp) => animals.filter((a) => a.group === group && re.test(a.mrecist ?? ''))
    .length;
  const dose = (d: GroupDetails) => [d.compound.trim(), d.dose === null ? '' : `${entered(d.dose)} ${d.unit}`.trim(),
    d.schedule.trim()].filter((x) => x).join(' ');
  const group = (g: string, i: number | null): CtdGroup => {
    const details = rowsOf(meta, g).filter((d) => d.compound.trim());
    const w = weighed.get(g);
    const get = (k: string) => i === null ? null : c.get(k, i) as number | null;
    return {name: g, control: i === null, given: details.map(dose).join(' + '),
      route: [...new Set(details.map((d) => d.route.trim()).filter((x) => x))].join(' + '),
      n: i === null ? v('n_control') : c.get('n', i), tgi: get('tgi'), tgiLower: get('tgi_lower'),
      tgiUpper: get('tgi_upper'), tc: get('tc'), tcLower: get('tc_lower'), tcUpper: get('tc_upper'),
      activity: i === null ? null : c.get('activity', i) as string, p: get('p_value'),
      regressions: scored ? {complete: count(g, /^m?CR$/i), partial: count(g, /^m?PR$/i)} : null,
      nadir: w?.mean_nadir ?? null, lost: w ? (w.removed_death ?? 0) + (w.removed_weight ?? 0) : null,
      ...nca && i !== null ? {pk: ctdPk(nca, g)} : {}};
  };
  const k = arms.length;
  const corr = v('correction');
  const a = meta?.animals;
  return {id: meta?.id ?? '', name, line: modelName(String(v('line'))), endpointDay: v('evaluation_day'),
    compounds: [...new Set([control, ...arms].flatMap((g) => rowsOf(meta, g)).map((d) => d.compound.trim())
      .filter((x) => x))],
    animals: [a?.species, a?.strain, a?.sex].map((x) => x?.trim()).filter((x) => x).join(', '),
    correction: k < 2 ? '' : corr === 'none' ? 'not corrected for multiple comparisons' :
      `adjusted by ${corr === 'bonferroni' ? 'Bonferroni' : 'Holm'} over ${k} comparisons`,
    groups: [group(control, null), ...arms.map((g, i) => group(g, i))],
    // indices among the treated arms, +1 for the control first
    combinations: findCombinations(arms).map((x) => ({ab: x.ab + 1, a: x.a + 1, b: x.b + 1,
      ...steps.has(arms[x.ab]) ? {step: steps.get(arms[x.ab])} : {}}))};
}

/** A group's exposure as NCA matched it to its dose (exposureByArm), or why there is none. */
function ctdPk(nca: NcaByArm, arm: string): CtdPk | {why: string} {
  const a = nca.arms.get(arm);
  if (!a)
    return {why: nca.drugs?.has(arm) ? 'a combination; each drug\'s exposure is in Analyze' : 'not in NCA\'s results'};
  if (!a.exposure)
    return {why: a.why};
  const x = a.exposure;
  const u = aucUsed(x);
  return {auc: u.auc, aucUnit: u.unit, aucKind: u.inf ? 'AUCinf' : 'AUClast', cmax: x.cmax, cmaxUnit: x.cmaxUnit,
    profiles: x.profiles, pooled: x.pooled};
}

const frameRows = (df: DG.DataFrame): {[k: string]: any}[] => Array.from({length: df.rowCount}, (_, i) =>
  Object.fromEntries(df.columns.names().map((n) => [n, df.get(n, i)])));

/** A saved study's rows of the table, as its last analysis kept them, under the ID the study has now; null: not
 *  analyzed since the table was first kept. */
export function ctdOfSaved(s: SavedStudy): CtdStudy | null {
  const c = s.exposure?.effect?.ctd;
  return c && Array.isArray(c.groups) ? {...c, id: s.meta?.id || c.id} : null;
}

// ICH M4S's order of species in a summary
const SPECIES = ['mouse', 'rat', 'hamster', 'rabbit', 'dog'];
const speciesRank = (s: CtdStudy) => {
  const i = SPECIES.indexOf(s.animals.split(',')[0].trim().toLowerCase());
  return i < 0 ? SPECIES.length : i;
};
/** The studies in the order M4S lists them: by species, then by model and study number. */
export const m4sOrder = (studies: CtdStudy[]) => [...studies].sort((a, b) => speciesRank(a) - speciesRank(b) ||
  a.line.localeCompare(b.line) || a.id.localeCompare(b.id));

/** The drug every one of the studies gives (their drugs, comma-separated), the first if several; null: none. */
export function commonDrug(drugs: string[]): string | null {
  const sets = drugs.map((d) => d.split(',').map((x) => x.trim()).filter((x) => x));
  return sets[0]?.find((x) => sets.every((s) => s.some((y) => y.toLowerCase() === x.toLowerCase()))) ?? null;
}

/** The summary's dialog rows: each study, in or not (`left`: analyzed before its table was kept), then what the
 *  user writes in. */
export function summaryItems(included: CtdStudy[], left: string[]): ArriveItem[] {
  let n = 0;
  const item = (name: string, done: boolean, says: string, asks: string): ArriveItem =>
    ({n: ++n, item: name, done, status: done ? 'met' : 'todo', says, asks});
  return [
    ...included.map((s) => item(s.id || s.name, true, `${s.line}, ${s.groups.length - 1} treated ` +
      `group${s.groups.length === 2 ? '' : 's'}, day ${s.endpointDay}`, 'The study\'s rows')),
    ...left.map((name) => item(name, false, 'Open it in Analyze once to include it',
      'The analysis keeps the study\'s rows when it is shown')),
    item('Testing facility and location', false, 'Write them in the document',
      'Where each study was done, and its volume and section in the dossier'),
  ];
}
