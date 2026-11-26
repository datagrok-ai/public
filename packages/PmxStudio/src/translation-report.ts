/* The report's section on exposure and patients (Analyze's report, Compare's): what was measured and where each value
   comes from, then what it does not say. Pure: the content as text, built into the report by report.ts. */
import {ArmExposure, HumanExposure, HumanMultiple, aucUsed, fuOk, meanCv, multipleText, referenceLine, sig as num,
  sourceText} from './exposure';
import {htmlTable} from './report';
import {esc} from './sheet';
import {PanelSection, PanelVerdict, objectPanel, pValue, percentRange, textLines} from './ui';
import {gain} from './study-table';
import {isFiniteNumber} from './engine';
import {plainArm} from './arm-name';

/** One treated group's line. */
export interface ExposureRow {
  group: string, tgi: number | null, x: ArmExposure | undefined, multiple: HumanMultiple | null,
  /** the patients' exposure of its drug (undefined: not given) */
  human: HumanExposure | undefined,
  /** Compare: the study and model it comes from; a TGI not from the analysis (1 − T/C of medians) on `day` */
  study?: string, model?: string, medians?: {day: number | null},
}

export interface ExposureSection {
  sentence: string | null, warning: string | null, head: string[], rows: string[][], curve: string | null,
  notes: string[], sources: string[],
}

/** Where the patients' exposure of a drug came from (exposure.ts sourceText). */
export const humanSource = (h: HumanExposure | undefined) => h ? sourceText(h) : 'Not given';

/** Total or unbound, as the values were given. */
export function boundText(humans: (HumanExposure | undefined)[]): string {
  const given = humans.filter((h): h is HumanExposure => !!h);
  const fu = (h: HumanExposure) => fuOk(h.fuMouse) && fuOk(h.fuHuman);
  return given.length && given.every(fu) ? 'Exposures are unbound: both free fractions were given.' :
    given.some(fu) ? 'Exposures are unbound where both free fractions were given, total otherwise.' :
      'Exposures are total, not free: no free fractions were given.';
}

/** The mouse AUC of a group, with its unit and which AUC. */
function aucCell(x: ArmExposure | undefined): string {
  if (!x?.exposure)
    return '—';
  const a = aucUsed(x.exposure);
  return `${num(a.auc)} ${a.unit} (${a.inf ? 'AUCinf' : 'AUClast'})`;
}

/** The multiple, its arrow and its reason. */
function multipleCell(m: HumanMultiple | null): string {
  if (!m)
    return '—';
  if (m.ratio === undefined)
    return `— (${m.why})`;
  // (no dosing at weekends: the daily AUC is the dose's times 5/7, so the AUC column over the patients' is not it)
  return `${multipleText(m.ratio)}${m.range === 'above' ? ' ↑ above the patients\' range' :
    m.range === 'below' ? ' ↓ below the patients\' range' : ''}${m.weekly ? ', averaged over the week' : ''}`;
}

/** The section's content; null when no group has exposure (no section). `curve`: the descriptive curve was drawn;
 *  `table`: the NCA table's name (Analyze). */
export function exposureSection(o: {sentence: string | null, rows: ExposureRow[], curve: boolean,
  table?: string | null, warning?: string | null}): ExposureSection | null {
  const rows = o.rows;
  if (!rows.some((r) => r.x?.exposure || r.human))
    return null;
  const compare = rows.some((r) => r.study);
  const head = [...compare ? ['Study', 'Model'] : [], 'Group', 'TGI from the analysis',
    'Mouse AUC per dose', '× patients', 'Patients\' exposure from'];
  const lines = rows.map((r) => [...compare ? [r.study ?? '', r.model ?? ''] : [], r.group,
    r.tgi === null ? '—' : `${Math.round(r.tgi)}%${r.medians ? ' *' : ''}`, aucCell(r.x), multipleCell(r.multiple),
    humanSource(r.human)]);
  // each drug's patients' exposure once (in Compare, once per study: two studies may hold two values)
  const humans = [...new Map(rows.filter((r) => r.human).map((r) => [`${r.study ?? ''}\u0000${r.human!.drug}`,
    {h: r.human!, study: r.study}])).values()];
  const projected = [...new Set(humans.filter(({h}) => h.source === 'projected').map(({h}) => h.drug))];
  const days = [...new Set(rows.filter((r) => r.medians && r.tgi !== null).map((r) => r.medians!.day))];
  const notes = ['The mice and the patients are set side by side: no dose for patients and no exposure to aim for ' +
    'is drawn from these data.',
  boundText(rows.map((r) => r.human)),
  ...projected.length ? [`The patients' exposure of ${projected.join(', ')} is a projection, not a measurement.`] :
    [],
  ...o.curve ? ['The descriptive curve holds only between the doses tested.'] : [],
  ...days.length ? [`* The study kept no analysis: 1 − T/C of the median volumes on ` +
      `day${days.length > 1 ? 's' : ''} ${days.map((d) => d ?? '?').join(', ')}.`] : []];
  // a multiple averaged over the week (no dosing at weekends): its line, as its hover gives it on screen
  const weekly = rows.filter((r) => r.multiple?.weekly).map((r) => `${r.group}: ${r.multiple!.weekly}`);
  const sources = [
    ...weekly,
    ...!rows.some((r) => r.x?.exposure) ? [] : o.table !== undefined ?
      [`Mouse exposure: NCA, table ${o.table || 'not named'}.`] : ['Mouse exposure: NCA kept with each study.'],
    ...humans.map(({h, study}) => `Patients' exposure of ${h.drug}${study ? ` (${study})` : ''}: ${humanSource(h)}` +
      `${h.cv !== null && h.cv !== undefined ? `, CV ${num(h.cv)}%` : ''}.`)];
  return {sentence: o.sentence, warning: o.warning ?? null, head, rows: lines,
    curve: o.curve ? 'Within the tested doses only.' : null, notes, sources};
}

/** The section as an element of the report (report.ts reportHtml's annex). */
export function sectionElement(s: ExposureSection): HTMLElement {
  const el = document.createElement('div');
  el.innerHTML = [s.sentence ? `<p>${esc(s.sentence)}</p>` : '', s.warning ? `<p><i>${esc(s.warning)}</i></p>` : '',
    htmlTable(s.head, s.rows), s.curve ? `<p>${esc(s.curve)}</p>` : '',
    `<h3>What this does not say</h3><ul>${s.notes.map((n) => `<li>${esc(n)}</li>`).join('')}</ul>`,
    `<h3>Sources</h3><p>${s.sources.map(esc).join(' ')}</p>`].join('');
  return el;
}

/** A treated group's details for the context panel: its name and colour, what its name says it received (`read`), a
 *  few facts and the exposure as headed blocks (folded). */
export interface GroupDetails {title: string, color?: string, read: string, verdict?: PanelVerdict | null,
  facts: {label: string, value: string, tip?: string}[], blocks: {heading: string, lines: string[]}[]}

/** One drug of a combination group, for its details: NCA's exposure, the patients' and the multiple. */
export interface DrugDetails {drug: string, x: ArmExposure | undefined, multiple: HumanMultiple | null,
  human: HumanExposure | undefined}

/** What Analyze knows of one treated group, from what it has already computed (no new numbers). `read`: what the
 *  name and details say (study-meta.ts parseLine); `live`: NCA run now (else saved with the study); `curve`: where the
 *  group sits among the doses the descriptive curve was drawn through (null: no curve); `drugs`: a combination's,
 *  each with its own exposure blocks (then `x`, `multiple` and `human` are not read). From the groups table, when
 *  known: `n` analyzed, `tc` (T/C, %, and its 95% range), the `result`, a combination's `step` against its drugs
 *  (its words and why) and the `responders` ("3 of 8"). */
export function groupDetails(o: {group: string, read: string, tgi: number | null, day: number | null,
  x: ArmExposure | undefined, table: string | null, live: boolean, multiple: HumanMultiple | null,
  human: HumanExposure | undefined, curve: 'lowest' | 'highest' | 'between' | 'off' | null,
  drugs?: DrugDetails[], color?: string, n?: number | null,
  tc?: {value: number, lower: number | null, upper: number | null} | null, result?: string | null,
  step?: {words: string, why: string} | null, responders?: string | null,
  tgiCi?: {lower: number | null, upper: number | null} | null, nadir?: number | null, p?: number | null,
  removed?: string | null}): GroupDetails {
  const facts: GroupDetails['facts'] = [];
  // the answer first, in the groups table's words (its dot's colour), with its p-value
  const slows = o.result === 'Slows growth';
  const verdict: PanelVerdict | null = !o.result ? null : {text: o.result + (o.tgi !== null ?
    `: TGI ${Math.round(o.tgi)}%${o.day !== null ? ` on day ${o.day}` : ''}` : ''),
  // (the side the result says: below, above or not clearly apart)
  sub: isFiniteNumber(o.p) ? `${slows ? 'Clearly below' : o.result === 'Grows faster' ? 'Clearly above' :
    'Not clearly different from'} the control (p ${pValue(o.p).startsWith('<') ? pValue(o.p) : `= ${pValue(o.p)}`})` :
    undefined, tone: slows ? 'ok' : 'warn'};
  if (o.tgi !== null) {
    const tc = o.tc ? ` Its tumors ended at ${Math.round(o.tc.value)}% of the control's size` +
      (o.tc.lower !== null && o.tc.upper !== null ?
        ` (95% range ${percentRange(o.tc.lower, o.tc.upper)}).` : '.') : '';
    const ci = o.tgiCi && Number.isFinite(o.tgiCi.lower) && Number.isFinite(o.tgiCi.upper) ?
      ` (${percentRange(o.tgiCi.lower!, o.tgiCi.upper!)})` : '';
    facts.push({label: 'TGI', value: `${Math.round(o.tgi)}%${ci}${o.day !== null ? ` on day ${o.day}` : ''}`,
      tip: 'How much less the group grew than the control: 100% means its tumors did not grow; above 100%, they ' +
        `shrank.${tc}`});
  }
  // T/C as the table names it, its range quiet
  if (o.tc) {
    facts.push({label: 'Size vs control', value: `${Math.round(o.tc.value)}%` + (o.tc.lower !== null &&
      o.tc.upper !== null ? ` (${percentRange(o.tc.lower, o.tc.upper)})` : ''),
    tip: 'T/C: the group’s tumors as a share of the control’s on the evaluation day'});
  }
  if (o.step)
    facts.push({label: 'Against its drugs', value: o.step.words, tip: o.step.why});
  if (o.responders) {
    facts.push({label: 'Responders', value: o.responders + (o.day !== null ? ` measured on day ${o.day}` : ''),
      tip: 'Animals whose tumor ends below the range it would reach untreated'});
  }
  // the animals, and why some are missing from what was measured
  if (o.n !== undefined && o.n !== null)
    facts.push({label: 'Animals', value: String(o.n) + (o.removed ? ` · ${o.removed}` : '')});
  if (o.nadir !== undefined && o.nadir !== null && Number.isFinite(o.nadir)) {
    facts.push({label: 'Weight', value: `${gain(o.nadir)} at its lowest`,
      tip: 'The group’s lowest mean body-weight change from day 0 (its nadir)'});
  }
  const blocks: {heading: string, lines: string[]}[] = [];
  for (const d of o.drugs ?? [{drug: '', x: o.x, multiple: o.multiple, human: o.human}])
    blocks.push(...exposureBlocks(d, o.drugs ? `, ${d.drug}` : '', o.table, o.live));
  if (o.curve) {
    blocks.push({heading: 'Descriptive curve', lines: [{
      lowest: 'At the lowest dose tested: the curve does not reach below it',
      highest: 'At the highest dose tested: the curve does not reach above it',
      between: 'Between the lowest and highest dose tested, where the curve holds',
      off: 'Not on the curve: it is drawn through another set of groups'}[o.curve]]});
  }
  return {title: o.group, color: o.color, read: o.read, verdict, facts, blocks};
}

/** A group's (or one of its drugs', `of` in the headings) exposure, the patients' and the multiple, as blocks. */
function exposureBlocks(o: DrugDetails, of: string, table: string | null, live: boolean):
{heading: string, lines: string[]}[] {
  const blocks: {heading: string, lines: string[]}[] = [];
  const e = o.x?.exposure;
  if (e) {
    // the AUC the multiple uses, and why, is in its own explanation (Against the patients)
    blocks.push({heading: `Mouse exposure${of}`, lines: [`AUClast: ${meanCv(e.auc, e.aucCv)} ${e.aucUnit}`,
      ...e.aucInf !== null && e.aucInf > 0 ? [`AUCinf: ${num(e.aucInf)} ${e.aucInfUnit || e.aucUnit}`] : [],
      ...e.cmax !== null ? [`Cmax: ${meanCv(e.cmax, e.cmaxCv)} ${e.cmaxUnit}`] : [],
      `${e.profiles} profile${e.profiles === 1 ? '' : 's'} at ${e.dose}, route ${e.route}` +
        (e.pooled ? ', pooled from animals sampled a few times each' : ''),
      live ? `NCA, table ${table || 'not named'}` :
        `NCA, saved with the study (table ${table || 'not named'})`]});
  } else if (o.x?.why)
    blocks.push({heading: `Mouse exposure${of}`, lines: [o.x.why]});
  const h = o.human;
  if (h && h.auc !== null) {
    blocks.push({heading: `Patients' exposure${of}`, lines: [`AUC: ${num(h.auc)} ${h.unit}·h` +
      (h.cmax !== null ? `; Cmax: ${num(h.cmax)} ${h.unit}` : ''),
    h.reference ? referenceLine(h) : `Source: ${humanSource(h)}`,
    ...h.reference ? [`Matrix: ${h.reference.row.matrix.replace(/_/g, ' ')}`] : [],
    h.cv !== null && h.cv !== undefined ? `CV ${num(h.cv)}%` : 'CV not given: the range is the 3-fold convention']});
  }
  const m = o.multiple;
  if (m) {
    blocks.push({heading: `Against the patients${of}`, lines: m.ratio === undefined ? [m.why] :
      [multipleCell(m), m.tip]});
  }
  return blocks;
}

/** The details as the context panel shows them: name, facts, then `first` (the group's animals) and Exposure folded.
 *  Data only: how the numbers are computed stays on the facts' hover and in the report. */
export function detailsElement(d: GroupDetails, first: PanelSection[] = []): HTMLElement {
  const exposure = () => {
    const el = document.createElement('div');
    for (const b of d.blocks) {
      const head = document.createElement('div');
      head.textContent = b.heading;
      head.style.cssText = 'font-weight:600;margin-top:8px';
      el.append(head, textLines(b.lines));
    }
    return el;
  };
  // its treatment in plain words, the codes it was read from on hover
  const subtitle = d.read && d.read !== d.title ? plainArm(d.read) : undefined;
  return objectPanel({title: d.title, color: d.color, subtitle, subtitleTip: d.read || undefined,
    verdict: d.verdict, facts: d.facts,
    sections: [...first, ...d.blocks.length ? [{title: 'Exposure', content: exposure}] : []]});
}
