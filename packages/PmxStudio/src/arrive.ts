/* The ARRIVE 2.0 Essential 10 (and, on request, the Recommended Set), for the study analyzed: each item, whether the
   study already records it and what it says, or what is missing and where to add it. Export > ARRIVE checklist… in
   Analyze, with the author checklist as a document. The checklist dialog the other documents share lives here. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';

import {Block, fill} from './docx';
import {HUMANE} from './engine';
import {GroupDetails, StudyMeta, animalsLine} from './study-meta';
import {closeNamed, num, plainError, showCentred} from './ui';

/** An item: met, to do, or known only once the study is analyzed (`later`); `says` one short line, `detail` the full
 *  text, `asks` what the guideline asks for; `fix`: where the app fills the gap (`details`: the study's
 *  details; `protocol`: copies the study protocol); `group`: the heading it sits under, if the checklist has some. */
export interface ArriveItem {n: number, item: string, done: boolean, status: 'met' | 'todo' | 'later', says: string,
  detail?: string, asks: string, fix?: 'details' | 'protocol', group?: string}

/** What the study records, for each item. `exclusions` undefined: not known yet (a plan); `analyzed` empty: not
 *  analyzed yet; `correction` as the analysis names it (none, bonferroni, holm). */
export interface ArriveFacts {control: string, arms: string[], plan: string | null, exclusions?: string | null,
  meta: StudyMeta | null, randomized: boolean, endpoint: number, evaluationDay: number, protocol: string,
  alpha?: number, correction?: string, analyzed: {group: string, n: number}[],
  /** the humane limits: tumor volume (mm³, default HUMANE) and weight loss (percent, default 20), not the endpoint */
  humane?: number, weightLimit?: number}

export const ARRIVE_URL = 'https://arriveguidelines.org/arrive-guidelines';
const LATER = 'Once the study is analyzed';

/** What each item asks for, in short. */
const ASKS: string[] = [
  'The groups compared, the control among them, and what one experimental unit is (an animal, a cage, a litter)',
  'How many animals per group and in all, and how that number was decided',
  'The rules, set beforehand, for leaving animals or data out; what was left out, and why',
  'How animals were allocated to groups, and how order or cage effects were kept out',
  'Who knew which animal was in which group, at each step: allocation, dosing, measuring, analysis',
  'Each outcome measured and how, and which one is primary',
  'The statistical tests for each outcome, and how their assumptions were checked',
  'Species, strain, sex, age, weight and source of the animals',
  'What each group received: compound, dose, route, schedule, and when',
  'Per group: the summary with its spread, and the effect with its confidence interval',
];

/** What the groups receive in one line: each group's drugs (compound, dose, route, schedule) joined by ' + ', the
 *  groups by ' · ', a combination of the two drugs given alone 'both' ('LD4172 20 mg/kg p.o. QD · anti-PD-1 100 µg
 *  i.p. Q3D · both'); `detail` the full list, a group a line ("group: its drugs"; a group whose name already says
 *  what it receives, as a plan's arms do, by its name alone). null: no group names a compound. */
export function receivesLine(groups: GroupDetails[]): {says: string, detail: string} | null {
  const given = groups.filter((g) => g.compound.trim() !== '');
  if (!given.length)
    return null;
  const dose = (g: GroupDetails) => g.dose === null ? '' : `${g.dose} ${g.unit}`.trim();
  // each drug with its dose, route and schedule (ARRIVE item 9 asks for all three)
  const full = (g: GroupDetails) => [g.compound.trim(), dose(g), g.route.trim(), g.schedule.trim()]
    .filter((x) => x).join(' ');
  const byGroup = new Map<string, GroupDetails[]>();
  given.forEach((g) => byGroup.set(g.group, [...byGroup.get(g.group) ?? [], g]));
  const alone = new Set([...byGroup.values()].filter((r) => r.length === 1).map((r) => full(r[0])));
  // a combination of the two drugs given alone, as given alone: 'both'; else each drug in full (QA: "Trametinib +
  // anti-PD-1" without its doses beside the drugs' other doses)
  const says = [...byGroup.values()].map((r) => r.length === 2 && alone.size === 2 &&
    r.every((g) => alone.has(full(g))) ? 'both' : r.map(full).join(' + ')).join(' · ');
  // (a decimal comma in a name is the dose's point: "1,5 mg/kg" says 1.5 mg/kg)
  const norm = (s: string) => s.toLowerCase().replace(/(\d),(\d)/g, '$1.$2').replace(/\s+/g, ' ').trim();
  // a group a line: its name when it says what it receives, else what it receives after it, a combination's drugs
  // joined (QA: the NTS repeated each arm, and the combination once per drug)
  return {says, detail: [...byGroup].map(([group, rows]) => rows.every((g) => norm(group).includes(norm(full(g)))) ?
    group : `${group}: ${rows.map(full).join(' + ')}`).join('\n')};
}

/** The animals, for ARRIVE item 8 and PREPARE topic 10 (`prefix` first): met once every field is filled; the line
 *  says what is, and what to add. */
export function animalsItem(m: StudyMeta | null, prefix = ''): {met: boolean, says: string} {
  const {filled, missing} = animalsLine(m?.animals);
  const said = [prefix, filled].filter((x) => x).join(' · ');
  return {met: !missing.length, says: missing.length ? `${said ? `${said}: ` : ''}add ${missing.join(', ')}` : said};
}

/** The Recommended Set (items 11 to 21): its name and what it asks for. */
const RECOMMENDED: [string, string][] = [
  ['Abstract', 'A summary of the objectives, species, key methods, principal findings and conclusions'],
  ['Background', 'Why the study was done, its scientific background, and why this species and model'],
  ['Objectives', 'The questions, objectives and hypotheses, primary and secondary'],
  ['Ethical statement', 'The ethical review body that approved the work, and its reference numbers'],
  ['Housing and husbandry', 'The housing, husbandry and any enrichment'],
  ['Animal care and monitoring', 'Care and monitoring, how pain and distress were reduced, adverse events, and ' +
    'the humane endpoints'],
  ['Interpretation', 'The results against the objectives and earlier work, and the study\'s limitations'],
  ['Generalisability and translation', 'How far the findings may apply to other species, or to humans'],
  ['Protocol registration', 'Whether a protocol was written beforehand, and where it was registered'],
  ['Data access', 'Where the data can be found'],
  ['Declaration of interests', 'Competing interests, and the funding'],
];

/** The Essential 10, then with `full` the Recommended Set, each under its heading; the Recommended Set is the
 *  author's to write, but for what the study already records. */
export function arriveItems(f: ArriveFacts, full = false): ArriveItem[] {
  const e10 = essentials(f);
  if (!full)
    return e10;
  const rs = (n: number, says: string, detail?: string): ArriveItem => ({n, item: RECOMMENDED[n - 11][0],
    done: false, status: 'todo', says, asks: RECOMMENDED[n - 11][1], group: 'Recommended Set', detail});
  return [...e10.map((x) => ({...x, group: 'Essential 10'})),
    rs(11, 'Write it once the results are in'),
    rs(12, 'Say why the study was done, and why this model'),
    rs(13, 'State the hypothesis and the primary outcome'),
    rs(14, 'Name the ethics committee and the licence number'),
    rs(15, 'Describe the housing, husbandry and enrichment'),
    rs(16, `Humane endpoints: tumor over ${f.humane ?? HUMANE} mm³ or weight loss over ${f.weightLimit ?? 20}%; ` +
      'add the monitoring and any adverse events'),
    rs(17, 'Weigh the results against the objectives; say the limitations'),
    rs(18, 'Say how far the model speaks for patients'),
    f.plan ? rs(19, `Planned beforehand (${f.plan}): say where the plan was registered`,
      'Plan > Export > Study plan (ARRIVE) gives the protocol to register') :
      rs(19, 'Say whether a protocol was written beforehand, and where it was registered'),
    rs(20, 'Say where the data are; Export gives them as CSV'),
    rs(21, 'Declare competing interests and funding'),
  ];
}

function essentials(f: ArriveFacts): ArriveItem[] {
  const m = f.meta;
  const receives = receivesLine(m?.groups ?? []);
  const animals = animalsItem(m);
  const item = (n: number, name: string, status: ArriveItem['status'], says: string,
    more: {detail?: string, fix?: 'details'} = {}): ArriveItem =>
    ({n, item: name, done: status === 'met', status, says, asks: ASKS[n - 1], ...more});
  const k = f.arms.length;
  const correction = k < 2 ? '' : f.correction === 'none' ? ', no correction (exploratory)' :
    f.correction ? `, ${f.correction === 'bonferroni' ? 'Bonferroni' : 'Holm'} over ${k} comparisons` : '';
  const analyzedN = f.analyzed.reduce((a, x) => a + x.n, 0);
  return [
    item(1, 'Study design', 'met', `${k + 1} groups, each animal is the unit`,
      {detail: `${f.control} (control), ${f.arms.join(', ')}`}),
    f.plan ? item(2, 'Sample size', 'met', `Sized in Plan (${f.plan})`,
      {detail: 'The plan\'s protocol gives the power calculation'}) :
      item(2, 'Sample size', 'todo', 'No plan linked: say how the number of animals was decided'),
    f.exclusions === undefined ? item(3, 'Inclusion and exclusion criteria', 'later', LATER) :
      item(3, 'Inclusion and exclusion criteria', 'met', f.exclusions ? 'Some data left out' :
        'No animal or measurement excluded', f.exclusions ? {detail: f.exclusions} : {}),
    f.randomized ? item(4, 'Randomisation', 'met', 'Balanced on tumor volume at day 0') :
      item(4, 'Randomisation', 'todo', 'Say how the animals were allocated to the groups'),
    m?.blinded ? item(5, 'Blinding', 'met', 'Measurements entered blind to the groups',
      {detail: 'The log says when the groups were revealed'}) :
      item(5, 'Blinding', 'todo', 'Not blinded: turn on blinded entry, or say why not', {fix: 'details'}),
    item(6, 'Outcome measures', 'met', `Tumor volume; endpoint ${f.endpoint} mm³, evaluated on day ${f.evaluationDay}`),
    f.protocol ? item(7, 'Statistical methods', 'met', `Welch's t test on growth rates${correction}` +
      `${f.alpha === undefined ? '' : `, α ${num(f.alpha, 2)}`}`, {detail: f.protocol}) :
      item(7, 'Statistical methods', 'later', LATER),
    animals.met ? item(8, 'Experimental animals', 'met', animals.says) :
      item(8, 'Experimental animals', 'todo', animals.says, {fix: 'details'}),
    receives ? item(9, 'Experimental procedures', 'met', receives.says, {detail: receives.detail}) :
      item(9, 'Experimental procedures', 'todo', 'Add what each group receives', {fix: 'details'}),
    f.analyzed.length ? item(10, 'Results', 'met', `${analyzedN} animals analyzed in ${f.analyzed.length} groups`,
      {detail: f.analyzed.map((a) => `${a.group}: ${a.n}`).join('\n')}) : item(10, 'Results', 'later', LATER),
  ];
}

const MARK = {met: ['✓', 'var(--green-2, #2e7d32)', 'Met'], todo: ['○', 'var(--orange-2, #b26a00)', 'To do'],
  later: ['⋯', 'var(--grey-4, #8a959f)', 'Once the study is analyzed']} as const;

/** What a checklist's call-to-action does: open the study's details, or copy the study protocol. */
export type ChecklistActions = {[fix in NonNullable<ArriveItem['fix']>]?: () => void};
const FIX_LINKS = {details: ['Study details', 'Opens the study\'s details, where this is added'],
  protocol: ['Copy protocol', 'Copies the study protocol, the text for the ethics application']} as const;

/** The checklist in a dialog, to copy: one row per item (under its `group`'s heading, if any), its status, what the
 *  study says (the full text on hover) and, where the app can fill the gap, a link that does it. `download` saves
 *  it as a document (Download (.docx) beside Copy); `options` sits above the rows. */
export function checklistDialog(spec: {title: string, guide: string, url: string, asker: string,
  download?: () => void | Promise<void>, options?: HTMLElement},
items: ArriveItem[], actions: ChecklistActions = {}): DG.Dialog {
  let lastGroup: string | undefined;
  const text = items.map((x) => {
    const head = x.group && x.group !== lastGroup ? `${x.group}\n` : '';
    lastGroup = x.group;
    return `${head}${x.n}. ${x.item} [${MARK[x.status][2].toLowerCase()}]: ` +
      `${x.detail ? `${x.says}. ${x.detail.split('\n').join('; ')}` : x.says}`;
  }).join('\n');
  lastGroup = undefined;
  const done = items.filter((x) => x.done).length;
  const grid = ui.div([], {style: {display: 'grid', gridTemplateColumns: '16px max-content minmax(0, 1fr)',
    gap: '6px 12px', alignItems: 'baseline', maxWidth: '640px'}});
  const d = ui.dialog(`${spec.title} · ${done} of ${items.length}`);
  for (const x of items) {
    if (x.group && x.group !== lastGroup) {
      grid.append(ui.span([x.group], {style: {gridColumn: '1 / -1', fontWeight: '600',
        marginTop: lastGroup ? '6px' : '0'}}));
    }
    lastGroup = x.group;
    const [icon, color, tip] = MARK[x.status];
    const mark = ui.span([icon], {style: {color, fontWeight: '600', textAlign: 'center'}});
    ui.tooltip.bind(mark, tip);
    const name = ui.span([`${x.n}. ${x.item}`], {style: {whiteSpace: 'nowrap'}});
    ui.tooltip.bind(name, `${spec.asker} asks for: ${x.asks}`);
    const says = ui.span([x.says], {style: {overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap',
      color: x.status === 'met' ? '' : 'var(--grey-5, #68727f)'}});
    // the full text on hover: the detail, or a line too long to show
    ui.tooltip.bind(says, () => x.detail ?? (says.scrollWidth > says.clientWidth ? x.says : null));
    const cell = ui.divH([says], {style: {gap: '8px', minWidth: '0', alignItems: 'baseline'}});
    const act = x.fix && actions[x.fix];
    if (x.fix && act) {
      const [label, linkTip] = FIX_LINKS[x.fix];
      const link = ui.link(label, () => {
        if (x.fix === 'details')
          d.close();
        act();
      }, linkTip);
      link.style.whiteSpace = 'nowrap';
      cell.append(link);
    }
    grid.append(mark, name, cell);
  }
  const guide = ui.link(spec.guide, spec.url, 'The guidelines\' website');
  guide.target = '_blank';
  guide.style.marginTop = '12px';
  if (spec.options)
    d.add(spec.options);
  closeNamed(d).add(grid).add(guide)
    .addButton('Copy', () => void navigator.clipboard.writeText(text).then(() => grok.shell.info('Checklist copied')));
  if (spec.download) {
    const save = spec.download;
    d.addButton('Download (.docx)', () => void Promise.resolve().then(save)
      .catch((e) => grok.shell.error(plainError(e, 'The document could not be made. Try again.'))));
  }
  return showCentred(d);
}

/** A checklist as a document: a title, a line under it, then a table of item, what is asked, its status, what the
 *  study records and a column for the author (`last`: its header and the grey field), each group under its heading,
 *  and the credit at the end. */
export function checklistBlocks(o: {title: string, sub: string, items: ArriveItem[], credit: string,
  last: [string, string]}): Block[] {
  const out: Block[] = [{heading: o.title}, {p: o.sub, small: true}];
  for (const g of new Set(o.items.map((x) => x.group))) {
    if (g)
      out.push({heading: g, level: 2});
    out.push({head: true, widths: [1.2, 2, 0.6, 2, 1], rows: [['Item', 'Recommendation', 'Status',
      'What the study records', o.last[0]], ...o.items.filter((x) => x.group === g).map((x) => [`${x.n}. ${x.item}`,
      x.asks, MARK[x.status][2], x.detail ? `${x.says}\n${x.detail}` : x.says, fill(o.last[1])])]});
  }
  out.push({p: o.credit, small: true});
  return out;
}

export const ARRIVE_CREDIT = 'ARRIVE 2.0: Percie du Sert N et al., PLoS Biol 2020;18(7):e3000410, ' +
  'arriveguidelines.org (CC0); recommendations in short. Prepared with PMX Studio.';

/** The ARRIVE checklist in a dialog, the Essential 10 or, with `full`, ARRIVE 2.0 whole; `details` opens the study's
 *  details, `download` saves it as a document, `options` the switch between the two. */
export const arriveDialog = (items: ArriveItem[], o: {full?: boolean, details?: () => void,
  download?: () => Promise<void>, options?: HTMLElement} = {}) =>
  checklistDialog({title: o.full ? 'ARRIVE 2.0' : 'ARRIVE Essential 10', guide: 'ARRIVE 2.0 guidelines',
    url: ARRIVE_URL, asker: 'ARRIVE', download: o.download, options: o.options}, items, {details: o.details});
