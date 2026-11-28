/* The PREPARE checklist (Smith et al., Lab Anim 2018), for a study being planned: its 15 topics under their three
   headings, which the plan already settles, which the study's details can record, and which are settled outside the
   app (with the facility, the ethics committee). Export > PREPARE checklist… in Plan, beside ARRIVE. */
import {ArriveItem, ChecklistActions, animalsItem, checklistDialog, receivesLine} from './arrive';
import type {StudyMeta} from './study-meta';

export const PREPARE_URL = 'https://norecopa.no/prepare/prepare-checklist';

/** What the plan says, for the topics it settles. `lastDay` the study's last day; `toImplant` the animals to
 *  implant, null when some arm reaches the target power with no group size. */
export interface PrepareFacts {line: string, protocol: string, humaneLimit: number, lastDay: number,
  toImplant: number | null, meta: StudyMeta | null}

const FORMULATION = 'Formulation of the study';
const DIALOGUE = 'Dialogue with the animal facility';
const QUALITY = 'Quality control of the components';

/** Each topic: its heading, a short name, what PREPARE asks for (its title first, then its points in short). */
const TOPICS: [string, string, string][] = [
  [FORMULATION, 'Literature searches', 'Literature searches: a clear hypothesis with primary and secondary ' +
    'outcomes, the searches and reviews to run, the species\' relevance, reproducibility and translatability'],
  [FORMULATION, 'Legal issues', 'Legal issues: the laws that apply (animal research, transport, health and ' +
    'safety) and their guidance documents'],
  [FORMULATION, 'Ethics and humane endpoints', 'Ethical issues, harm-benefit assessment and humane endpoints: a ' +
    'lay summary, the 3Rs, a harm-benefit assessment, a severity class, clear humane endpoints'],
  [FORMULATION, 'Design and statistics', 'Experimental design and statistical analysis: power and significance, ' +
    'the experimental unit and animal numbers, randomization, observer bias, inclusion and exclusion criteria'],
  [DIALOGUE, 'Timescale, funding and staff', 'Objectives and timescale, funding and division of labour: meetings ' +
    'with all staff, a timescale, the costs, who does what'],
  [DIALOGUE, 'Facility', 'Facility evaluation: an inspection of the building and equipment, staffing at times ' +
    'of extra risk'],
  [DIALOGUE, 'Education and training', 'Education and training: the staff\'s competence and any training ' +
    'needed before the study'],
  [DIALOGUE, 'Health risks and waste', 'Health risks, waste disposal and decontamination: a risk assessment with ' +
    'the facility, containment, decontamination and disposal'],
  [QUALITY, 'Test substances', 'Test substances and procedures: all that is known about the substances, and ' +
    'whether the procedures are feasible and valid'],
  [QUALITY, 'Experimental animals', 'Experimental animals: the characteristics essential for the study and its ' +
    'report; no surplus animals'],
  [QUALITY, 'Quarantine and health', 'Quarantine and health monitoring: the animals\' health status, transport, ' +
    'quarantine, monitoring'],
  [QUALITY, 'Housing and husbandry', 'Housing and husbandry: the animals\' needs, acclimatization, housing and ' +
    'any limits the study sets on them'],
  [QUALITY, 'Experimental procedures', 'Experimental procedures: refined ways to dose, sample, sedate and ' +
    'operate'],
  [QUALITY, 'Humane killing', 'Humane killing, release, reuse or rehoming: the primary and emergency methods, ' +
    'and who is competent to use them'],
  [QUALITY, 'Necropsy', 'Necropsy: a plan for every stage, where, and how animals and samples are identified'],
];

export function prepareItems(f: PrepareFacts): ArriveItem[] {
  const m = f.meta;
  const item = (n: number, status: ArriveItem['status'], says: string,
    more: {detail?: string, fix?: ArriveItem['fix']} = {}): ArriveItem => {
    const [group, name, asks] = TOPICS[n - 1];
    return {n, item: name, done: status === 'met', status, says, asks, group, ...more};
  };
  const receives = receivesLine(m?.groups ?? []);
  const animals = animalsItem(m, f.line);
  const surplus = f.toImplant === null ? '' : `The plan implants ${f.toImplant}, the loss it expects included`;
  return [
    item(1, 'todo', 'Search the literature; state the hypothesis'),
    item(2, 'todo', 'Check the laws and guidance that apply'),
    // the humane limit is one point of the topic; the rest is the ethics application's
    item(3, 'todo', `Humane limit ${f.humaneLimit} mm³ set; harm-benefit and severity to do`),
    // (an arm no group size powers is a design still to settle, though the protocol says it)
    !f.protocol ? item(4, 'todo', 'The plan is not ready yet') : f.toImplant === null ?
      item(4, 'todo', 'No group size reaches the target power for an arm: change its effect or the power',
        {detail: f.protocol, fix: 'protocol'}) :
      item(4, 'met', 'In the study protocol', {detail: f.protocol, fix: 'protocol'}),
    // day 0 to the last day, both counted: the days themselves, not a count that is one short
    item(5, 'todo', `${m?.start ? `From ${m.start} (day 0)` : 'From day 0'} to day ${f.lastDay}; agree costs and ` +
      'who does what'),
    item(6, 'todo', 'Inspect the rooms and equipment with the facility'),
    item(7, 'todo', 'Check that the staff are trained for each task'),
    item(8, 'todo', 'Assess risks and plan waste disposal with the facility'),
    receives ? item(9, 'met', receives.says, {detail: receives.detail}) :
      item(9, 'todo', 'Add what each group receives', {fix: 'details'}),
    animals.met ? item(10, 'met', animals.says, surplus ? {detail: surplus} : {}) :
      item(10, 'todo', animals.says, {fix: 'details'}),
    item(11, 'todo', 'Agree health status, transport and quarantine with the facility'),
    item(12, 'todo', 'Agree housing, acclimatization and care with the facility'),
    item(13, 'todo', 'Settle how dosing, measuring and sampling are done'),
    item(14, 'todo', 'Set the killing method and who performs it'),
    item(15, 'todo', 'Plan the necropsy and how samples are labelled'),
  ];
}

/** CC BY 4.0 asks for the source, the licence and what was changed. */
export const PREPARE_CREDIT = 'PREPARE: Smith AJ, Clutton RE, Lilley E, Hansen KEA, Brattelid T. Lab Anim ' +
  '2018;52(2):135-141, norecopa.no/prepare, CC BY 4.0; topics in short, not the original wording. Prepared with ' +
  'PMX Studio.';

/** The PREPARE checklist in a dialog; `details` opens the study's details, `protocol` copies the study protocol,
 *  `download` saves it as a document. */
export const prepareDialog = (items: ArriveItem[], actions: ChecklistActions, download?: () => Promise<void>) =>
  checklistDialog({title: 'PREPARE', guide: 'PREPARE guidelines', url: PREPARE_URL, asker: 'PREPARE', download},
    items, actions);
