/* The Overview's demo studies: a small lab's portfolio beside the user's own, each simulated from the atlas when it
   is opened (Crux, as the stages' example), so that every stage reads it as a real study. The MC38 study is there at
   its three moments, from one simulation: its plan (Plan), day 9 with today a measurement day (Monitor) and all its
   data (Analyze); the same animals, the same plan, the same day 0. Nothing is saved: no file is written. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {ExampleSpec, MC38_EXAMPLE, exampleDesign, exampleInputs, exampleStudy} from './analyze';
import {analyzeView} from './analyze-view';
import {LibraryEntry} from './library';
import {monitorView} from './monitor-view';
import {PLAN_TAG} from './plan-file';
import {planView} from './plan-view';
import {showStage} from './stages';
import {LOOKS_TAG, WEIGHT_TAG} from './study-file';
import {GroupDetails, StudyMeta, setMetaOn, treatmentOf} from './study-meta';

/** A demo study: what is simulated (`spec`), its details, and the stage it is listed at; `startIn`: day 0 from today
 *  for a plan (a followed or finished study's day 0 follows from its days). */
interface Demo {
  spec: ExampleSpec & {name: string};
  id: string;
  stage: 1 | 2 | 3;
  /** the treated groups' details, in the order of the spec's arms */
  treated: Omit<GroupDetails, 'group'>[];
  tags: string[];
  description: string;
  /** weeks are counted: day 0 is the Monday on or after `startIn` days */
  startIn?: number;
}

const cpdA = (dose: number) => ({compound: 'Cpd-A', dose, unit: 'mg/kg', route: 'i.p.', schedule: 'Q3D×4'});
const MC38: Omit<Demo, 'stage'> = {id: 'PMX-2026-014', tags: ['demo', 'dose response'],
  spec: {...MC38_EXAMPLE, name: 'MC38 · Cpd-A dose response'},
  treated: [cpdA(10), cpdA(30)],
  description: 'Cpd-A at two doses against vehicle in subcutaneous MC38: does the higher dose slow growth enough?'};

const DEMOS: Demo[] = [
  {...MC38, stage: 1},
  {...MC38, stage: 2},
  {...MC38, stage: 3},
  {id: 'PMX-2026-016', stage: 1, startIn: 10, tags: ['demo', 'KRAS'],
    spec: {name: 'CT26 · Cpd-B oral', line: 'CT26', seed: 2026016, arms: [
      {name: 'Cpd-B 25 mg/kg', therapy: 'small_molecule', metric: 'TGI', effect: 50},
      {name: 'Cpd-B 50 mg/kg', therapy: 'small_molecule', metric: 'TGI', effect: 75}]},
    treated: [25, 50].map((dose) => ({compound: 'Cpd-B', dose, unit: 'mg/kg', route: 'p.o.', schedule: 'QD×14'})),
    description: 'An oral KRAS inhibitor at two doses in CT26, before the study runs.'},
  {id: 'PMX-2026-012', stage: 2, tags: ['demo', 'IO'],
    spec: {name: 'B16F10 · anti-PD-1', line: 'B16F10', seed: 2026012, arms: [
      {name: 'anti-PD-1', therapy: 'immunotherapy', metric: 'TGI', effect: 60}]},
    treated: [{compound: 'anti-PD-1', dose: 10, unit: 'mg/kg', route: 'i.p.', schedule: 'Q3D×4'}],
    description: 'Checkpoint blockade in B16F10 melanoma, being followed.'},
  {id: 'PMX-2026-009', stage: 3, tags: ['demo', 'breast'],
    spec: {name: '4T1 · Cpd-C', line: '4T1', seed: 2026009, arms: [
      {name: 'Cpd-C 20 mg/kg', therapy: 'small_molecule', metric: 'TGI', effect: 65}]},
    treated: [{compound: 'Cpd-C', dose: 20, unit: 'mg/kg', route: 'p.o.', schedule: 'QD×21'}],
    description: 'A finished study of Cpd-C in orthotopic 4T1, ready to analyze.'},
];

const DAY_MS = 86_400_000;
const isoDay = (offset: number) => new Date(Date.now() + offset * DAY_MS).toISOString().slice(0, 10);
/** Today's weekday, Monday 0: the plans start on a Monday, so a study's day today keeps it. */
const weekday = () => (new Date().getDay() + 6) % 7;
/** A followed study's day today: the one nearest half its planned end that falls on today's weekday (day 0 a
 *  Monday, as the plan says). */
const halfway = (end: number) => {
  const wd = weekday();
  return Math.max(wd, Math.min(wd + 7 * Math.round((end / 2 - wd) / 7), Math.max(wd, end - 1)));
};
/** The Monday `weeks` weeks from this week's Monday. */
const monday = (weeks: number) => isoDay(7 * weeks - weekday());

const metaOf = (d: Demo, start: string | null): StudyMeta => ({id: d.id, description: d.description, owner: 'demo',
  project: 'Demo portfolio', tags: d.tags, start, notes: [], blinded: false,
  groups: [{group: 'Vehicle', compound: '', dose: null, unit: '', route: d.treated[0].route,
    schedule: d.treated[0].schedule}, ...d.spec.arms.map((a, k) => ({group: a.name, ...d.treated[k]}))]});

/** Opens a demo in its stage (followed: up to `day`), the simulation said in the task bar while it runs. */
async function open(d: Demo, day: number, start: string | null): Promise<void> {
  const progress = DG.TaskBarProgressIndicator.create(`Preparing the demo study ${d.spec.name}...`);
  try {
    const meta = metaOf(d, start);
    if (d.stage === 1) {
      showStage('plan', await planView({name: d.spec.name, study: meta, inputs: exampleInputs(d.spec)}));
      return;
    }
    const whole = await exampleStudy(d.spec);
    const df = d.stage === 3 ? whole : whole.clone(DG.BitSet.create(whole.rowCount, (i) => whole.get('day', i) <= day));
    df.name = d.stage === 3 ? d.spec.name : `${d.spec.name}, day ${day}`;
    df.setTag(PLAN_TAG, whole.getTag(PLAN_TAG)!);
    df.setTag(WEIGHT_TAG, 'body_weight');
    setMetaOn(df, meta);
    if (d.stage === 2) {
      // today a planned review, as Monitor's example
      df.setTag(LOOKS_TAG, String(day));
      showStage('monitor', await monitorView(df));
    } else
      showStage('analyze', await analyzeView({study: df}));
  } catch (e: any) {
    grok.shell.error(`The demo study could not be prepared: ${e?.message ?? e}`);
  } finally {
    progress.close();
  }
}

/** The demo studies as the Overview lists them, each with its design's end (a design per study, cached by the
 *  engine): a plan starting soon, followed studies on their day, finished ones. */
export async function demoEntries(): Promise<LibraryEntry[]> {
  return Promise.all(DEMOS.map(async (d, k): Promise<LibraryEntry> => {
    const end = (await exampleDesign(d.spec)).studyEnd;
    // a followed study is on its halfway day today; a finished one ended a week past its end, a week ago; the MC38
    // study's three moments share the followed one's day 0
    const chain = d.id === MC38.id;
    const day = d.stage === 1 ? 0 : d.stage === 2 ? halfway(end) : end + 7;
    const start = chain ? isoDay(-halfway(end)) : d.stage === 1 ? monday(Math.ceil((d.startIn ?? 7) / 7)) :
      d.stage === 2 ? isoDay(-day) : monday(-Math.ceil((day + 7) / 7));
    const name = chain && d.stage === 1 ? `${d.spec.name} (its plan)` :
      chain && d.stage === 3 ? `${d.spec.name} (all its data)` : d.spec.name;
    const meta = metaOf(d, start);
    return {kind: d.stage === 1 ? 'plan' : 'study', stage: d.stage, file: null, text: '', id: d.id, name,
      line: d.spec.line, compounds: d.treated.map((t) => t.compound).join(', '), owner: meta.owner, tags: d.tags,
      updated: Date.now() - k * 3 * DAY_MS, treatment: treatmentOf(meta), start, day: d.stage === 1 ? null : day,
      next: null, end, open: () => open(d, day, start)};
  }));
}
