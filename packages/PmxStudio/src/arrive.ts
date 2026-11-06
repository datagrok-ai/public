/* The ARRIVE 2.0 Essential 10, for the study analyzed: each item, whether the study already records it and what it
   says, or what is missing. Analyze's menu (Analyze > ARRIVE checklist…). */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import type {StudyMeta} from './study-meta';

export interface ArriveItem {n: number, item: string, done: boolean, says: string}

/** What the study records, for each item. */
export interface ArriveFacts {control: string, arms: string[], plan: string | null, exclusions: string | null,
  meta: StudyMeta | null, randomized: boolean, endpoint: number, evaluationDay: number, protocol: string,
  analyzed: {group: string, n: number}[]}

export function arriveItems(f: ArriveFacts): ArriveItem[] {
  const m = f.meta;
  const given = (m?.groups ?? []).filter((g) => g.compound.trim() !== '');
  const item = (n: number, name: string, done: boolean, says: string) => ({n, item: name, done, says});
  return [
    item(1, 'Study design', true, `${f.arms.length + 1} groups: ${f.control} (control), ${f.arms.join(', ')}; ` +
      'each animal is the experimental unit'),
    item(2, 'Sample size', !!f.plan, f.plan ? `Sized in Plan (${f.plan}): its protocol gives the power ` +
      'calculation' : 'No plan linked: say how the number of animals was decided'),
    item(3, 'Inclusion and exclusion criteria', true, f.exclusions ?? 'No animal or measurement excluded'),
    item(4, 'Randomisation', f.randomized, f.randomized ? 'Allocated to groups balanced on tumor volume at ' +
      'randomization (day 0)' : 'Say how the animals were allocated to the groups'),
    item(5, 'Blinding', !!m, !m ? 'Say who knew the groups when' : m.blinded ? 'Measurements entered blinded to ' +
      'the groups (the log says when they were revealed)' : 'Not blinded: say so, and why'),
    item(6, 'Outcome measures', true, `Tumor volume (mm³); endpoint ${f.endpoint} mm³, evaluated on day ` +
      `${f.evaluationDay}`),
    item(7, 'Statistical methods', f.protocol !== '', f.protocol || 'Analyze the study to write its methods'),
    item(8, 'Experimental animals', !!m?.description.trim(), m?.description.trim() ? m.description :
      'Add species, strain, sex, age, weight and source in the study’s details (its name)'),
    item(9, 'Experimental procedures', given.length > 0, given.length ? given.map((g) => `${g.group}: ` +
      `${g.compound}${g.dose === null ? '' : ` ${g.dose} ${g.unit}`}${g.route ? ` ${g.route}` : ''}` +
      `${g.schedule ? ` ${g.schedule}` : ''}`).join('; ') :
      'Add what each group receives in the study’s details (its name)'),
    item(10, 'Results', f.analyzed.length > 0, f.analyzed.map((a) => `${a.group}: ${a.n} analyzed`).join('; ')),
  ];
}

/** The checklist in a dialog, to copy. */
export function arriveDialog(items: ArriveItem[]): void {
  const text = items.map((x) => `${x.n}. ${x.item} [${x.done ? 'done' : 'to do'}]: ${x.says}`).join('\n');
  const done = items.filter((x) => x.done).length;
  const table = ui.table(items, (x) => [`${x.n}. ${x.item}`, x.done ? '✓' : 'To do', x.says],
    ['Item', '', 'What the study says']);
  table.style.maxWidth = '720px';
  ui.dialog(`ARRIVE Essential 10: ${done} of ${items.length} recorded`)
    .add(table)
    .addButton('Copy', () => void navigator.clipboard.writeText(text).then(() => grok.shell.info('Checklist copied')))
    .show();
}
