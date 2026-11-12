// read from the file alone (no engine, no index, nothing computed)
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {parsePlan} from './plan-file';
import type {SavedPlan} from './plan-file';
import {groupsShownIn, parseStudy} from './study-file';
import type {SavedStudy} from './study-file';
import {parseMeta} from './study-meta';

/** `n` with its noun, singular or plural: "1 animal", "4 groups". */
const count = (n: number, noun: string) => `${n} ${noun}${n === 1 ? '' : 's'}`;

/** What a plan says, in a few lines: its model and end, each group with its animals, who saved it. */
export function planSummary(p: SavedPlan): string[] {
  const pr = p.predicted;
  // null means more than the search went to; a count absent (an older file) or not a number is not said
  const animals = (n: unknown) => n === null ? ', over 40 analyzable animals' :
    typeof n === 'number' && Number.isFinite(n) ? `, ${n} analyzable animals` : '';
  const saved = typeof p.created === 'string' ? `Saved on ${p.created.slice(0, 10)}` : 'Saved';
  return [
    `${pr.line} · endpoint ${pr.endpoint} mm³ · the study ends on day ${pr.studyEnd}`,
    `Control${animals(pr.nControl)}`,
    ...pr.arms.map((a) => `${a.name}: ${a.metric === 'TC' ? 'T/C' : 'TGI'} ${a.effect}% expected${animals(a.n)}`),
    `${saved}${typeof p.author === 'string' && p.author ? ` by ${p.author}` : ''}`,
  ];
}

/** What a study holds so far, from its saved rows: its animals, days measured, the values marked and its plan; its
 *  groups (names, number, the plan's name) only when the groups may be shown. */
export function studySummary(s: SavedStudy): string[] {
  const shown = groupsShownIn(s.meta, s.log);
  const days = [...new Set(s.rows.map((r) => r.day).filter((d): d is number => typeof d === 'number'))];
  const groups = [...new Set(s.rows.map((r) => r.group).filter((g): g is string => !!g))];
  const animalsOf = (g: string) => new Set(s.rows.filter((r) => r.group === g).map((r) => r.animal)).size;
  const [flagged, excluded] = [s.rows.filter((r) => r.flagged).length, s.rows.filter((r) => r.excluded).length];
  const animals = count(new Set(s.rows.map((r) => r.animal)).size, 'animal');
  return [
    ...parseMeta(s.meta)?.id ? [parseMeta(s.meta)!.id] : [],
    `${animals}${shown && groups.length ? ` in ${count(groups.length, 'group')}` : ''}, measured on ` +
      `${count(days.length, 'day')}${days.length ? `, the last day ${Math.max(...days)}` : ''}`,
    ...!shown ? ['Blinded: the groups are shown in Monitor once someone reveals them'] :
      groups.map((g) => `${g}: ${count(animalsOf(g), 'animal')}`),
    ...s.rows.some((r) => typeof r.weight === 'number') ? ['Body weight measured'] : [],
    ...flagged || excluded ? [`${count(flagged, 'value')} flagged to look at, ${excluded} excluded`] : [],
    !s.plan ? 'No plan attached' : shown && typeof s.plan.name === 'string' && s.plan.name ? `Plan: ${s.plan.name}` :
      'A plan attached',
  ];
}

/** The preview of `file`: its summary and the button that opens it, or a sentence on why it cannot be shown. */
export async function filePreview(file: DG.FileInfo, kind: 'plan' | 'study',
  open: (text: string) => unknown): Promise<DG.View> {
  const view = DG.View.create();
  view.name = file.name;
  let text: string;
  let lines: string[];
  try {
    text = await file.readAsString();
  } catch {
    view.append(ui.divText('The file could not be read.'));
    return view;
  }
  try {
    lines = kind === 'plan' ? planSummary(parsePlan(text)) : studySummary(parseStudy(text));
  } catch (e: any) {
    // the parsers say it plainly (damaged, newer, not ours); anything else is a file they let through
    const said = e instanceof TypeError ? null : String(e?.message ?? '');
    view.append(ui.divText(said || `This ${kind} file could not be read: open it in PMX Studio to see why.`));
    return view;
  }
  view.append(ui.divV([
    ui.h1(file.name.replace(/\.pmx(plan|study)$/i, '')),
    ...lines.map((l) => ui.divText(l)),
    ui.button(kind === 'plan' ? 'Open in Plan' : 'Open in Monitor', () => Promise.resolve().then(() => open(text))
      .catch((e) => grok.shell.error(`This ${kind} could not be opened: ${e?.message ?? e}`))),
  ], {style: {gap: '6px', padding: '12px'}}));
  return view;
}
