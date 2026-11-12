/* A plan or a study file selected in Datagrok's Files: what it holds, read from the file alone (no engine, no
   index, nothing computed), and a button to open it in its stage. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {parsePlan} from './plan-file';
import type {SavedPlan} from './plan-file';
import {REVEALED, parseStudy} from './study-file';
import type {SavedStudy} from './study-file';

/** What a plan says, in a few lines: its model and end, each group with its animals, who saved it. */
export function planSummary(p: SavedPlan): string[] {
  const pr = p.predicted;
  const animals = (n: number | null | undefined) => n === null || n === undefined ? 'over 40' : `${n}`;
  return [
    `${pr.line} · endpoint ${pr.endpoint} mm³ · the study ends on day ${pr.studyEnd}`,
    `Control: ${animals(pr.nControl)} analyzable animals`,
    ...pr.arms.map((a) => `${a.name}: ${a.metric === 'TC' ? 'T/C' : 'TGI'} ${a.effect}% expected, ` +
      `${animals(a.n)} analyzable animals`),
    `Saved on ${p.created.slice(0, 10)}${p.author ? ` by ${p.author}` : ''}`,
  ];
}

/** What a study holds so far, from its saved rows: its animals, groups (not named while blinded), days measured,
 *  the values marked, its plan. */
export function studySummary(s: SavedStudy): string[] {
  const hidden = !!s.meta?.blinded && !s.log.some((e) => e.what === REVEALED);
  const days = [...new Set(s.rows.map((r) => r.day).filter((d): d is number => typeof d === 'number'))];
  const groups = [...new Set(s.rows.map((r) => r.group).filter((g): g is string => !!g))];
  const animalsOf = (g: string) => new Set(s.rows.filter((r) => r.group === g).map((r) => r.animal)).size;
  const [flagged, excluded] = [s.rows.filter((r) => r.flagged).length, s.rows.filter((r) => r.excluded).length];
  return [
    ...s.meta?.id ? [s.meta.id] : [],
    `${new Set(s.rows.map((r) => r.animal)).size} animals in ${groups.length} groups, measured on ${days.length} ` +
      `days${days.length ? `, the last day ${Math.max(...days)}` : ''}`,
    ...hidden ? ['Blinded: the groups are named in Monitor once someone reveals them'] :
      groups.map((g) => `${g}: ${animalsOf(g)} animals`),
    ...s.rows.some((r) => typeof r.weight === 'number') ? ['Body weight measured'] : [],
    ...flagged || excluded ? [`${flagged} values flagged to look at, ${excluded} excluded`] : [],
    s.plan ? `Plan: ${s.plan.name}` : 'No plan attached',
  ];
}

/** The preview of `file`: its summary and the button that opens it, or what keeps it from being read. */
export async function filePreview(file: DG.FileInfo, kind: 'plan' | 'study',
  open: (text: string) => unknown): Promise<DG.View> {
  const view = DG.View.create();
  view.name = file.name;
  const text = await file.readAsString();
  let lines: string[];
  try {
    lines = kind === 'plan' ? planSummary(parsePlan(text)) : studySummary(parseStudy(text));
  } catch (e: any) {
    view.append(ui.divText(String(e?.message ?? e))); // (the parsers say it plainly: damaged, newer, not ours)
    return view;
  }
  view.append(ui.divV([
    ui.h1(file.name.replace(/\.pmx(plan|study)$/i, '')),
    ...lines.map((l) => ui.divText(l)),
    ui.button(kind === 'plan' ? 'Open in Plan' : 'Open in Monitor', () => open(text)),
  ], {style: {gap: '6px', padding: '12px'}}));
  return view;
}
