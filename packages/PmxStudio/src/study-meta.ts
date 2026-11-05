/* The study as a thing of its own, from its plan to its report: an ID (PMX-2026-014), a description, who owns it, its
   project and tags, the day it starts (day 0, randomization), and what each group receives (compound, dose, route,
   schedule). It is written once in Plan and travels with the study: in the .pmxplan (SavedPlan.study), on the
   measurements table (META_TAG) and in the .pmxstudy (SavedStudy.meta). The Overview lists studies by it. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** What one group receives. */
export interface GroupDetails {
  group: string, compound: string, dose: number | null, unit: string, route: string, schedule: string,
}

export interface StudyMeta {
  /** PMX-<year>-<number>, given when the study is first saved or followed (empty until then) */
  id: string;
  description: string;
  /** the login of who runs it */
  owner: string;
  project: string;
  tags: string[];
  /** day 0 (randomization), YYYY-MM-DD; null: not fixed yet */
  start: string | null;
  groups: GroupDetails[];
}

/** The tag of a measurements table that carries its study's details (JSON). */
export const META_TAG = 'pmx.meta';
/** The tag that carries the plan (plan-file.ts's PLAN_TAG; here to read its details without a cycle). */
const PLAN_TAG = 'pmx.plan';

export const emptyMeta = (): StudyMeta => ({id: '', description: '', owner: grok.shell.user?.login ?? '',
  project: '', tags: [], start: null, groups: []});

/** Details read from a file or a tag someone may have edited: each field checked, the rest from emptyMeta. */
export function parseMeta(x: unknown): StudyMeta | null {
  if (!x || typeof x !== 'object')
    return null;
  const m = x as Partial<StudyMeta>;
  const text = (v: unknown) => typeof v === 'string' ? v : '';
  return {...emptyMeta(), id: text(m.id), description: text(m.description), owner: text(m.owner), project:
    text(m.project), tags: Array.isArray(m.tags) ? m.tags.filter((t) => typeof t === 'string') : [],
  start: typeof m.start === 'string' && /^\d{4}-\d{2}-\d{2}$/.test(m.start) ? m.start : null,
  groups: Array.isArray(m.groups) ? m.groups.filter((g) => g && typeof g.group === 'string').map((g) => ({
    group: g.group, compound: text(g.compound), dose: typeof g.dose === 'number' && Number.isFinite(g.dose) ?
      g.dose : null, unit: text(g.unit), route: text(g.route), schedule: text(g.schedule)})) : []};
}

/** A table's study details: its own, else its plan's (null: none yet). */
export function metaOf(df: DG.DataFrame | null | undefined): StudyMeta | null {
  const read = (tag: string, pick: (x: any) => unknown) => {
    try {
      return parseMeta(pick(JSON.parse(df?.getTag(tag) || 'null')));
    } catch {
      return null;
    }
  };
  return read(META_TAG, (x) => x) ?? read(PLAN_TAG, (p) => p?.study);
}

export const setMetaOn = (df: DG.DataFrame, m: StudyMeta) => df.setTag(META_TAG, JSON.stringify(m));

/** The next study ID of this user this year: PMX-2026-001, -002, ... ponytail: the counter is the user's own
 *  (grok.userSettings); a lab-wide sequence needs a shared counter (package settings) once studies are shared. */
export function nextStudyId(): string {
  const year = new Date().getFullYear();
  const key = `lastStudyId-${year}`;
  const n = (Number(grok.userSettings.getValue('PmxStudio', key)) || 0) + 1;
  grok.userSettings.add('PmxStudio', key, String(n));
  return `PMX-${year}-${String(n).padStart(3, '0')}`;
}

/** `m` with an ID (a new one when it has none). */
export const withId = (m: StudyMeta): StudyMeta => m.id ? m : {...m, id: nextStudyId()};

/** The groups' details in the order of `groups` (the study's groups now), each kept from `m` by name. */
export function groupsOf(m: StudyMeta, groups: string[]): GroupDetails[] {
  return groups.map((group) => m.groups.find((g) => g.group === group) ??
    {group, compound: '', dose: null, unit: 'mg/kg', route: '', schedule: ''});
}

/** One line on what the groups receive, for a list: the compounds, each once. */
export const compoundsOf = (m: StudyMeta | null) =>
  [...new Set((m?.groups ?? []).map((g) => g.compound.trim()).filter((c) => c))].join(', ');

/** The study's details to edit: its fields, and its groups as a table (one row per group, typed or pasted into like
 *  any Datagrok grid). `save` gets them on OK, with an ID. */
export function detailsDialog(name: string, m: StudyMeta, groups: string[], save: (m: StudyMeta) => void): void {
  const meta = withId(m);
  const id = ui.input.string('ID', {value: meta.id});
  const description = ui.input.textArea('Description', {value: meta.description});
  const owner = ui.input.string('Owner', {value: meta.owner});
  const project = ui.input.string('Project', {value: meta.project});
  const tags = ui.input.string('Tags', {value: meta.tags.join(', '),
    tooltipText: 'Comma-separated, e.g. syngeneic, IO, combination'});
  const start = ui.input.date('Day 0', {tooltipText: 'Randomization, the first day of treatment'});
  if (meta.start)
    start.stringValue = meta.start;
  const rows = groupsOf(meta, groups);
  const table = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('group', rows.map((g) => g.group)),
    DG.Column.fromStrings('compound', rows.map((g) => g.compound)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'dose', rows.map((g) => g.dose)),
    DG.Column.fromStrings('unit', rows.map((g) => g.unit)),
    DG.Column.fromStrings('route', rows.map((g) => g.route)),
    DG.Column.fromStrings('schedule', rows.map((g) => g.schedule))]);
  table.col('route')!.meta.choices = ['p.o.', 'i.p.', 'i.v.', 's.c.', 'i.t.'];
  const grid = DG.Viewer.grid(table, {allowEdit: true, showAddNewRowIcon: false, showRemoveRowIcon: false});
  grid.root.style.width = '640px';
  grid.root.style.height = `${Math.min(8, rows.length + 1) * 30 + 8}px`;
  ui.dialog(`Study details: ${name}`)
    .add(ui.inputs([id, description, owner, project, tags, start]))
    .add(ui.divText('What each group receives (the schedule as written in the protocol, e.g. Q3D×4 or QD×14):'))
    .add(grid.root)
    .onOK(() => {
      const day = (start.value as unknown as {format(f: string): string} | null)?.format('YYYY-MM-DD') ?? null;
      const cell = (c: string, i: number) => table.col(c)!.isNone(i) ? '' : String(table.get(c, i)).trim();
      save({id: id.value.trim() || meta.id, description: description.value.trim(), owner: owner.value.trim(),
        project: project.value.trim(), tags: tags.value.split(',').map((t) => t.trim()).filter((t) => t),
        start: day, groups: rows.map((g, i) => ({group: g.group, compound: cell('compound', i),
          dose: table.col('dose')!.isNone(i) ? null : table.get('dose', i), unit: cell('unit', i),
          route: cell('route', i), schedule: cell('schedule', i)}))});
    }).show();
}
