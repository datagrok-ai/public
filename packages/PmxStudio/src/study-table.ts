/* A study's measurements live in Datagrok's own table view, as any table does; PMX Studio adds to it what it knows
   of a study. What the outlier rules find (outlier-rules.ts) is painted on the cells, with its reasons on hover; the
   grid's context menu gets PMX Studio's items (detectors.js) to flag a value or an animal, a value to look at that
   stays in the analysis, or to exclude it, left out of the analysis with its reason; the view's menu gets the rules
   and a paste that adds rows. The marks are columns of the table itself (Excluded, Flagged, Mark reason), visible,
   filterable and editable like any other, so they travel with the study file; each one set from the menu is written
   to the study's log (who, when, what). Nothing is ever deleted: an excluded row stays, greyed. A table is a study once
   Monitor or Analyze has run on it: its columns' roles are kept on it (STUDY_TAG), which is what the menu and the
   painting look for. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {button, divH, span} from '@datagrok-libraries/u2';

import {DEFAULT_RULES, Finding, OutlierRules, RulesSource, canEditLabRules, findOutliers, labRules, rulesInForce,
  saveLabRules, saveUserRules, userRules} from './outlier-rules';
import {StudyInputs as StudyNames} from './explore';
import {ROLES, WEIGHT} from './open-study';
import {EXCLUDED, FLAGGED, MARK_REASON, addLog, ticked} from './study-file';

/** The columns' roles of a study table, as JSON: what makes a table a study (detectors.js reads it too). */
export const STUDY_TAG = 'pmx.study';
/* The marks (study-file.ts's EXCLUDED, FLAGGED, MARK_REASON) are visible columns created on first use; a tick typed
   straight into the grid counts as well, and so does a study file's own Excluded column (true, yes, 1, x). */
export {EXCLUDED, FLAGGED, MARK_REASON};
export type Mark = 'flagged' | 'excluded';

const COLORS = {finding: DG.Color.fromHtml('#fff1c2'), flagged: DG.Color.fromHtml('#ffd8a8'),
  excluded: DG.Color.fromHtml('#9aa5b1')};

/** Keeps the columns' roles on the study table (a run of Monitor or Analyze) and decorates its open table views. */
export function tagStudy(n: StudyNames | null): void {
  if (!n)
    return;
  const roles = JSON.stringify({animal: n.animal, group: n.group, day: n.day, volume: n.volume,
    bodyWeight: n.bodyWeight});
  if (n.study.getTag(STUDY_TAG) !== roles)
    n.study.setTag(STUDY_TAG, roles);
  for (const v of grok.shell.tableViews)
    decorateStudyView(v);
}

/** A table as a study: its columns by role (null: not a study, or a role's column gone). */
export function studyOf(df: DG.DataFrame | null | undefined): StudyNames | null {
  try {
    const r = JSON.parse(df?.getTag(STUDY_TAG) || 'null');
    return df && r && [r.animal, r.group, r.day, r.volume].every((c) => df.col(c)) ?
      {study: df, ...r, bodyWeight: r.bodyWeight && df.col(r.bodyWeight) ? r.bodyWeight : null} : null;
  } catch {
    return null;
  }
}

/** `x` (a call's inputs) with the excluded rows of its study left out: the engine never sees them. The columns among
 *  the inputs are taken from the kept table by name. */
export function withoutExcluded<T extends object>(x: T): T {
  const study = (x as {study?: unknown}).study;
  if (!(study instanceof DG.DataFrame))
    return x;
  const excluded = study.col(EXCLUDED);
  const out = DG.BitSet.create(study.rowCount, (i) => ticked(excluded, i));
  if (out.trueCount === 0)
    return x;
  const kept = study.clone(out.invert());
  kept.name = study.name;
  return Object.fromEntries(Object.entries(x).map(([k, v]) => [k, v === study ? kept :
    v instanceof DG.Column && kept.col(v.name) ? kept.col(v.name) : v])) as T;
}

/** What was excluded and why, for the report (null: nothing): a whole animal once, else each measurement. */
export function exclusionsOf(n: StudyNames | null): string | null {
  const df = n?.study;
  const marks = df ? marksOf(df) : [];
  if (!n || !df || !marks.some((m) => m.state === 'excluded'))
    return null;
  const byAnimal = new Map<string, number[]>();
  for (let i = 0; i < df.rowCount; i++) {
    const a = String(df.get(n.animal, i));
    byAnimal.set(a, [...(byAnimal.get(a) ?? []), i]);
  }
  const parts: string[] = [];
  byAnimal.forEach((rows, a) => {
    const out = rows.filter((i) => marks[i].state === 'excluded');
    const why = (i: number) => marks[i].reason ? ` (${marks[i].reason})` : '';
    if (out.length === rows.length && new Set(out.map((i) => marks[i].reason)).size === 1)
      parts.push(`animal ${a}, all its measurements${why(out[0])}`);
    else
      out.forEach((i) => parts.push(`${a} on day ${df.get(n.day, i)}${why(i)}`));
  });
  return `Left out by the investigator, each with its reason (who and when in the study's log): ${parts.join('; ')}`;
}

/** A row's mark (excluded over flagged) and its reason. */
function markOf(df: DG.DataFrame, i: number): {state: Mark | '', reason: string} {
  const r = df.col(MARK_REASON);
  return {state: ticked(df.col(EXCLUDED), i) ? 'excluded' : ticked(df.col(FLAGGED), i) ? 'flagged' : '',
    reason: !r || r.isNone(i) ? '' : String(r.get(i))};
}

/** The marks of a table, row by row. */
function marksOf(df: DG.DataFrame): {state: Mark | '', reason: string}[] {
  return Array.from({length: df.rowCount}, (_, i) => markOf(df, i));
}

/** Marks rows (or clears them, `state` null) and writes it to the study's log. Flagging leaves an exclusion and its
 *  reason as they are: only clearing brings a row back. */
export function mark(df: DG.DataFrame, rows: number[], state: Mark | null, reason: string, what: string): void {
  // the three columns together, in this order, the first time a row is marked
  const [e, f, r] = state === null ? [df.col(EXCLUDED), df.col(FLAGGED), df.col(MARK_REASON)] :
    [df.columns.getOrCreate(EXCLUDED, DG.TYPE.BOOL), df.columns.getOrCreate(FLAGGED, DG.TYPE.BOOL),
      df.columns.getOrCreate(MARK_REASON, DG.TYPE.STRING)];
  for (const i of rows) {
    const wasOut = ticked(e, i);
    if (state === null) {
      e?.set(i, false, false);
      f?.set(i, false, false);
      r?.set(i, '', false);
    } else {
      (state === 'excluded' ? e : f)!.set(i, true, false);
      if (reason && !(state === 'flagged' && wasOut))
        r!.set(i, reason, false);
    }
  }
  df.fireValuesChanged();
  addLog(df, state ? `${state === 'excluded' ? 'Excluded' : 'Flagged'} ${what}${reason ? `: ${reason}` : ''}` :
    `Cleared the mark on ${what}`);
}

/** A context menu's item (detectors.js): flags or excludes the rows, asking why (required to exclude), or clears
 *  their marks. */
export function markRows(df: DG.DataFrame, rows: number[], action: string, what: string): void {
  if (action === 'clear') {
    mark(df, rows, null, '', what);
    return;
  }
  const exclude = action === 'exclude';
  const reason = ui.input.textArea('Reason', {value: ''});
  const d = ui.dialog(`${exclude ? 'Exclude' : 'Flag'} ${what}`)
    .add(ui.divText(exclude ? 'Left out of the analysis and the forecasts, and reported with its reason.' :
      'Marked to look at: it stays in the analysis.'))
    .add(reason)
    .onOK(() => {
      const text = (reason.value ?? '').trim();
      if (exclude && text === '')
        grok.shell.warning('Write why: an exclusion is reported with its reason');
      else
        mark(df, rows, exclude ? 'excluded' : 'flagged', text, what);
    }).show();
  // an exclusion is reported with its reason: OK waits for one (the dialog would close on OK, losing what was typed)
  if (exclude) {
    const ok = d.getButton('OK');
    const sync = () => ok.disabled = (reason.value ?? '').trim() === '';
    reason.input.addEventListener('input', sync);
    sync();
  }
}

/** What the rules in force find in a study, row by row. */
export function findingsOf(n: StudyNames): Map<number, Finding[]> {
  const values = (name: string | null) => {
    const c = name ? n.study.col(name) : null;
    return c ? Array.from({length: n.study.rowCount}, (_, i) => c.isNone(i) ? null : c.get(i)) : null;
  };
  const text = (name: string) => (values(name) ?? []).map((a) => a === null ? null : String(a));
  const found = new Map<number, Finding[]>();
  for (const f of findOutliers({animal: text(n.animal), group: text(n.group), day: values(n.day) ?? [],
    volume: values(n.volume) ?? [], weight: values(n.bodyWeight)}, rulesInForce().rules))
    found.set(f.row, [...(found.get(f.row) ?? []), f]);
  return found;
}

/** A new study typed by hand: its columns, no rows. */
export function emptyStudy(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.string('animal', 0), DG.Column.string('group', 0),
    DG.Column.float('day', 0), DG.Column.float('volume', 0), DG.Column.float('body_weight', 0)]);
  df.name = 'New study';
  return df;
}

/** Rows pasted from a spreadsheet, added at the end of `df` (the grid's own paste writes over the cells there):
 *  with a header row, each pasted column goes to the table's column of its name (or of its role: animal, group, day,
 *  volume, body weight); without one, in the table's order. Returns how many rows were added. */
export function appendPasted(df: DG.DataFrame, text: string): number {
  const lines = text.replace(/\r/g, '').split('\n').filter((l) => l.trim() !== '').map((l) => l.split('\t'));
  if (lines.length === 0)
    return 0;
  const visible = df.columns.names().filter((c) => !c.startsWith('~'));
  const roles: RegExp[] = [...ROLES.map(([, re]) => re), WEIGHT];
  const target = (cell: string): string | null => {
    const name = cell.trim();
    const exact = visible.find((c) => c.toLowerCase() === name.toLowerCase());
    const role = roles.find((re) => re.test(name));
    return exact ?? (role ? visible.find((c) => role.test(c)) ?? null : null);
  };
  const header = lines[0].map(target);
  const hasHeader = header.filter((h) => h !== null).length >= Math.min(2, lines[0].length) &&
    lines[0].every((c) => c.trim() === '' || Number.isNaN(Number(c.replace(',', '.'))));
  const columns = hasHeader ? header : visible.slice(0, lines[0].length);
  const body = hasHeader ? lines.slice(1) : lines;
  for (const cells of body) {
    const row = df.rows.addNew(null, false).idx;
    columns.forEach((c, k) => {
      const v = (cells[k] ?? '').trim();
      const col = c ? df.col(c) : null;
      if (col && v !== '')
        col.set(row, col.isNumerical ? Number(v.replace(',', '.')) : v, false);
    });
  }
  df.fireValuesChanged();
  return body.length;
}

// whoever shows what the rules find, told when the rules in force change
const rulesWatchers = new Set<() => void>();
function rulesChanged(): void {
  rulesWatchers.forEach((fn) => fn());
}
export function onRulesChanged(fn: () => void): () => void {
  rulesWatchers.add(fn);
  return () => rulesWatchers.delete(fn);
}

/** The rules dialog: each rule on or off with its threshold; saved as the user's own, or for the whole lab by who may,
 *  and back to the lab's. */
export async function rulesDialog(): Promise<void> {
  const {rules, from} = rulesInForce();
  const canLab = await canEditLabRules().catch(() => false);
  const impossible = ui.input.bool('Impossible values', {value: rules.impossible,
    tooltipText: 'A negative volume or body weight; an animal measured twice on one day'});
  const jumpOn = ui.input.bool('Jumps between measurements', {value: rules.jump.on});
  const fold = ui.input.float('Times larger or smaller', {value: rules.jump.fold, min: 1.5, max: 20,
    tooltipText: 'A volume this many times its previous one (a caliper typo, a misread animal)'});
  const groupOn = ui.input.bool('Far from its group', {value: rules.group.on});
  const mads = ui.input.float('Robust deviations', {value: rules.group.mads, min: 2, max: 10,
    tooltipText: 'How far from the group’s median on that day, in robust deviations (MAD) of log volume'});
  const minAnimals = ui.input.int('In groups of at least', {value: rules.group.minAnimals, min: 3, max: 50});
  const weightOn = ui.input.bool('Weight falls fast', {value: rules.weightDrop.on});
  const percent = ui.input.float('Percent since the last weighing', {value: rules.weightDrop.percent, min: 1,
    max: 50});
  const read = (): OutlierRules => ({impossible: impossible.value, jump: {on: jumpOn.value, fold: fold.value ?? 3},
    group: {on: groupOn.value, mads: mads.value ?? 3.5, minAnimals: minAnimals.value ?? 5},
    weightDrop: {on: weightOn.value, percent: percent.value ?? 10}});
  const set = (r: OutlierRules) => {
    impossible.value = r.impossible;
    jumpOn.value = r.jump.on;
    fold.value = r.jump.fold;
    groupOn.value = r.group.on;
    mads.value = r.group.mads;
    minAnimals.value = r.group.minAnimals;
    weightOn.value = r.weightDrop.on;
    percent.value = r.weightDrop.percent;
  };
  const whose: {[s in RulesSource]: string} = {user: 'Your own rules', lab: 'The lab’s rules',
    default: 'PMX Studio’s default rules'};
  const d = ui.dialog('Outlier rules');
  d.add(ui.divText(`${whose[from]} are in force. A value they find is only painted, for you to look at: nothing ` +
    'leaves the analysis unless you exclude it.'));
  d.add(ui.inputs([impossible, jumpOn, fold, groupOn, mads, minAnimals, weightOn, percent]));
  d.addButton('Defaults', () => set(DEFAULT_RULES));
  if (userRules()) {
    d.addButton('Use the lab’s', () => {
      saveUserRules(null);
      rulesChanged();
      d.close();
    });
  }
  if (canLab) {
    d.addButton('Save for the lab', () => {
      saveLabRules(read()).then(() => {
        grok.shell.info('The lab’s rules are saved: everyone without rules of their own follows them');
        rulesChanged();
        d.close();
      }, (e) => grok.shell.error(`The lab’s rules were not saved: ${e?.message ?? e}`));
    });
  } else if (labRules())
    d.add(ui.divText('The lab’s rules are changed by its responsible group.'));
  d.onOK(() => {
    saveUserRules(read());
    rulesChanged();
  });
  d.show();
}

/** A study's table view as PMX Studio shows it: what the rules find painted, the marks, the reasons on hover, and
 *  PMX Studio's menu (paste rows, the rules, back to the study's views). Once per view; a table that is not a study is
 *  left alone. */
export function decorateStudyView(view: DG.TableView): void {
  const root = view.root;
  const df = view.dataFrame;
  if (root.dataset.pmxStudy || !studyOf(df))
    return;
  root.dataset.pmxStudy = '1';
  const grid = view.grid;
  let findings = new Map<number, Finding[]>();
  let timer = 0;
  const rerun = () => {
    const n = studyOf(df);
    findings = n ? findingsOf(n) : new Map();
    grid.invalidate();
  };
  const later = () => {
    window.clearTimeout(timer);
    timer = window.setTimeout(rerun, 300);
  };
  const kindOf = (column: string): Finding['on'] | null => {
    const n = studyOf(df);
    return !n ? null : column === n.volume ? 'volume' : column === n.bodyWeight ? 'weight' :
      [n.animal, n.group, n.day].includes(column) ? 'row' : null;
  };
  const subs = [grid.onCellPrepare((gc) => {
    if (!gc.isTableCell || gc.tableRowIndex === null || !gc.tableColumn)
      return;
    const i = gc.tableRowIndex;
    const state = markOf(df, i).state;
    if (state === 'excluded')
      gc.style.textColor = COLORS.excluded;
    const kind = kindOf(gc.tableColumn.name);
    if (state === 'flagged')
      gc.style.backColor = COLORS.flagged;
    else if (kind && findings.get(i)?.some((f) => f.on === kind || f.on === 'row'))
      gc.style.backColor = COLORS.finding;
  }), df.onDataChanged.subscribe(later)];
  grid.onCellTooltip((gc, x, y) => {
    if (!gc.isTableCell || gc.tableRowIndex === null)
      return false;
    const m = markOf(df, gc.tableRowIndex);
    const lines = [...(m.state ? [`${m.state === 'excluded' ? 'Excluded from the analysis' : 'Flagged'}` +
      `${m.reason ? `: ${m.reason}` : ''}`] : []), ...(findings.get(gc.tableRowIndex) ?? []).map((f) => f.reason)];
    if (lines.length === 0)
      return false;
    ui.tooltip.show(ui.divV(lines.map((l) => ui.divText(l))), x, y);
    return true;
  });
  const stopRules = onRulesChanged(rerun);
  const removed = grok.events.onViewRemoved.subscribe((v) => {
    if (v.root !== root)
      return;
    window.clearTimeout(timer);
    subs.forEach((s) => s.unsubscribe());
    stopRules();
    removed.unsubscribe();
  });

  const paste = () => navigator.clipboard.readText().then((text) => {
    const n = appendPasted(df, text);
    grok.shell.info(n ? `${n} row${n === 1 ? '' : 's'} added from the clipboard` : 'Nothing to paste');
  }, () => grok.shell.warning('The browser did not give the clipboard: allow it, then try again'));
  const back = () => {
    const pmx = [...grok.shell.views].filter((v) => v.name.startsWith('PMX Studio ·')).pop();
    if (pmx)
      setTimeout(() => grok.shell.v = pmx);
    else
      grok.shell.warning('No PMX Studio view is open');
  };
  view.ribbonMenu.group('PMX Studio')
    .item('Paste rows at the end', () => void paste())
    .item('Outlier rules…', () => void rulesDialog())
    .item('Back to PMX Studio', back);
  rerun();
}

/** Where Monitor and Analyze show the study's data: one line, what the rules found and the marks, and the way to its
 *  table view. `refresh` after each result. */
export function studyDataLine(names: () => StudyNames | null): {root: HTMLElement, refresh: () => void,
  dispose: () => void} {
  const text = span('');
  const open = () => {
    const df = names()?.study;
    if (!df)
      return;
    const view = [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id) ?? grok.shell.addTableView(df);
    decorateStudyView(view);
    // after the click: the platform makes the clicked view current again once the click is handled (stages.ts)
    setTimeout(() => grok.shell.v = view);
  };
  const root = divH([text, button('Open the table', open)], 'pmx-study-line');
  const refresh = () => {
    const n = names();
    root.hidden = !n;
    if (!n)
      return;
    const findings = findingsOf(n).size;
    const marks = marksOf(n.study);
    const flagged = marks.filter((m) => m.state === 'flagged').length;
    const excluded = marks.filter((m) => m.state === 'excluded').length;
    const whose = {user: 'your rules', lab: 'the lab’s rules', default: 'the default rules'}[rulesInForce().from];
    text.textContent = `Study data: ${n.study.rowCount} measurements · ${findings} to look at (${whose}) · ` +
      `${flagged} flagged · ${excluded} excluded. Right-click a value in the table to flag or exclude it.`;
  };
  return {root, refresh, dispose: onRulesChanged(refresh)};
}
