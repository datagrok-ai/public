/* The study's own table in Monitor and Analyze: Datagrok's grid on the measurements, to type, correct, add or paste
   them, with what the outlier rules find painted on the cells (outlier-rules.ts) and a person's marks on them. A mark
   is one of two: flagged, a value to look at that stays in the analysis, or excluded, a value or a whole animal left
   out of it, always with its reason. Both live in the table itself, in two hidden columns (~flag, ~flag_reason), so
   they travel with the study file, and each is written to the study's log (who, when, what). Files open here too:
   .csv, .tsv, .txt, .xlsx (PowerPack's reader) and .pmxstudy; an .xls is asked to be saved as .xlsx. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {button, div, divH, span} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {DEFAULT_RULES, Finding, OutlierRules, RulesSource, canEditLabRules, findOutliers, labRules, rulesInForce,
  saveLabRules, saveUserRules, userRules} from './outlier-rules';
import {StudyInputs as StudyNames} from './explore';
import {ROLES, WEIGHT} from './open-study';
import {STUDY_EXT, addLog, parseStudy, tableOf} from './study-file';

export const FLAG = '~flag';
export const FLAG_REASON = '~flag_reason';
export type Mark = 'flagged' | 'excluded';

const COLORS = {finding: DG.Color.fromHtml('#fff1c2'), flagged: DG.Color.fromHtml('#ffd8a8'),
  excluded: DG.Color.fromHtml('#9aa5b1')};

/** `x` (a call's inputs) with the excluded rows of its study left out: the engine never sees them. The columns among
 *  the inputs are taken from the kept table by name. */
export function withoutExcluded<T extends object>(x: T): T {
  const study = (x as {study?: unknown}).study;
  const flag = study instanceof DG.DataFrame ? study.col(FLAG) : null;
  if (!(study instanceof DG.DataFrame) || !flag || !flag.toList().includes('excluded'))
    return x;
  const kept = study.clone(DG.BitSet.create(study.rowCount, (i) => flag.get(i) !== 'excluded'));
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

/** The marks of a table, row by row. */
function marksOf(df: DG.DataFrame): {state: Mark | '', reason: string}[] {
  const f = df.col(FLAG);
  const r = df.col(FLAG_REASON);
  return Array.from({length: df.rowCount}, (_, i) => ({state: (f?.get(i) ?? '') as Mark | '',
    reason: r?.get(i) ?? ''}));
}

/** Marks rows (or clears them, `state` null) and writes it to the study's log. */
function mark(df: DG.DataFrame, rows: number[], state: Mark | null, reason: string, what: string): void {
  const f = df.columns.getOrCreate(FLAG, DG.TYPE.STRING);
  const r = df.columns.getOrCreate(FLAG_REASON, DG.TYPE.STRING);
  for (const i of rows) {
    f.set(i, state ?? '', false);
    r.set(i, state ? reason : '', false);
  }
  df.fireValuesChanged();
  addLog(df, state ? `${state === 'excluded' ? 'Excluded' : 'Flagged'} ${what}${reason ? `: ${reason}` : ''}` :
    `Cleared the mark on ${what}`);
}

/** Asks why, then calls `done` with the reason (required to exclude, optional to flag). */
function askReason(title: string, required: boolean, done: (reason: string) => void): void {
  const reason = ui.input.textArea('Reason', {value: ''});
  ui.dialog(title).add(reason).onOK(() => {
    const text = (reason.value ?? '').trim();
    if (required && text === '')
      grok.shell.warning('Write why: an exclusion is reported with its reason');
    else
      done(text);
  }).show();
}

/** A table's rows as the rules read them. */
function columnsOf(n: StudyNames) {
  const values = (name: string | null) => {
    const c = name ? n.study.col(name) : null;
    return c ? Array.from({length: n.study.rowCount}, (_, i) => c.isNone(i) ? null : c.get(i)) : null;
  };
  const text = (name: string) => (values(name) ?? []).map((a) => a === null ? null : String(a));
  return {animal: text(n.animal), group: text(n.group), day: values(n.day) ?? [], volume: values(n.volume) ?? [],
    weight: values(n.bodyWeight)};
}

/** A new study typed by hand: its columns, no rows (the grid's last row adds one). */
export function emptyStudy(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.string('animal', 0), DG.Column.string('group', 0),
    DG.Column.float('day', 0), DG.Column.float('volume', 0), DG.Column.float('body_weight', 0)]);
  df.name = 'New study';
  return df;
}

/** A study file as a table: .csv, .tsv, .txt, .xlsx (the sheet that looks most like a study) or .pmxstudy. */
export async function tableOfFile(file: File): Promise<DG.DataFrame | null> {
  const ext = (file.name.split('.').pop() ?? '').toLowerCase();
  let df: DG.DataFrame | null;
  if (ext === 'xls') {
    grok.shell.warning(`${file.name}: save it as .xlsx (in Excel, File, Save As) and open that`);
    return null;
  } else if (ext === 'csv' || ext === 'txt')
    df = DG.DataFrame.fromCsv(await file.text());
  else if (ext === 'tsv')
    df = grok.data.parseCsv(await file.text(), {delimiter: '\t'});
  else if (ext === STUDY_EXT)
    return tableOf(parseStudy(await file.text()));
  else if (ext === 'xlsx') {
    const sheets: DG.DataFrame[] = await grok.functions.call('PowerPack:xlsxFileHandler',
      {bytes: new Uint8Array(await file.arrayBuffer())});
    const fit = (t: DG.DataFrame) => ROLES.filter(([, re]) => t.columns.names().some((c) => re.test(c))).length;
    df = [...sheets].sort((a, b) => fit(b) - fit(a))[0] ?? null;
  } else {
    grok.shell.warning(`${file.name}: open a .csv, .tsv, .txt, .xlsx or .${STUDY_EXT} file`);
    return null;
  }
  if (df)
    df.name = file.name.replace(/\.[^.]+$/, '');
  return df;
}

/** Rows pasted from a spreadsheet, added at the end of `df`: with a header row, each pasted column goes to the table's
 *  column of its name (or of its role: animal, group, day, volume, body weight); without one, in the table's order.
 *  Returns how many rows were added. */
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

/** The rules dialog: each rule on or off with its threshold; saved as the user's own, or for the whole lab by who may,
 *  and back to the lab's. */
async function rulesDialog(changed: () => void): Promise<void> {
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
      changed();
      d.close();
    });
  }
  if (canLab) {
    d.addButton('Save for the lab', () => {
      saveLabRules(read()).then(() => {
        grok.shell.info('The lab’s rules are saved: everyone without rules of their own follows them');
        changed();
        d.close();
      }, (e) => grok.shell.error(`The lab’s rules were not saved: ${e?.message ?? e}`));
    });
  } else if (labRules())
    d.add(ui.divText('The lab’s rules are changed by its responsible group.'));
  d.onOK(() => {
    saveUserRules(read());
    changed();
  });
  d.show();
}

export interface StudyTableHost {
  /** the study table and its columns as the form holds them (null: none yet) */
  names: () => StudyNames | null;
  /** fires when the form's study may have changed (each result) */
  result: ReadonlySignal<unknown>;
  /** makes a table the view's study (an imported file, a new table) */
  open: (df: DG.DataFrame) => void;
  /** runs `fn` again whenever a signal it reads changes, for as long as the view lives */
  effect: (fn: () => void) => void;
}

/** The panel: its bar (what the rules found and the marks, and the actions) over the study's grid. */
export function studyTablePanel(host: StudyTableHost): {content: HTMLElement[], dispose: () => void} {
  const summary = span('');
  const fileInput = document.createElement('input');
  fileInput.type = 'file';
  fileInput.accept = `.csv,.tsv,.txt,.xlsx,.xls,.${STUDY_EXT}`;
  fileInput.hidden = true;
  fileInput.addEventListener('change', () => {
    const f = fileInput.files?.[0];
    fileInput.value = '';
    if (f) {
      tableOfFile(f).then((df) => df && host.open(df),
        (e) => grok.shell.error(`${f.name} could not be read: ${e?.message ?? e}`));
    }
  });
  const gridBox = div([], 'pmx-study-grid');

  let grid: DG.Grid | null = null;
  let subs: {unsubscribe(): void}[] = [];
  let findings = new Map<number, Finding[]>();
  let current: StudyNames | null = null;
  let refresh = 0;

  const recount = () => {
    if (!current) {
      summary.textContent = 'Open the study’s table, import a file or start a new one.';
      return;
    }
    const marks = marksOf(current.study);
    const flagged = marks.filter((m) => m.state === 'flagged').length;
    const excluded = marks.filter((m) => m.state === 'excluded').length;
    const whose = {user: 'your rules', lab: 'the lab’s rules', default: 'the default rules'}[rulesInForce().from];
    summary.textContent = `${findings.size} value${findings.size === 1 ? '' : 's'} to look at (${whose}) · ` +
      `${flagged} flagged · ${excluded} excluded from the analysis`;
  };
  const rerunRules = () => {
    if (!current)
      return;
    findings = new Map();
    for (const f of findOutliers(columnsOf(current), rulesInForce().rules))
      findings.set(f.row, [...(findings.get(f.row) ?? []), f]);
    recount();
    grid?.invalidate();
  };
  const later = () => {
    window.clearTimeout(refresh);
    refresh = window.setTimeout(rerunRules, 300);
  };

  // what a cell is about: the volume, the body weight, or the row (animal, group, day)
  const kindOf = (column: string): Finding['on'] | null => !current ? null : column === current.volume ? 'volume' :
    column === current.bodyWeight ? 'weight' : [current.animal, current.group, current.day].includes(column) ?
      'row' : null;

  const menuAt = (n: StudyNames, i: number) => {
    const df = n.study;
    const animal = String(df.get(n.animal, i));
    const selected = df.selection.get(i) && df.selection.trueCount > 1 ?
      Array.from(df.selection.getSelectedIndexes()) : null;
    const rowsOf = (a: string) => Array.from({length: df.rowCount}, (_, k) => k)
      .filter((k) => String(df.get(n.animal, k)) === a);
    const one = selected ? `${selected.length} selected measurements` : `${animal} on day ${df.get(n.day, i)}`;
    const these = selected ?? [i];
    const m = DG.Menu.popup();
    m.item(`Flag ${one}…`, () => askReason(`Flag ${one}`, false, (r) => mark(df, these, 'flagged', r, one)));
    m.item(`Exclude ${one}…`, () => askReason(`Exclude ${one}`, true,
      (r) => mark(df, these, 'excluded', r, one)));
    m.separator();
    m.item(`Flag animal ${animal}…`, () => askReason(`Flag animal ${animal}`, false,
      (r) => mark(df, rowsOf(animal), 'flagged', r, `animal ${animal}`)));
    m.item(`Exclude animal ${animal}…`, () => askReason(`Exclude animal ${animal}`, true,
      (r) => mark(df, rowsOf(animal), 'excluded', r, `animal ${animal}`)));
    const marks = marksOf(df);
    if (these.some((k) => marks[k].state !== ''))
      m.item(`Clear the mark on ${one}`, () => mark(df, these, null, '', one));
    if (rowsOf(animal).some((k) => marks[k].state !== ''))
      m.item(`Clear the marks on animal ${animal}`, () => mark(df, rowsOf(animal), null, '', `animal ${animal}`));
    m.show();
  };

  const stop = () => {
    subs.forEach((s) => s.unsubscribe());
    subs = [];
    grid?.detach();
    grid = null;
  };
  const show = (now: StudyNames | null) => {
    // the same table and roles: the grid stays, only what the rules find may change
    if (grid && now && current && now.study === current.study) {
      const roles = (n: StudyNames) => [n.animal, n.group, n.day, n.volume, n.bodyWeight].join('\u0000');
      const same = roles(now) === roles(current);
      current = now;
      if (!same)
        rerunRules();
      return;
    }
    stop();
    current = now;
    gridBox.replaceChildren();
    if (!now) {
      recount();
      return;
    }
    const df = now.study;
    const g = DG.Viewer.grid(df, {allowEdit: true, showAddNewRowIcon: true, showRemoveRowIcon: true});
    grid = g;
    g.root.style.width = '100%';
    g.root.style.height = '100%';
    gridBox.append(g.root);
    subs.push(g.onCellPrepare((gc) => {
      if (!gc.isTableCell || gc.tableRowIndex === null || !gc.tableColumn)
        return;
      const i = gc.tableRowIndex;
      const state = df.col(FLAG)?.get(i);
      if (state === 'excluded')
        gc.style.textColor = COLORS.excluded;
      const kind = kindOf(gc.tableColumn.name);
      if (state === 'flagged')
        gc.style.backColor = COLORS.flagged;
      else if (kind && findings.get(i)?.some((f) => f.on === kind || f.on === 'row'))
        gc.style.backColor = COLORS.finding;
    }));
    g.onCellTooltip((gc, x, y) => {
      if (!gc.isTableCell || gc.tableRowIndex === null)
        return false;
      const i = gc.tableRowIndex;
      const m = marksOf(df)[i];
      const lines = [...(m.state ? [`${m.state === 'excluded' ? 'Excluded from the analysis' : 'Flagged'}` +
        `${m.reason ? `: ${m.reason}` : ''}`] : []), ...(findings.get(i) ?? []).map((f) => f.reason)];
      if (lines.length === 0)
        return false;
      ui.tooltip.show(ui.divV(lines.map((l) => ui.divText(l))), x, y);
      return true;
    });
    const onMenu = (e: MouseEvent) => {
      const box = g.root.getBoundingClientRect();
      const cell = g.hitTest(e.clientX - box.left, e.clientY - box.top);
      if (!current || !cell?.isTableCell || cell.tableRowIndex === null)
        return;
      e.preventDefault();
      e.stopImmediatePropagation();
      menuAt(current, cell.tableRowIndex);
    };
    g.root.addEventListener('contextmenu', onMenu, true);
    subs.push({unsubscribe: () => g.root.removeEventListener('contextmenu', onMenu, true)});
    subs.push(df.onDataChanged.subscribe(later));
    rerunRules();
  };
  host.effect(() => {
    void host.result.value;
    show(host.names());
  });

  const paste = () => {
    const df = current?.study;
    if (!df) {
      grok.shell.warning('Open the study’s table or start a new one first');
      return;
    }
    navigator.clipboard.readText().then((text) => {
      const n = appendPasted(df, text);
      grok.shell.info(n ? `${n} row${n === 1 ? '' : 's'} added from the clipboard` : 'Nothing to paste');
    }, () => grok.shell.warning('The browser did not give the clipboard: allow it, or paste into the grid'));
  };
  const openInDatagrok = () => {
    const df = current?.study;
    if (!df)
      return;
    const open = Array.from(grok.shell.tableViews).find((v) => v.dataFrame === df);
    if (open)
      grok.shell.v = open;
    else
      grok.shell.addTableView(df);
  };
  const bar = divH([summary,
    divH([button('Import file…', () => fileInput.click()), button('Paste rows', paste),
      button('New table', () => host.open(emptyStudy())), button('Rules…', () => void rulesDialog(rerunRules)),
      button('Open in Datagrok', openInDatagrok)], 'pmx-study-actions')], 'pmx-study-bar');
  return {content: [bar, fileInput, gridBox], dispose: () => {
    window.clearTimeout(refresh);
    stop();
  }};
}
