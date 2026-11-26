// the patients' exposure, typed by drug (Analyze, a group's details, Compare): a description beside each group's
// exposure, kept with the study and never analyzed
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {StudyMeta, armParts} from './study-meta';
import {HUMAN_UNITS, HumanExposure, ReferenceRow, STAT_TEXT, cvOk, fromReference, referenceKept, referenceLine,
  referenceRow}
  from './exposure';
import {humanReference} from './priors';
import {closeNamed, dialogLine, okNamed, showCentred} from './ui';
import {esc} from './sheet';

// (what is optional is said on its header's hover: a wrapped header left "( optional)" on a line of its own)
const COLS = {drug: 'Drug', auc: 'AUC0–24 at steady state (unit·h)', cmax: 'Cmax at steady state',
  unit: 'Concentration unit', fuMouse: 'Unbound fraction, mouse', fuHuman: 'Unbound fraction, patients',
  source: 'Measured or projected', cv: 'Patients\' CV % of the AUC'};
const SOURCES = ['projected', 'measured'];

/** Each drug of the groups once (each group's details first, then its name), by its key, as first written (the
 *  drugs of a combination too: their own groups may be in the study). */
export function drugsOfArms(arms: string[], meta: StudyMeta | null = null): Map<string, string> {
  const out = new Map<string, string>();
  for (const p of arms.flatMap((a) => armParts(a, meta))) {
    if (!out.has(p.key) && p.key !== 'radiation')
      out.set(p.key, p.drug);
  }
  return out;
}

/** What the field says once the human exposure is typed: "Cpd-A (projected)", or empty. */
export function humanExposureSummary(m: Map<string, HumanExposure>): string {
  return [...m.values()].filter((h) => h.auc !== null).map((h) => `${h.drug} (${h.source})`).join(', ');
}

/** PMX Studio's reference table as a Datagrok table, to read where each patients' value comes from. */
export function openReferenceTable(rows: ReferenceRow[]): void {
  const t = DG.DataFrame.fromObjects(rows.map((r) => ({
    'Drug': r.drug, 'Also called': r.synonyms.join(', '), 'Regimen': r.regimen, 'Measured as': r.analyte ?? '',
    'AUC (ng·h/mL)': r.auc, 'AUC covers': r.aucWindow, 'AUC CV %': r.aucCv, 'Cmax (ng/mL)': r.cmax,
    'Cmax CV %': r.cmaxCv, 'Unbound fraction': r.fu, 'Matrix': r.matrix.replace('_', ' '),
    'Statistic': STAT_TEXT[r.stat] ?? '', 'Source': r.source})))!;
  t.name = 'Patients\' exposure (PMX Studio reference)';
  grok.shell.addTableView(t);
}

/** The human exposure of each drug of `arms`, in an editable table (typed or pasted into like any Datagrok grid),
 *  started from `now`, else from PMX Studio's reference table when the drug is in it; `done` gets it, by drug key, on
 *  OK. */
export async function humanExposureDialog(arms: string[], now: Map<string, HumanExposure>,
  done: (m: Map<string, HumanExposure>) => void, meta: StudyMeta | null = null): Promise<void> {
  const drugs = [...drugsOfArms(arms, meta)];
  if (drugs.length === 0) {
    showCentred(closeNamed(ui.dialog('Patients\' exposure')).add(ui.divText('The groups name no drug.')));
    return;
  }
  const rows = await humanReference().then((t) => t.drugs).catch(() => [] as ReferenceRow[]);
  const inTable = drugs.map(([k]) => referenceRow(k, rows));
  const refs = drugs.map(([k], i) => now.has(k) ? now.get(k)!.reference?.row : inTable[i]);
  const start = drugs.map(([k, d], i) => now.get(k) ?? (refs[i] ? fromReference(d, refs[i]!) : undefined));
  const table = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(COLS.drug, drugs.map(([, d]) => d)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, COLS.auc, start.map((h) => h?.auc ?? null)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, COLS.cmax, start.map((h) => h?.cmax ?? null)),
    DG.Column.fromStrings(COLS.unit, start.map((h) => h?.unit ?? HUMAN_UNITS[0])),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, COLS.fuMouse, start.map((h) => h?.fuMouse ?? null)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, COLS.fuHuman, start.map((h) => h?.fuHuman ?? null)),
    DG.Column.fromStrings(COLS.source, start.map((h) => h?.source ?? SOURCES[0])),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, COLS.cv, start.map((h) => h?.cv ?? null))]);
  table.col(COLS.unit)!.meta.choices = HUMAN_UNITS;
  table.col(COLS.source)!.meta.choices = SOURCES;
  table.col(COLS.drug)!.setTag('editable', 'false');
  // the detail of each column, on its header
  const tips: [string, string][] = [
    [COLS.auc, 'Daily AUC at steady state, in the concentration unit times hours (ng/mL: ng·h/mL). The mouse\'s ' +
      'AUC is shown as NCA gives it, in mg/L·h (1 mg/L·h = 1000 ng·h/mL): the two are converted to compare them'],
    [COLS.cmax, 'Optional: if known'],
    [COLS.fuMouse, 'Optional. With both unbound fractions, exposures are compared unbound'],
    [COLS.fuHuman, 'Optional. With both unbound fractions, exposures are compared unbound'],
    [COLS.source, 'Measured: a label, a clinical study. Projected: your DMPK team\'s projection'],
    [COLS.cv, 'Optional. A label\'s "mean (CV%)": marks a group outside the patients\' range; without it, ' +
      'beyond 3-fold']];
  tips.forEach(([c, tip]) => table.col(c)!.meta.description = tip);
  const grid = DG.Viewer.grid(table, {allowEdit: true, showAddNewRowIcon: false, showRemoveRowIcon: false});
  // never wider than the window (the grid scrolls its columns), and its lines wrap to it: a long source line
  // widened the dialog past the window's edge
  const width = 'min(760px, calc(100vw - 96px))';
  grid.root.style.width = width;
  grid.root.style.height = `${Math.min(8, drugs.length + 1) * 30 + 8}px`;
  const num = (c: string, i: number) => {
    const col = table.col(c)!;
    return col.isNone(i) || !Number.isFinite(col.get(i)) ? null : col.get(i) as number;
  };
  // a row as typed, its reference marked edited once a value differs from the table's
  const read = (i: number): HumanExposure => {
    const cv = num(COLS.cv, i);
    const h: HumanExposure = {drug: drugs[i][1], auc: num(COLS.auc, i), cmax: num(COLS.cmax, i),
      unit: String(table.get(COLS.unit, i) || HUMAN_UNITS[0]), fuMouse: num(COLS.fuMouse, i),
      fuHuman: num(COLS.fuHuman, i), source: table.get(COLS.source, i) === 'measured' ? 'measured' : 'projected',
      cv: cvOk(cv) ? cv : null};
    const r = refs[i];
    return r ? {...h, reference: {row: r, edited: !referenceKept(h, r)}} : h;
  };
  const lines = ui.divV([]);
  Object.assign(lines.style, {color: 'var(--grey-4)', maxWidth: width});
  // a drug typed by hand before that the table holds: its values on a click
  const useTable = (i: number) => {
    const h = fromReference(drugs[i][1], inTable[i]!);
    const values: [string, unknown][] = [[COLS.auc, h.auc], [COLS.cmax, h.cmax], [COLS.unit, h.unit],
      [COLS.fuHuman, h.fuHuman], [COLS.source, h.source], [COLS.cv, h.cv]];
    refs[i] = inTable[i];
    values.forEach(([c, v]) => table.set(c, i, v));
    showLines();
  };
  // a CV that cannot be one, said on its drug's line (it is ignored on Save)
  const badCv = (i: number) => {
    const cv = num(COLS.cv, i);
    return cv !== null && !cvOk(cv) ? [ui.divText(`⚠ ${drugs[i][1]}: CV in %, from 5 to under 300`,
      {style: {color: 'var(--orange-2, #b26a00)'}})] : [];
  };
  const showLines = () => lines.replaceChildren(...drugs.flatMap(([, drug], i) => [refs[i] ?
    ui.divText(referenceLine(read(i))) : inTable[i] ? ui.divH([ui.divText(`${drug}: in the reference table`),
      ui.button('Use table values', () => useTable(i))]) :
      ui.divText(`${drug}: not in the reference table, enter your source's values`), ...badCv(i)]));
  showLines();
  const sub = table.onDataChanged.subscribe(showLines);
  const d = ui.dialog('Patients\' exposure').add(ui.divV([
    dialogLine('Each drug\'s exposure in patients, at their regimen. Empty: skipped.', 'Measured (a label, a ' +
      'clinical study) or projected (your DMPK team\'s projection). A drug in PMX Studio\'s reference table is ' +
      'filled from it, editable. Hover a column\'s header for its detail.'),
    grid.root, lines,
    ...rows.length ? [ui.link('See PMX Studio\'s reference table', () => openReferenceTable(rows),
      'All the drugs it holds, with the regimen, statistic and source of each value')] : []]));
  d.onClose.subscribe(() => sub.unsubscribe());
  d.onOK(() => {
    const out = new Map<string, HumanExposure>();
    const badCv: string[] = [];
    drugs.forEach(([key, drug], i) => {
      const h = read(i);
      if (h.auc === null || h.auc <= 0)
        return;
      const cv = num(COLS.cv, i);
      if (cv !== null && !cvOk(cv))
        badCv.push(drug);
      out.set(key, h);
    });
    if (badCv.length) {
      // (a toast reads HTML: the drugs come from the groups' names)
      grok.shell.warning(`CV of ${esc(badCv.join(', '))}: enter the CV in %, e.g. 45 (from 5 to under 300); the ` +
        'value was ignored');
    }
    done(out);
  });
  okNamed(d, 'Save');
  showCentred(d);
}
