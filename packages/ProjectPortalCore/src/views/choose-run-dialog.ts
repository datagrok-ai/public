/** The choose-run dialog (F-15): the program's runs a criterion admits, two or more of them side
 * by side (inputs and outputs), and the one chosen with a rationale. Selection instead of
 * approval: a run chosen before is superseded, not erased. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {badge, divV, divH, span} from '@datagrok-libraries/u2';
import {CriterionRuns, runComparison, SavedRunItem, selectRunForCriterion} from '../model-runs';
import {describeRunFilter, describeSelection, RunFilter} from '../model/run-selection';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const slug = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');

export interface ChooseRunDialogOptions {
  programId: string;
  criterion: {id: string, code: string, name: string, runFilter?: RunFilter | null};
  runs: CriterionRuns;
  onSelected?: () => void | Promise<void>;
}

const titleOf = (i: SavedRunItem) => i.run.title ?? i.run.run_id ?? i.run.model_name;

/** Automation names from run titles, made unique: two runs may share a title. */
function candidateNames(items: SavedRunItem[]): string[] {
  const seen = new Map<string, number>();
  return items.map((i) => {
    const base = `candidate-${slug(titleOf(i))}`;
    const n = (seen.get(base) ?? 0) + 1;
    seen.set(base, n);
    return n === 1 ? base : `${base}-${n}`;
  });
}

async function comparisonTable(items: SavedRunItem[]): Promise<HTMLElement> {
  const rows = await runComparison(items);
  const table = document.createElement('table');
  table.className = 'pp-card-table pp-compare';
  const head = table.insertRow();
  for (const h of ['', ...items.map(titleOf)])
    head.append(Object.assign(document.createElement('th'), {textContent: h}));
  let section = '';
  for (const r of rows) {
    if (r.section !== section) {
      section = r.section;
      const sub = table.insertRow();
      const cell = Object.assign(document.createElement('th'), {textContent: section === 'input' ? 'Inputs' : 'Outputs'});
      cell.colSpan = items.length + 1;
      sub.append(cell);
    }
    const tr = table.insertRow();
    if (r.differs)
      tr.className = 'pp-differs';
    tr.dataset.u2Name = `compare-${slug(r.name)}`;
    tr.insertCell().textContent = r.label;
    for (const v of r.values)
      tr.insertCell().textContent = v;
  }
  return table;
}

export async function openChooseRunDialog(opts: ChooseRunDialogOptions): Promise<DG.Dialog> {
  installStyles();
  const {runs, criterion} = opts;
  const items = runs.candidates;
  const names = candidateNames(items);
  const filterNames: Record<string, string> = {};
  for (const i of items) {
    if (i.unitOperationCode && i.unitOperation)
      filterNames[i.unitOperationCode] = i.unitOperation;
    if (i.runTypeCode && i.runType)
      filterNames[i.runTypeCode] = i.runType;
    if (i.model)
      filterNames[i.model.code] = i.model.name;
  }
  const filterLine = named(span(describeRunFilter(criterion.runFilter, filterNames) +
    (runs.excluded ? ` · ${runs.excluded} other run${runs.excluded === 1 ? '' : 's'} of the program do not match` : ''),
  'pp-muted'), 'run-filter');

  const checked = new Set<number>();
  const comparison = named(divV([span('Tick two or more runs to compare them', 'pp-muted')], 'pp-compare-box'),
    'run-comparison');
  let drawn = 0;
  const redraw = async () => {
    const mine = ++drawn;
    const picked = [...checked].sort((a, b) => a - b).map((k) => items[k]);
    if (picked.length < 2) {
      comparison.replaceChildren(span('Tick two or more runs to compare them', 'pp-muted'));
      return;
    }
    comparison.replaceChildren(span('Loading…', 'pp-muted'));
    const table = await comparisonTable(picked).catch((e) => span(`Cannot compare: ${e?.message ?? e}`, 'pp-reason'));
    if (mine === drawn)
      comparison.replaceChildren(table);
  };

  const rows = items.map((i, k) => {
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.onchange = () => {
      if (box.checked)
        checked.add(k);
      else
        checked.delete(k);
      void redraw();
    };
    const when = (i.run.run_started ?? i.run.created_on)?.format?.('YYYY-MM-DD') ?? '';
    const label = named(document.createElement('label'), names[k]);
    label.className = 'pp-candidate';
    label.append(box, divV([
      span(titleOf(i)),
      divH([
        span(i.model?.name ?? i.run.model_name, 'pp-muted'),
        ...(i.unitOperation ? [span(i.unitOperation, 'pp-muted')] : []),
        ...(i.runType ? [span(i.runType, 'pp-muted')] : []),
        ...(i.run.purpose ? [badge(i.run.purpose)] : []),
        ...(i.run.id === runs.selectedRunId ? [badge('selected', {variant: 'success'})] : []),
        span(when, 'pp-muted'),
      ], 'pp-run-meta'),
    ], 'pp-run'));
    return label;
  });

  const pick = document.createElement('select');
  pick.className = 'pp-select';
  pick.replaceChildren(...[{v: '', t: 'Choose…'}, ...items.map((i, k) => ({v: String(k), t: titleOf(i)}))]
    .map((o) => Object.assign(document.createElement('option'), {value: o.v, textContent: o.t})));
  const pickRow = named(divH([span('Select', 'pp-field-label'), pick], 'pp-field'), 'selected-run-input');
  const rationale = named(document.createElement('textarea'), 'selection-rationale-input');
  rationale.placeholder = 'Why this run over the others';
  rationale.rows = 2;
  const errors = named(divV([], 'pp-warnings'), 'choose-run-errors');

  const history = runs.superseded.length ? [span('Superseded', 'pp-heading'),
    named(divV(runs.superseded.map((e) => span(describeSelection(e), 'pp-history-row'))), 'run-selection-history')] : [];

  const dlg = ui.dialog({title: `Choose a run: ${criterion.name}`});
  dlg.root.dataset.u2Name = 'choose-run-dialog';
  const ok = () => dlg.getButton('SELECT');
  const refresh = () => {
    const problems = [...runs.refused, ...(pick.value === '' ? ['Choose the run to select'] : []),
      ...(rationale.value.trim() ? [] : ['Say why this run over the others'])];
    errors.replaceChildren(...problems.map((p) => span(p, 'pp-muted')));
    ok()?.classList.toggle('disabled', problems.length > 0);
    return problems;
  };
  pick.onchange = refresh;
  rationale.oninput = refresh;

  dlg.add(named(divV([
    filterLine,
    ...(runs.selected ? [named(span(`Selected now: ${describeSelection(runs.selected)}`), 'current-selection')] : []),
    items.length ? named(divV(rows, 'pp-runs'), 'run-candidates') :
      named(span('No run saved to this program matches. Save one from the program page or after a run.', 'pp-muted'),
        'no-candidates'),
    span('Compare', 'pp-heading'),
    comparison,
    pickRow,
    rationale,
    errors,
    ...history,
  ], 'pp-choose-run'), 'choose-run-form'));

  dlg.addButton('SELECT', async () => {
    if (refresh().length > 0 || ok()?.classList.contains('pp-busy'))
      return;
    ok()?.classList.add('pp-busy', 'disabled');
    try {
      const run = items[Number(pick.value)];
      await selectRunForCriterion(opts.programId, criterion.id, run.run.id, rationale.value.trim());
      dlg.close();
      grok.shell.info(`${titleOf(run)} selected for ${criterion.name}`);
      await opts.onSelected?.();
    } catch (e: any) {
      errors.replaceChildren(named(span(e?.message ?? String(e), 'pp-reason'), 'choose-run-error'));
      ok()?.classList.remove('pp-busy', 'disabled');
    }
  });
  dlg.show({width: 620});
  named(ok(), 'select-run-button');
  refresh();
  return dlg;
}
