/** The choose-run dialog (F-15): the program's runs a criterion admits, two or more of them side
 * by side (inputs and outputs), and the one chosen with a rationale. Selection instead of
 * approval: a run chosen before is superseded, not erased. Every control carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import {badge, computed, divV, divH, span} from '@datagrok-libraries/u2';
import {CriterionRuns, runComparison, SavedRunItem, selectRunForCriterion} from '../model-runs';
import {describeRunFilter, describeSelection, RunFilter} from '../model/run-selection';
import {area, check, choice, choiceValue, ErrorRegion, named, openDialog, PortalDialog} from './u2-controls';

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

export async function openChooseRunDialog(opts: ChooseRunDialogOptions): Promise<PortalDialog> {
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

  const history = runs.superseded.length ? [span('Superseded', 'pp-heading'),
    named(divV(runs.superseded.map((e) => span(describeSelection(e), 'pp-history-row'))), 'run-selection-history')] :
    [];
  const errors = new ErrorRegion('choose-run-errors', 'choose-run-error');
  let chosen: SavedRunItem;

  return openDialog({title: `Choose a run: ${criterion.name}`, name: 'choose-run-dialog', width: 620, errors}, () => {
    const rows = items.map((i, k) => {
      const box = check(`${names[k]}-box`, undefined, false, {onChanged: (on) => {
        if (on)
          checked.add(k);
        else
          checked.delete(k);
        void redraw();
      }});
      const when = (i.run.run_started ?? i.run.created_on)?.format?.('YYYY-MM-DD') ?? '';
      const label = named(document.createElement('label'), names[k]);
      label.className = 'pp-candidate';
      label.append(box.root, divV([
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
    const pick = choice('selected-run-input', 'Select', items.map((i, k) => ({code: String(k), name: titleOf(i)})), '',
      {blank: 'Choose…'});
    const rationale = area('selection-rationale-input', undefined,
      {placeholder: 'Why this run over the others', rows: 2});
    const problems = computed(() => [...runs.refused, ...(choiceValue(pick) === '' ? ['Choose the run to select'] : []),
      ...(rationale.value.value.trim() ? [] : ['Say why this run over the others'])]);
    rationale.effect(() => errors.root.replaceChildren(...problems.value.map((p) => span(p, 'pp-muted'))));
    return {
      content: named(divV([
        filterLine,
        ...(runs.selected ? [named(span(`Selected now: ${describeSelection(runs.selected)}`), 'current-selection')] :
          []),
        items.length ? named(divV(rows, 'pp-runs'), 'run-candidates') :
          named(span('No run saved to this program matches. Save one from the program page or after a run.',
            'pp-muted'),
          'no-candidates'),
        span('Compare', 'pp-heading'),
        comparison,
        pick.root,
        rationale.root,
        errors.root,
        ...history,
      ], 'pp-choose-run pp-save-run'), 'choose-run-form'),
      buttons: [{label: 'SELECT', name: 'select-run-button', primary: true,
        enabled: computed(() => problems.value.length === 0),
        run: async () => {
          chosen = items[Number(choiceValue(pick))];
          await selectRunForCriterion(opts.programId, criterion.id, chosen.run.id, rationale.value.value.trim());
        },
        then: async () => {
          grok.shell.info(`${titleOf(chosen)} selected for ${criterion.name}`);
          await opts.onSelected?.();
        }}],
    };
  });
}
