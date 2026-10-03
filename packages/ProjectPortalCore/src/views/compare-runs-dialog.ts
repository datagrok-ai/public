/** Compare a criterion's candidate runs in Compute2 (F-29): the runs by model, what to add in
 * Compare Runs, and a button that opens it. Compare Runs takes no runs as input, so the portal
 * cannot seed it; the dialog stays open beside the tool as the list to work from. Every control
 * carries a data-u2-name. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divH, divV, span} from '@datagrok-libraries/u2';
import {CriterionRuns} from '../model-runs';
import {comparePlanFrame, criterionComparePlan, openCompareRuns} from '../compare-runs';
import {describePlan} from '../model/compare-runs';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

export interface CompareRunsDialogOptions {
  program: string;
  criterion: {code: string, name: string};
  runs: CriterionRuns;
  /** Falls back to the portal's own comparison (the choose-run dialog). */
  compareHere?: () => void;
}

export function openCompareRunsDialog(o: CompareRunsDialogOptions): DG.Dialog {
  installStyles();
  const plan = criterionComparePlan(o.runs);
  const dlg = ui.dialog({title: `Compare runs · ${o.criterion.name}`});
  dlg.root.dataset.u2Name = 'compare-runs-dialog';
  const groups = named(divV(plan.groups.map((g, i) => named(divV([
    span(g.modelName, 'pp-heading'),
    ...g.runs.map((r) => span(`${r.title}${r.on ? ` · ${r.on.slice(0, 10)}` : ''}${r.selected ? ' · selected' : ''}` +
      `${r.runId ? '' : ' · no saved call'}`, r.runId ? '' : 'pp-muted')),
  ], 'pp-units'), `compare-runs-model-${i}`)), 'pp-units'), 'compare-runs-groups');
  const open = named(ui.button('Open Compare Runs', () => void openCompareRuns()
    .catch((e) => grok.shell.error(e?.message ?? String(e)))), 'compare-runs-open');
  open.disabled = plan.unavailable !== null;
  open.title = plan.unavailable ?? 'Compute2\'s run comparison, in a view of its own; this list stays open';
  const table = named(ui.button('Open as table', () => grok.shell.addTableView(comparePlanFrame(plan,
    `${o.program} · ${o.criterion.code} runs`))), 'compare-runs-table');
  table.title = 'The candidates with their run ids';
  const here = o.compareHere ? named(ui.button('Compare here', () => {
    dlg.close();
    o.compareHere!();
  }), 'compare-runs-here') : null;
  if (here)
    here.title = 'The portal\'s side-by-side of inputs and outputs, for runs of one model';
  dlg.add(named(divV([
    named(span(`${o.program} · ${o.criterion.name}: ${describePlan(plan)}${plan.crossModel ? ', across models' : ''}`),
      'compare-runs-summary'),
    groups,
    ...(plan.unavailable ? [named(span(plan.unavailable, 'pp-gap'), 'compare-runs-unavailable')] :
      [named(divV(plan.steps.map((s, i) => span(`${i + 1}. ${s}`))), 'compare-runs-steps')]),
    ...(plan.notAddable ? [span(`${plan.notAddable} run(s) have no saved call; Compare Runs cannot add them`, 'pp-muted')] : []),
    named(divH([open, table, ...(here ? [here] : [])], 'pp-actions'), 'compare-runs-actions'),
  ]), 'compare-runs-content'));
  dlg.show({width: 560});
  return dlg;
}
