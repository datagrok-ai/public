// Opens Run Comparison views and tracks the open ones, so snapshots can be added to them.

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {filter, take} from 'rxjs/operators';
import dayjs from 'dayjs';
import type {RibbonItemBase, RibbonPanelItem} from '@datagrok-libraries/webcomponents-vue';
import {RunComparisonApp} from '../../apps/RunComparisonApp';
import {isRoleOnlyModelFilter} from '../../package-instance';
import {setVueAppOptions} from '../../utils';
import {ComparisonEntry} from './types';
import {COMPARISON_SERVICE, ComparisonService} from './comparison-service';

export interface ComparisonView {
  view: DG.ViewBase;
  service: ComparisonService;
}

export const openComparisonViews = Vue.shallowRef<ComparisonView[]>([]);

/** Opens a new Run Comparison view, loading the runs by id, adding the entries and selecting the model. */
export function openRunComparison({runIds = [], entries = [], model}: {
  runIds?: string[], entries?: ComparisonEntry[], model?: DG.Func,
} = {}) {
  const view = new DG.ViewBase();
  view.name = 'Run Comparison';
  view.helpUrl = '/help/compute/run-comparison.md';
  view.setRibbonPanels([[
    ui.iconFA('question',
      () => window.open('https://datagrok.ai/help/compute/run-comparison', '_blank'),
      'Open the Run Comparison documentation'),
  ]]);
  view.root.classList.remove('ui-panel');
  const handle: ComparisonView = {view, service: new ComparisonService()};
  const app = Vue.createApp(RunComparisonApp, {roleOnlyFilter: isRoleOnlyModelFilter()});
  app.provide(COMPARISON_SERVICE, handle.service);
  setVueAppOptions(app);
  app.mount(view.root);
  openComparisonViews.value = [...openComparisonViews.value, handle];
  handle.service.add({runIds, entries, model, selectModel: true});

  grok.events.onViewRemoved.pipe(
    filter((closedView) => closedView === view),
    take(1),
  ).subscribe(() => {
    openComparisonViews.value = openComparisonViews.value.filter((item) => item !== handle);
    app.unmount();
  });

  grok.shell.addView(view);
  return handle;
}

const NEW_COMPARISON = 'New comparison';

/** Adds the entry to the open view the user picks, the most recent one by default, or to a new view. */
export function showAddToComparisonDialog(entry: ComparisonEntry, model?: DG.Func) {
  const views = openComparisonViews.value;
  const labels = views.map((item, i) => {
    const names = item.service.entryNames;
    return `${i + 1}. ${item.view.name}: ${names.length ? names.join(', ') : 'no runs yet'}`;
  });
  const input = ui.input.radio('', {value: labels[labels.length - 1], items: [NEW_COMPARISON, ...labels]});
  ui.dialog(`Add '${entry.name}' to comparison`)
    .add(input)
    .onOK(() => {
      const target = views[labels.indexOf(input.value ?? '')];
      if (!target) {
        openRunComparison({entries: [entry], model});
        return;
      }
      target.service.add({entries: [entry], model});
      grok.shell.v = target.view;
    })
    // modal, so the listed views can't be closed while the dialog is open
    .show({center: true, modal: true});
}

/** Runs the action, asking first when the snapshot has outdated outputs. */
export function confirmOutdated(entry: ComparisonEntry, action: () => void) {
  if (!entry.outdatedSteps?.length) {
    action();
    return;
  }
  ui.dialog('Outdated results')
    .add(ui.divText(`Outputs of ${entry.outdatedSteps.join(', ')} are outdated: ` +
      'the inputs changed after the last run, or the step never ran. Compare anyway?'))
    // a centered dialog opened from onOK leaks its wrapper, so the action waits a tick
    .onOK(() => setTimeout(action))
    .show({center: true, modal: true});
}

export const snapshotName = (name: string) => `${name} (current, ${dayjs().format('HH:mm')})`;

/** Always a new view: the original run and the snapshot of its edited copy. */
export const compareWithOriginal = (snapshot: ComparisonEntry, originalRunId: string, model?: DG.Func) =>
  confirmOutdated(snapshot, () => openRunComparison({runIds: [originalRunId], entries: [snapshot], model}));

/** A new view when none is open, otherwise the user picks one. */
export const addToComparison = (snapshot: ComparisonEntry, model?: DG.Func) =>
  confirmOutdated(snapshot, () => openComparisonViews.value.length > 0 ?
    showAddToComparisonDialog(snapshot, model) :
    openRunComparison({entries: [snapshot], model}));

type ItemGuard = Partial<Omit<RibbonItemBase, 'onClick' | 'tooltip'>>;

export const compareRibbonIcon = (compare: () => void, guard: ItemGuard = {}): RibbonPanelItem => ({
  icon: 'exchange',
  tooltip: 'Add the current state to a comparison',
  onClick: compare,
  ...guard,
});
