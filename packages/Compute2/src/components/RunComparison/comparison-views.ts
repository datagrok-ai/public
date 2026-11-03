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
import {COMPARISON_SERVICE, ComparisonAddition, ComparisonService} from './comparison-service';
import {withEntryName} from './selection';

export interface ComparisonView {
  view: DG.ViewBase;
  service: ComparisonService;
}

export const openComparisonViews = Vue.shallowRef<ComparisonView[]>([]);

/** Opens a new Run Comparison view, loading the runs by id, adding the entries and selecting the model. */
export function openRunComparison({runIds = [], runNames, entries = [], model}: ComparisonAddition = {}) {
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
  handle.service.add({runIds, runNames, entries, model, selectModel: true});

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

/** A saved run to add by id, or a snapshot. A run with no name left in the dialog keeps its saved name. */
export type ComparisonItem = {runId: string, name?: string} | {entry: ComparisonEntry};

export interface AddToComparisonOptions {
  model?: DG.Func;
  // the choice starts on a new comparison, not on the most recent one
  preferNew?: boolean;
  // adds to this view, with no choice of where
  into?: ComparisonService;
  onAdded?: () => void;
}

/** Every way of comparing goes through here: each item gets an editable name, and the user picks where to add them
 * when other comparisons are open. */
export function showAddToComparisonDialog(items: ComparisonItem[],
  {model, preferNew = false, into, onAdded}: AddToComparisonOptions = {}) {
  const views = openComparisonViews.value;
  const labels = views.map((item, i) => {
    const names = item.service.entryNames;
    return `${i + 1}. ${item.view.name}: ${names.length ? names.join(', ') : 'no runs yet'}`;
  });
  const target = ui.input.radio('', {items: [NEW_COMPARISON, ...labels],
    value: preferNew || labels.length === 0 ? NEW_COMPARISON : labels[labels.length - 1]});
  const rows = items.map((item) => ({
    item, input: ui.input.string('Name', {value: 'entry' in item ? item.entry.name : item.name ?? ''}),
  }));

  const addition = (): ComparisonAddition => {
    const runIds: string[] = [];
    const runNames = new Map<string, string>();
    const entries: ComparisonEntry[] = [];
    for (const {item, input} of rows) {
      const name = input.value.trim();
      if ('entry' in item)
        entries.push(name ? withEntryName(item.entry, name) : item.entry);
      else {
        runIds.push(item.runId);
        if (name)
          runNames.set(item.runId, name);
      }
    }
    return {runIds, runNames, entries, model};
  };

  const dialog = ui.dialog('Add to comparison')
    .add(ui.divV(rows.flatMap(({item, input}) => {
      const outdated = 'entry' in item ? item.entry.outdatedSteps : undefined;
      return outdated?.length ? [input.root, ui.divText(`Outputs of ${outdated.join(', ')} are outdated: the inputs ` +
        'changed after the last run, or the step never ran.', 'c2-comparison-add-outdated')] : [input.root];
    }), 'c2-comparison-add-names'));
  if (!into && views.length > 0)
    dialog.add(target);
  dialog.onOK(() => {
    const picked = views[labels.indexOf(target.value ?? '')];
    if (into)
      into.add(addition());
    else if (picked) {
      picked.service.add(addition());
      grok.shell.v = picked.view;
    } else
      openRunComparison(addition());
    onAdded?.();
  })
    // modal, so the listed views can't be closed while the dialog is open
    .show({center: true, modal: true});
}

// seconds tell apart snapshots taken within one minute
export const snapshotName = (name: string) => `${name} (current, ${dayjs().format('HH:mm:ss')})`;

/** The original run and the snapshot of its edited copy. */
export const compareWithOriginal = (snapshot: ComparisonEntry, originalRunId: string, originalName?: string,
  model?: DG.Func) =>
  showAddToComparisonDialog([{runId: originalRunId, name: originalName}, {entry: snapshot}], {model, preferNew: true});

export const addToComparison = (snapshot: ComparisonEntry, model?: DG.Func) =>
  showAddToComparisonDialog([{entry: snapshot}], {model});

export const addRunToComparison = (runId: string, name?: string, model?: DG.Func) =>
  showAddToComparisonDialog([{runId, name}], {model});

type ItemGuard = Partial<Omit<RibbonItemBase, 'onClick' | 'tooltip'>>;

export const compareRibbonIcon = (compare: () => void, guard: ItemGuard = {}): RibbonPanelItem => ({
  icon: 'exchange',
  tooltip: 'Add the current state to a comparison',
  onClick: compare,
  ...guard,
});
