import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {
  signal, computed, Component, Scope, VirtualList, TextInput,
  Splitter, divV, divH, span, button, dot, loader, tableFromMap, iconButton,
} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import {dapiPager, chip} from '@datagrok-libraries/u2/src/dg/index.js';

const PAGE_SIZE = 30;
const ROW_HEIGHT = 44;

/** The view content plus the chrome that belongs to the shell (see the shell-integration
 * recipe): ribbon groups for the main view controls, a signal for the per-view status bar. */
export interface ReportsBrowserParts {
  content: Component;
  ribbon: HTMLElement[][];
  status: ReadonlySignal<string>;
}

function created(date: dayjs.Dayjs): string {
  return date && date.isValid() ? date.format('MMM D, YYYY HH:mm') : '';
}

function reportRow(r: DG.UserReport): HTMLElement {
  return divV([
    divH([span(r.description ?? '', 'u2demo-report-text'),
      span(r.reporter?.friendlyName ?? '', 'u2demo-dim')], 'u2demo-report-line'),
    divH([span(created(r.createdOn), 'u2demo-dim'),
      dot(r.isResolved ? 'success' : 'warning')], 'u2demo-report-line'),
  ], 'u2demo-report');
}

function reportDetails(r: DG.UserReport): HTMLElement {
  return tableFromMap({
    'Description': r.description ?? '',
    'Created': created(r.createdOn),
    'Jira ticket': r.jiraTicket ?? '',
    'Resolved': String(r.isResolved),
    'Reporter': r.reporter ? chip(r.reporter).root : '',
    'Assignee': r.assignee ? chip(r.assignee).root : '',
  });
}

export function buildReportsBrowser(): ReportsBrowserParts {
  let ribbon!: HTMLElement[][];
  let status!: ReadonlySignal<string>;
  const content = Component.build(() => {
    const scope = Scope.ambient!;
    const query = signal('');
    // reporter/assignee are stubs unless asked for; description is not a filterable field on
    // this endpoint (probed 2026-08-14), so the box passes the smart filter through verbatim
    const pager = dapiPager<DG.UserReport>(() => grok.dapi.reports.include('reporter,assignee'), {
      pageSize: PAGE_SIZE, order: 'createdOn', desc: true,
      filter: () => query.value.trim() || undefined,
    });
    scope.own(() => pager.dispose());

    const list = new VirtualList<DG.UserReport>({itemHeight: ROW_HEIGHT, keyOf: (r) => r.id, render: reportRow});
    list.setItems(pager.items);
    const onScroll = () => {
      const root = list.root;
      if (root.scrollTop + root.clientHeight > root.scrollHeight - 5 * ROW_HEIGHT)
        pager.loadMore();
    };
    list.root.addEventListener('scroll', onScroll);
    list.own(() => list.root.removeEventListener('scroll', onScroll));

    const stateArea = divV([], 'u2demo-report-status');
    const retry = divH([span('Could not load the reports.'), button('Retry', () => pager.loadMore())], 'u2demo-row');
    scope.effect(() => {
      const state = pager.state.value;
      if (state === 'loading')
        stateArea.replaceChildren(loader('Loading reports…'));
      else if (state === 'error')
        stateArea.replaceChildren(retry);
      else if (state === 'done' && pager.items.value.length === 0)
        stateArea.replaceChildren(span('No reports match the filter.', 'u2demo-hint'));
      else
        stateArea.replaceChildren();
    });

    // the chips in the details are platform components: each rendering gets its own scope, so the
    // previous report's chips are disposed before the next ones are built
    const current = computed(() => pager.items.value[list.selectedIndex.value] ?? null);
    const details = divV([], 'u2demo-report-details');
    let shown: Scope | undefined;
    scope.own(() => shown?.dispose());
    scope.effect(() => {
      const r = current.value;
      shown?.dispose();
      shown = new Scope();
      details.replaceChildren(Scope.runWith(shown, () =>
        r ? reportDetails(r) : span('Select a report to see its details.', 'u2demo-hint')));
      if (r)
        grok.shell.o = r;
    });

    const search = new TextInput({search: true, inline: true, placeholder: 'Filter, e.g. isResolved = false',
      bind: query, onChanged: () => pager.reset()});
    ribbon = [[search.root], [iconButton('sync', () => pager.reset(), {tooltip: 'Reload from the server'})]];
    status = computed(() => {
      const total = pager.total.value ?? '…';
      return pager.state.value === 'loading' ? 'Loading reports…'
        : `${pager.items.value.length} of ${total} reports`;
    });

    const split = new Splitter([divV([list, stateArea], 'u2demo-report-list'), details],
      {direction: 'horizontal', sizes: [40, 60]});
    split.root.classList.add('u2demo-report-split');

    pager.reset();
    return divV([split], 'u2demo-root');
  });
  return {content, ribbon, status};
}
