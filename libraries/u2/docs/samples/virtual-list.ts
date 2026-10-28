import {Control, VirtualList, component, divH, divV, dot, notify, rowActions, span, timestamp} from '@datagrok-libraries/u2';
import type {Action} from '@datagrok-libraries/u2';
import './virtual-list.css';

interface Report {
  id: string;
  title: string;
  author: string;
  created: Date;
  resolved: boolean;
}

function actionsOf(r: Report): Action[] {
  return [
    {text: 'Copy title', icon: 'copy', run: () => {
      notify.info(`Copied "${r.title}"`);
    }},
    {text: r.resolved ? 'Reopen' : 'Resolve', icon: r.resolved ? 'undo' : 'check', run: () => {
      notify.info(`${r.id} updated`);
    }},
    {text: 'Copy id', run: () => {
      notify.info(`Copied ${r.id}`);
    }},
  ];
}

export function reportList(reports: Report[]): VirtualList<Report> {
  const list = new VirtualList<Report>({itemHeight: 44, keyOf: (r) => r.id, contextActions: actionsOf,
    render: (r) => divV([
      divH([span(r.title, 'report-title'), span(r.author, 'report-author u2-p2')], {gap: 'm', cls: 'report-line'}),
      divH([timestamp(r.created), divH([rowActions(actionsOf(r)), dot(r.resolved ? 'success' : 'warning')])], {gap: 'm', cls: 'report-line'}),
    ], 'report'),
  });
  list.root.classList.add('u2-adaptive', 'report-list');
  list.setItems(reports);
  return list;
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = component(() => reportList(Array.from({length: 10000}, (_, i) => ({
    id: `R-${i + 1}`, title: `Viewer ${['tooltip', 'legend', 'axis', 'filter'][i % 4]} misses a value after reload (#${i + 1})`,
    author: ['Ada Lovelace', 'Alan Turing', 'Grace Hopper'][i % 3], created: new Date(2026, 9, 26 - i % 30, 9 + i % 8),
    resolved: i % 5 === 0}))));
  host.append(demo.root);
  return demo;
}
