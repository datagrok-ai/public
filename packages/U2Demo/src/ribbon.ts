/* The `crosshairs` toggle needs an ambient scope, so the caller wraps this in `content.run(...)`. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {dropDownButton, iconButton, notify} from '@datagrok-libraries/u2';
import {leakReport} from '@datagrok-libraries/u2/src/dg/index.js';
import type {DemoShell} from './demo';
import {pageShown} from './nav';

/** What the drop-down opened, so `Close demo tables` closes those and not the user's own work. */
const opened: DG.DataFrame[] = [];

function addTable(table: DG.DataFrame): void {
  opened.push(table);
  grok.shell.addTableView(table);
}

export function demoRibbon(shell: DemoShell): HTMLElement[] {
  const inspect = iconButton('crosshairs', () => {}, {toggle: shell.inspect,
    tooltip: 'Inspect: click any control to see its properties in the context panel'});
  inspect.classList.add('u2demo-inspect');
  return [
    iconButton('code', () => pageShown(shell.current.peek(), true),
      {tooltip: 'Show this demo\'s source in the context panel'}),
    inspect,
    iconButton('sync', () => {
      shell.rebuild();
      notify.info(`Rebuilt — ${leakReport().liveScopes} live scopes`);
    }, {tooltip: 'Rebuild this demo — disposes its scope and builds it again'}),
    dropDownButton('Demo tools', (m) => m
      .item('Add demog table', () => addTable(grok.data.demo.demog(100)))
      .item('Add molecules table', () => addTable(grok.data.demo.molecules(30)))
      .separator()
      .item('Close demo tables', () => {
        const names = new Set(grok.shell.tables.map((t) => t.name));
        for (const t of opened)
          if (names.has(t.name))
            grok.shell.closeTable(t);
        opened.length = 0;
      })
      .separator()
      .item('Leak report', () => {
        const r = leakReport();
        notify.info(`${r.liveScopes} live scopes · ${r.liveWidgets} widgets`);
      })),
  ];
}
