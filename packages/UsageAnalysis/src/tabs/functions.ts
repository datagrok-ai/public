import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView, Filter} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';


export class FunctionsView extends UaView {
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox, 'function', true);
    this.name = 'Functions';
  }

  async initViewers(): Promise<void> {
    const functionsViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Functions',
      'FunctionsUsage',
      (t: DG.DataFrame) => {
        t.onSelectionChanged.subscribe(async () => {
          if (!t.selection.anyTrue) return;
          let df = t.clone(t.selection);
          const gen = t.rows[Symbol.iterator]();
          const dateMin = df.getCol('time_start').stats.min;
          const dateMax = df.getCol('time_end').stats.max;
          const dateFrom = new Date(dateMin / 1000);
          const dateTo = new Date(dateMax / 1000);
          const packages: string[] = df.getCol('pid').categories;
          const functions: string[] = df.getCol('function').categories;
          t.selection.init((i) => {
            const row = gen.next().value as DG.Row;
            return dateFrom <= row.time_start && row.time_start < dateTo &&
            packages.includes(row.pid) && functions.includes(row.function);
          }, false);
          const cp = DG.Accordion.create();
          df = t.clone(t.selection);
          const users: string[] = df.getCol('uid').categories;
          const filter: Filter = {time_start: dateMin / 1000000, time_end: dateMax / 1000000,
            users: users, packages: packages, functions: functions};
          const dateFromLS = dateFrom.toLocaleString('es-pa', {hour12: false}).replace(',', '');
          const dateToLS = dateTo.toLocaleString('es-pa', {hour12: false}).replace(',', '');
          cp.addPane('Time interval', () => ui.tableFromMap({'From': dateFromLS, 'To': dateToLS}), true);
          cp.addPane('Users', () => ui.divV(users.map((u) => ui.render(`#{x.${u}}`))), true);
          // cp.addPane('Packages', () => ui.divV(packages.map((p) => ui.render(`#{x.${p}}`))), true);
          this.getFunctionPane(cp, filter);
          // this.getLogsPane(cp, filter);
          grok.shell.o = cp.root;
        });

        t.onCurrentRowChanged.subscribe(async () => {
          t.selection.setAll(false);
          const rowValues = Array.from(t.currentRow.cells).map((c) => c.value);
          const row = Object.fromEntries(t.columns.names().map((k, i) => [k, rowValues[i]]));
          row.time_start = row.time_start.a;
          row.time_end = row.time_end.a;
          const filter: Filter = {time_start: row.time_start / 1000, time_end: row.time_end / 1000,
            users: [row.uid], packages: [row.pid], functions: [row.function]};
          const cp = DG.Accordion.create();
          const dateFrom = new Date(row.time_start).toLocaleString('es-pa', {hour12: false}).replace(',', '');
          const dateTo = new Date(row.time_end).toLocaleString('es-pa', {hour12: false}).replace(',', '');
          cp.addPane('Details', () => {
            return ui.tableFromMap({'User': ui.render(`#{x.${row.uid}}`),
              'Package': ui.render(`#{x.${row.pid}}`),
              'From': dateFrom, 'To': dateTo});
          }, true);
          this.getFunctionPane(cp, filter);
          grok.shell.o = cp.root;
        });
      }, null, null, this.viewer);

    this.viewers.push(functionsViewer);
    this.root.append(functionsViewer.root);
  }

  async getFunctionPane(cp: DG.Accordion, filter: Filter) {
    const df = await grok.data.query('UsageAnalysis:FunctionsContextPane', filter);
    const data: {[key: string]: [string, any, string][]} = {};
    for (const r of df.rows) {
      const key = r.pid + ':' + r.function;
      if (!data[key]) data[key] = [];
      data[key].push([r.rid, r.time, r.run]);
    }
    const keysSorted = Object.keys(data).sort((a, b) => a[a.indexOf(':') + 1].localeCompare(b[b.indexOf(':') + 1]));
    const expand = keysSorted.length === 1;
    const packages: {[key: string]: DG.Accordion} = {};
    for (const k of keysSorted) {
      const [p, n] = k.split(':');
      if (!packages[p]) packages[p] = DG.Accordion.create();
      packages[p].addPane(n, () => {
        return ui.wait(async () => {
          // data[k].sort((a, b) => a[2].localeCompare(b[2]));
          // return ui.divV(data[k].map((i) => i[0]));
          const t = ui.table(data[k].sort((a, b) => a[2].localeCompare(b[2])),
            (i) => [ui.render(`#{x.${i[0]}."${i[2]}"}`), i[1].format('HH:mm:ss MM/DD/YYYY')]);
          t.classList.add('ua-table');
          return t;
        });
      }, expand);
    }
    Object.keys(packages).forEach((k) => {
      const pane = cp.addPane('', () => packages[k].root);
      pane.root.querySelector('.d4-accordion-pane-header')?.prepend(ui.render(`#{x.${k}}`));
    });
  }
}
