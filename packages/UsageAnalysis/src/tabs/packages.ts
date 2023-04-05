import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView, Filter} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {ViewHandler} from '../view-handler';


export class PackagesView extends UaView {
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox, 'package');
    this.name = 'Packages';
  }

  async initViewers(): Promise<void> {
    const packagesViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Packages',
      'PackagesUsage',
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
          t.selection.init((i) => {
            const row = gen.next().value as DG.Row;
            return dateFrom <= row.time_start && row.time_start < dateTo &&
            packages.includes(row.pid);
          }, false);
          const cp = DG.Accordion.create();
          df = t.clone(t.selection);
          const groups: string[] = df.getCol('ugid').categories;
          const users: string[] = df.getCol('uid').categories;
          const filter: Filter = {time_start: dateMin / 1000000, time_end: dateMax / 1000000,
            groups: groups, users: users, packages: packages};
          const dateFromLS = dateFrom.toLocaleString('es-pa', {hour12: false}).replace(',', '');
          const dateToLS = dateTo.toLocaleString('es-pa', {hour12: false}).replace(',', '');
          cp.addPane('Time interval', () => ui.tableFromMap({'From': dateFromLS, 'To': dateToLS}), true);
          cp.addPane('Users', () => ui.divV(users.map((u) => ui.render(`#{x.${u}}`))), true);
          cp.addPane('Packages', () => ui.divV(packages.map((p) => ui.render(`#{x.${p}}`))), true);
          this.getFunctionsPane(cp, filter, [dateFromLS, dateToLS], df.getCol('package').categories);
          this.getLogsPane(cp, filter);
          grok.shell.o = cp.root;
        });

        t.onCurrentRowChanged.subscribe(async () => {
          t.selection.setAll(false);
          const rowValues = Array.from(t.currentRow.cells).map((c) => c.value);
          const row = Object.fromEntries(t.columns.names().map((k, i) => [k, rowValues[i]]));
          row.time_start = row.time_start.a;
          row.time_end = row.time_end.a;
          const filter: Filter = {time_start: row.time_start / 1000, time_end: row.time_end / 1000,
            groups: [row.ugid], users: [row.uid], packages: [row.pid]};
          const cp = DG.Accordion.create();
          const dateFrom = new Date(row.time_start).toLocaleString('es-pa', {hour12: false}).replace(',', '');
          const dateTo = new Date(row.time_end).toLocaleString('es-pa', {hour12: false}).replace(',', '');
          cp.addPane('Details', () => {
            return ui.tableFromMap({'User': ui.render(`#{x.${row.uid}}`),
              'Package': ui.render(`#{x.${row.pid}}`),
              'From': dateFrom, 'To': dateTo});
          }, true);
          this.getFunctionsPane(cp, filter, [dateFrom, dateTo], [row.package]);
          this.getLogsPane(cp, filter);
          grok.shell.o = cp.root;
        });
      }, null, null, this.viewer);

    this.viewers.push(packagesViewer);
    this.root.append(packagesViewer.root);
  }

  async getFunctionsPane(cp: DG.Accordion, filter: Filter, date: string[], packageNames: string[]) {
    const button = ui.button('Details', async () => {
      this.uaToolbox.dateFromDD.value = date[0];
      this.uaToolbox.dateToDD.value = date[1];
      this.uaToolbox.usersDD.value = filter.users.length === 1 ?
        (await grok.dapi.users.find(filter.users[0])).friendlyName : `${filter.users.length} users`;
      this.uaToolbox.packagesDD.value = packageNames.length === 1 ?
        packageNames[0] : `${packageNames.length} packages`;
      ViewHandler.getView('Functions').getScatterPlot()
        .reloadViewer({date: `${date[0]}-${date[1]}`, groups: filter.groups, packages: packageNames});
      ViewHandler.changeTab('Functions');
      this.uaToolbox.drilldown = ViewHandler.getCurrentView();
    });
    button.classList.add('ua-details-button');
    const fPane = cp.addPane('Functions', () => {
      return ui.wait(async () => {
        console.log('started query');
        const df: DG.DataFrame = await grok.data.query('UsageAnalysis:PackagesContextPaneFunctions', filter);
        console.log('ended query');
        const data: {[key: string]: number} = {};
        console.log('started cycle');
        for (const r of df.rows) {
          const d = data[r.package + ':' + r.name];
          data[r.package + ':' + r.name] = r.count + (d ?? 0);
        }
        const info = fPane.root.querySelector('#info') as HTMLElement;
        info.textContent = df.getCol('count').stats.sum.toString();
        info.style.removeProperty('display');
        info.after(button);
        if (!Object.keys(data).length) return ui.divText('No data');
        const data1: {[key: string]: HTMLElement | string} = {};
        for (const k of Object.keys(data)) {
          const el = ui.render(`#{x.${k}}`);
          try {
            await awaitCheck(() => el.innerHTML !== '<span data-markup-ready="false"></span>', '', 200);
          } catch (e: any) {}
          data1[k] = el.innerHTML !== '<span></span>' && el.innerText !== 'error' ? el : k.split(':')[1];
        }
        console.log('ended cycle');
        return ui.table(Object.keys(data).sort((a, b) =>
          data[b] - data[a]), (k) => [data1[k], data[k]]);
      });
    }, true);
  }

  async getLogsPane(cp: DG.Accordion, filter: Filter) {
    const lPane = cp.addPane('Log events summary', () => {
      return ui.wait(async () => {
        const df = await grok.data.query('UsageAnalysis:PackagesContextPaneLogs', filter);
        const data: {[key: string]: number} = {};
        for (const r of df.rows) data[r.source] = r.count + (data[r.source] ?? 0);
        const info = lPane.root.querySelector('#info') as HTMLElement;
        info.textContent = df.getCol('count').stats.sum.toString();
        info.style.removeProperty('display');
        return Object.keys(data).length ? ui.tableFromMap(data) : ui.divText('No data');
      });
    }, true);
  }
}
