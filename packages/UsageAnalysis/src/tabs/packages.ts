import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {ViewHandler} from '../view-handler';


interface Filter {
  time_start: number;
  time_end: number;
  groups: string[];
  users: string[];
  packages: string[];
}


export class PackagesView extends UaView {
  // static viewName = 'Packages';

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
          let dateFrom: Date | string = new Date(df.getCol('time_start').stats.min / 1000);
          let dateTo: Date | string = new Date(df.getCol('time_end').stats.max / 1000);
          let packages: string[] = df.getCol('pid').categories;
          t.selection.init((i) => {
            const row = gen.next().value as DG.Row;
            return dateFrom <= row.time_start && row.time_start < dateTo &&
            packages.includes(row.pid);
          }, false);
          const cp = DG.Accordion.create();
          df = t.clone(t.selection);
          const dateMin = df.getCol('time_start').stats.min;
          const dateMax = df.getCol('time_end').stats.max;
          dateFrom = new Date(dateMin / 1000);
          dateTo = new Date(dateMax / 1000);
          const groups: string[] = df.getCol('ugid').categories;
          const users: string[] = df.getCol('uid').categories;
          packages = df.getCol('pid').categories;
          const filter: Filter = {time_start: dateMin / 1000000, time_end: dateMax / 1000000,
            groups: groups, users: users, packages: packages};
          dateFrom = dateFrom.toLocaleString('es-pa', {hour12: false}).replace(',', '');
          dateTo = dateTo.toLocaleString('es-pa', {hour12: false}).replace(',', '');
          cp.addPane('Time interval', () => ui.tableFromMap({'From': dateFrom, 'To': dateTo}), true);
          cp.addPane('Users', () => ui.divV(users.map((u) => ui.render(`#{x.${u}}`))), true);
          cp.addPane('Packages', () => ui.divV(packages.map((p) => ui.render(`#{x.${p}}`))), true);
          this.getFunctionsPane(cp, filter, [dateFrom, dateTo], df.getCol('package').categories);
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
    button.classList.add('ua-button');
    const fPane = cp.addPane('Functions', () => {
      return ui.wait(async () => {
        const packageFunctions = [];
        const df: DG.DataFrame = await grok.data.query('UsageAnalysis:PackagesContextPaneFunctions', filter);
        for (const p of packageNames)
          packageFunctions.push(...await grok.dapi.functions.filter(`package.name = "${p}"`).list());
        const data: {[key: string]: {element: Element, count: number}} = {};
        for (const r of df.rows) {
          const f = packageFunctions.find((f) => f.id === r.id);
          const d = data[r.package + r.name];
          data[r.package + r.name] = {element: f ? ui.render(f) :
            (d?.element ?? ui.divText(r.name)), count: r.count + (d?.count ?? 0)};
        }
        const info = fPane.root.querySelector('#info') as HTMLElement;
        info.textContent = df.getCol('count').stats.sum.toString();
        info.style.removeProperty('display');
        info.after(button);
        return Object.keys(data).length ? ui.table(Object.keys(data).sort((a, b) =>
          data[b].count - data[a].count), (k) => [data[k].element, data[k].count]) :
          ui.divText('No data');
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
