import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView, Filter} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {awaitCheck} from '@datagrok-libraries/utils/src/test';
import {ViewHandler} from '../view-handler';
import {getTime} from "../utils";


export class PackagesView extends UaView {
  expanded: {[key: string]: boolean} = {f: true, l: true};

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Packages';
  }

  async initViewers(): Promise<void> {
    const packagesViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Packages',
      queryName: 'PackagesUsage',
      viewerFunction: (t: DG.DataFrame) => {
        t.onSelectionChanged.subscribe(async () => {
          PackagesView.showSelectionContextPanel(t, this.uaToolbox, this.expanded, 'Packages');
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
          const dateFrom = new Date(row.time_start);
          const dateTo = new Date(row.time_end);
          cp.addPane('Details', () => {
            return ui.tableFromMap({'User': ui.render(`#{x.${row.uid}}`),
              'Package': ui.render(`#{x.${row.pid}}`),
              'From': getTime(dateFrom),
              'To': getTime(dateTo)});
          }, true);
          PackagesView.getFunctionsPane(cp, filter, [dateFrom, dateTo], [row.package], this.uaToolbox, this.expanded.f, 'Packages');
          PackagesView.getLogsPane(cp, filter);
          PackagesView.getAuditPane(cp, filter);
          grok.shell.o = cp.root;
        });
        const viewer = DG.Viewer.scatterPlot(t, {
          x: 'time_start',
          y: 'package',
          size: 'count',
          color: 'user',
          jitterSize: 5,
          markerMinSize: 10,
          markerMaxSize: 30,
          showColorSelector: false,
          showSizeSelector: false,
          showXSelector: false,
          showYSelector: false,
        });
        return viewer;
      }});

    this.viewers.push(packagesViewer);
    this.root.append(packagesViewer.root);
  }

  static showSelectionContextPanel(t: DG.DataFrame, uaToolbox: UaToolbox, expanded: {[key: string]: boolean}, backToView: string) {
    if (!t.selection.anyTrue && !t.filter.anyTrue)
      return;
    let df = t.clone(t.selection.and(t.filter));
    const gen = t.rows[Symbol.iterator]();
    const dateMin = df.getCol('time_start').stats.min;
    const dateMax = df.getCol('time_end').stats.max;
    const dateFrom = new Date(dateMin / 1000);
    const dateTo = new Date(dateMax / 1000);
    const packages: string[] = df.getCol('pid').categories.filter((c) => c !== '');
    console.log(packages);
    const users: string[] = df.getCol('uid').categories;
    df.selection.init((i) => {
      const row = gen.next().value as DG.Row;
      return dateFrom <= row.time_start && row.time_start < dateTo &&
        packages.includes(row.pid) && users.includes(row.uid);
    }, false);
    const cp = DG.Accordion.create();
    df = t.clone(t.selection);
    const groups: string[] = df.getCol('ugid').categories;
    const usersHistogram = df
      .groupBy(['uid'])
      .sum('count')
      .aggregate();

    const filter: Filter = {
      time_start: dateMin / 1000000, time_end: dateMax / 1000000,
      groups: groups, users: users, packages: packages
    };
    cp.addPane('Time interval', () => ui.tableFromMap({
      'From': getTime(dateFrom),
      'To': getTime(dateTo),
    }), true);
    const data: {[key: string]: number} = {};
    const data1: {[key: string]: HTMLElement | string} = {};
    for (const r of usersHistogram.rows) {
      const d = data[r.uid];
      data[r.uid] = r['sum(count)'] + (d ?? 0);
      data1[r.uid] = ui.render(`#{x.${r.uid}}`);
    }
    cp.addPane('Users', () => ui.table(Object.keys(data).sort((a, b) =>
      data[b] - data[a]), (k) => [data1[k], data[k]]), true);
    //cp.addPane('Users', () => ui.divV(users.map((u) => ui.render(`#{x.${u}}`))), true);
    cp.addPane('Packages', () => ui.divV(packages.map((p) => ui.render(`#{x.${p}}`))), true);
    PackagesView.getFunctionsPane(cp, filter, [dateFrom, dateTo], df.getCol('package').categories, uaToolbox, expanded.f, backToView);
    PackagesView.getLogsPane(cp, filter);
    PackagesView.getAuditPane(cp, filter);
    grok.shell.o = cp.root;
  }

  static async getFunctionsPane(cp: DG.Accordion, filter: Filter, date: Date[], packageNames: string[], uaToolbox: UaToolbox, expanded: boolean, backToView: string): Promise<void> {
    const button = ui.button('Details', async () => {
      uaToolbox.dateFromDD.value = getTime(date[0]);
      uaToolbox.dateToDD.value = getTime(date[1]);
      uaToolbox.backToView = backToView;
      uaToolbox.usersDD.value = filter.users.length === 1 ?
        (await grok.dapi.users.find(filter.users[0])).friendlyName : `${filter.users.length} users`;
      uaToolbox.packagesDD.value = packageNames.length === 1 ?
        packageNames[0] : `${packageNames.length} packages`;
      ViewHandler.getView('Functions').getScatterPlot()
        .reloadViewer({date: `${getTime(date[0], 'es-pa')}-${getTime(date[1], 'es-pa')}`,
          groups: filter.groups, packages: packageNames});
      ViewHandler.changeTab('Functions');
      uaToolbox.drilldown = ViewHandler.getCurrentView();
    });
    button.classList.add('ua-details-button');
    const fPane = cp.addPane('Functions', () => {
      return ui.wait(async () => {
        // console.log('started query');
        const df: DG.DataFrame = await grok.data.query('UsageAnalysis:PackagesContextPaneFunctions', filter);
        // console.log('ended query');
        const data: {[key: string]: number} = {};
        // console.log('started cycle');
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
        // console.log('ended cycle');
        return ui.table(Object.keys(data).sort((a, b) =>
          data[b] - data[a]), (k) => [data1[k], data[k]]);
      });
    }, expanded);
  }

  static async getLogsPane(cp: DG.Accordion, filter: Filter): Promise<void> {
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

  static async getAuditPane(cp: DG.Accordion, filter: Filter): Promise<void> {
    const pane = cp.addPane('Audit summary', () => {
      return ui.wait(async () => {
        const df = await grok.data.query('UsageAnalysis:PackagesContextPaneAudit', filter);
        const data: {[key: string]: number} = {};
        for (const r of df.rows) data[r.name] = r.count + (data[r.name] ?? 0);
        const info = pane.root.querySelector('#info') as HTMLElement;
        info.textContent = df.getCol('count').stats.sum.toString();
        info.style.removeProperty('display');
        return Object.keys(data).length ? ui.tableFromMap(data) : ui.divText('No data');
      });
    }, true);
  }
}
