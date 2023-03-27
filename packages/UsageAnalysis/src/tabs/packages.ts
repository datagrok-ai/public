import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {ViewHandler} from '../view-handler';


export class PackagesView extends UaView {
  viewName = 'Packages';

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox, 'package');
  }

  async initViewers(): Promise<void> {
    const packagesViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Packages',
      'PackagesUsage',
      (t: DG.DataFrame) => {
        const FunctionsView = ViewHandler.getView('Functions');
        t.onCurrentRowChanged.subscribe(async () => {
          const rowValues = Array.from(t.currentRow.cells).map((c) => c.value);
          const row = Object.fromEntries(t.columns.names().map((k, i) => [k, rowValues[i]]));
          row.time_start = row.time_start.a;
          row.time_end = row.time_end.a;
          const filter = {time_start: row.time_start / 1000, time_end: row.time_end / 1000,
            user: row.uid, package: row.package};
          const cp = DG.Accordion.create();
          const dateFrom = new Date(row.time_start).toLocaleString();
          const dateTo = new Date(row.time_end).toLocaleString();
          const user = await grok.dapi.users.find(row.uid);
          const package_ = (await grok.dapi.packages.filter(`package.name = "${row.package}"`).list())[0];
          cp.addPane('Details', () => {
            return ui.tableFromMap({'User': ui.render(user), 'Package': ui.render(package_),
              'From': dateFrom, 'To': dateTo});
          }, true);
          const button = ui.button('Details', async () => {
            this.uaToolbox.dateFromDD.value = dateFrom;
            this.uaToolbox.dateToDD.value = dateTo;
            this.uaToolbox.usersDD.value = row.user;
            this.uaToolbox.packagesDD.value = row.package;
            FunctionsView.viewers[0].reloadViewer({date: 'this week', groups: [row.ugid], packages: [row.package]});
            ViewHandler.changeView(FunctionsView.viewName);
          });
          button.style.height = '100%';
          button.style.verticalAlign = 'middle';
          button.style.marginLeft = '7px';
          const fPane = cp.addPane('Functions', () => {
            return ui.wait(async () => {
              const df: DG.DataFrame = await grok.data.query('UsageAnalysis:PackagesContextPaneFunctions', filter);
              const packageFunctions = await grok.dapi.functions.filter(`package.name = "${row.package}"`).list();
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
          const lPane = cp.addPane('Log events summary', () => {
            return ui.wait(async () => {
              const df = await grok.data.query('UsageAnalysis:PackagesContextPaneLogs', filter);
              const data: {[key: string]: number} = {};
              for (const r of df.rows) data[r.source] = r.count;
              const info = lPane.root.querySelector('#info') as HTMLElement;
              info.textContent = df.getCol('count').stats.sum.toString();
              info.style.removeProperty('display');
              return Object.keys(data).length ? ui.tableFromMap(data) : ui.divText('No data');
            });
          }, true);
          grok.shell.o = cp.root;
        });
      }, null, null, this.viewer,
    );
    this.viewers.push(packagesViewer);
    this.root.append(packagesViewer.root);
  }
}
