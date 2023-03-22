import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
// import {UaDataFrameQueryViewer} from '../viewers/ua-data-frame-query-viewer';
// import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';
// import {PropertyPanel} from '../property-panel';
import {ViewHandler} from '../view-handler';
import {FunctionsView} from './functions';


export class PackagesView extends UaView {
  static viewName = 'Packages';

  constructor(uaToolbox: UaToolbox, view: DG.View) {
    super(uaToolbox, view);
  }

  async initViewers(): Promise<void> {
    const packagesViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Packages',
      'PackagesUsage',
      (t: DG.DataFrame) => {
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
        t.onCurrentRowChanged.subscribe(async () => {
          const rowValues = Array.from(t.currentRow.cells).map((c) => c.value);
          const row = Object.fromEntries(t.columns.names().map((k, i) => [k, rowValues[i]]));
          row.time_start = row.time_start.a;
          row.time_end = row.time_end.a;
          const filter = {time_start: row.time_start / 1000, time_end: row.time_end / 1000, user: row.uid,
            package: row.package, ...this.uaToolbox.filterStream.getValue()};
          const cp = DG.Accordion.create();
          const user = await grok.dapi.users.find(row.uid);
          const package_ = (await grok.dapi.packages.filter(`package.name = "${row.package}"`).list())[0];
          cp.addPane('Details', () => {
            return ui.tableFromMap({'User': ui.render(user), 'Package': ui.render(package_),
              'From': new Date(row.time_start).toLocaleString(),
              'To': new Date(row.time_end).toLocaleString()});
          }, true);
          const button = ui.button('Details', async () => {
            this.uaToolbox.groupsInput.choices.removeActiveItems(-1);
            const userGroup = (await grok.dapi.groups.find(row.ugid)).friendlyName;
            this.uaToolbox.groupsInput.choices._addItem({value: row.ugid, label: userGroup});
            this.uaToolbox.packagesInput.choices.removeActiveItems(-1);
            this.uaToolbox.packagesInput.choices._addItem({value: row.package, label: row.package});
            this.uaToolbox.applyFilter();
            ViewHandler.changeTab(FunctionsView.viewName);
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
        return viewer.root;
      },
    );
    this.viewers.push(packagesViewer);
    this.root.append(packagesViewer.root);
  }
}
