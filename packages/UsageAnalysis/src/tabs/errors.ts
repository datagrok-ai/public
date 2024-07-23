import {UaView} from "./ua";
import {UaToolbox} from "../ua-toolbox";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";

const filtersStyle = {
  columnNames: ['event_time', 'user', 'error_message', 'is_reported'],
};

const users: {[_: string]: any} = {};

export class ErrorsView extends UaView {
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Errors';
  }

  async initViewers(path?: string): Promise<void> {
    (await grok.dapi.users.list()).forEach((user) => {
      users[user.friendlyName] = {
        'avatar': user.picture,
        'name': user.friendlyName,
        'data': user,
      };
    });

    const filters = ui.box();
    filters.style.maxWidth = '250px';

    const errorViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Event Errors',
      queryName: 'EventErrors',
      processDataFrame: (t: DG.DataFrame) => {
        t.onSelectionChanged.subscribe(async () => {
          await this.showErrorContextPanel(t);
        });
        t.onCurrentRowChanged.subscribe(async () => {
          t.selection.setAll(false);
          t.selection.set(t.currentRowIdx, true);
        });
        return t;
      },
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t, {
          'showColumnLabels': false,
          'showRowHeader': false,
          'showColumnGridlines': false,
          'allowRowSelection': false,
          'allowBlockSelection': false,
          'showCurrentCellOutline': false,
          'defaultCellFont': '13px monospace'
        });
        filters.appendChild(DG.Viewer.filters(t, filtersStyle).root);

        viewer.col('user')!.cellType = 'html';
        viewer.col('user')!.width = 30;
        viewer.col('id')!.visible = false;

        viewer.onCellPrepare(async function(gc) {
          if (gc.gridColumn.name === 'event_time') {
            gc.style.textColor = 0xFFB8BAC0;
            gc.style.font = '13px Roboto';
          }

          if (gc.gridColumn.name === 'user') {
            gc.gridColumn.tooltipType = 'Columns';
            const user = users[gc.cell.value];
            if (!user) return;
            const img = ui.div();
            img.style.width = '20px';
            img.style.height = '20px';
            img.style.backgroundSize = 'contain';
            img.style.margin = '5px 0 0 5px';
            img.style.borderRadius = '100%';
            if (gc.cell.value != 'Test')
              img.style.backgroundImage = user.avatar.style.backgroundImage;
            else
              img.style.backgroundImage = 'url(/images/entities/grok.png);';
            img.addEventListener('click', () => {
              grok.shell.o = user.data;
            });
            gc.style.element = ui.tooltip.bind(img, user.name);
          }
        });

        return viewer;
      },
    });

    const topErrors = new UaFilterableQueryViewer(
      {
        filterSubscription: this.uaToolbox.filterStream,
        name: 'Top Errors',
        queryName: 'TopErrors',
        createViewer: (t: DG.DataFrame) => {
          const viewer =  DG.Viewer.barChart(t, {
            'valueColumnName': 'count',
            'valueAggrType': 'sum',
            'barSortType': 'by value',
            'barSortOrder': 'desc',
            'showValueAxis': false,
            'showValueSelector': false,
            'splitColumnName': 'error',
            'showCategoryValues': false,
            'showCategorySelector': false,
            'stackColumnName': '',
            'showStackSelector': false,
            'title': 'Top errors'
          });

          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            const df: DG.DataFrame | undefined = errorViewer.viewer?.dataFrame;
            if (df) {
              df.filter.handleClick((i) => {
                const column = df.getCol('error_message');
                return column.get(i) == args.args.options.categories[0];
              }, new MouseEvent(''));
            }
          });
          return viewer;
        }
      }
    );

    errorViewer.root.classList.add('ui-panel');
    this.viewers.push(errorViewer);
    this.viewers.push(topErrors);
    this.root.append(ui.splitV([
      ui.splitH([
        filters,
        errorViewer.root
      ]),
      ui.box(topErrors.root, {style: {maxHeight: '250px'}})
    ]));
  }

  showErrorContextPanel(table: DG.DataFrame): void {
    if (!table.selection.anyTrue) return;
    let df = table.clone(table.selection);
    const eventId = df.getCol('id').get(0);
    if (!eventId) return;
    const accordion = DG.Accordion.create();
    const properties = ui.div([accordion.root]);

    accordion.addPane('Details', () => ui.wait(async () => {
      const entity: DG.LogEvent = await grok.dapi.log.find(eventId);
      return ui.tableFromMap({
        'Error message': df.getCol('error_message').get(0),
        'Stack trace': df.getCol('error_stack_trace').get(0),
        'Handled': entity.parameters.find((p) => p.parameter.name === 'handled')?.value,
        'Source': entity.parameters.find((p) => p.parameter.name === 'source')?.value,
        'User': users[df.getCol('user').get(0)]?.data,
        'Reported': df.getCol('is_reported').get(0),
      });
    }), true);

    accordion.addPane('Statistics', () => ui.wait(async () => {
      const detailsButton = ui.button('Details', async () => {
        const ev = this.uaToolbox.viewHandler.getView('Reports');
        const viewer = ev.viewers[0];
        viewer.reloadViewer({'event_id': eventId});
        if (!viewer.activated)
          viewer.activated = true;
        this.uaToolbox.viewHandler.changeTab('Reports');
        this.uaToolbox.drilldown = this.uaToolbox.viewHandler.getCurrentView();
      });
      detailsButton.classList.add('ua-details-button');
      const promises: Promise<any>[] = [
        grok.functions.call('UsageAnalysis:ReportsCount', {'event_id': eventId}),
        grok.functions.call('UsageAnalysis:SameErrors', {'event_id': eventId}),
      ];
      const results = await Promise.all(promises);
      let div: HTMLDivElement;
      // const createReportButton = ui.button('Report', async () => {
      //   const prog = DG.TaskBarProgressIndicator.create('Reporting...');
      //   await window.grok_Dapi_Admin_Post_Event_Report(this.dart, eventId)(eventId);
      //   const ev = ViewHandler.getView('Reports');
      //   ev.viewers[0].reloadViewer();
      //   while (div.hasChildNodes())
      //     div.removeChild(div.lastChild!);
      //   div.appendChild(ui.span([1]));
      //   div.appendChild(detailsButton);
      //   prog.close();
      // });
      div = ui.divH([ui.span([results[0]]), results[0] > 0 ? detailsButton : null]);
      div.style.alignItems = 'center';
      const map = {'Reports': div, 'Same errors': results[1]};
      return ui.tableFromMap(map);
    }));
    grok.shell.o = properties;
  }
}
