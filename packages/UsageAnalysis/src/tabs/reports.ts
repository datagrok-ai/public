import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import * as grok from 'datagrok-api/grok';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import * as DG from 'datagrok-api/dg';
import {DetailedLog} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

const filtersStyle = {
  columnNames: ['error', 'reporter', 'report_time'],
};

export class ReportsView extends UaView {
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Reports';
  }

  async initViewers(): Promise<void> {
    const users: {[_: string]: any} = {};
    (await grok.dapi.users.list()).forEach((user) => {
      users[user.friendlyName] = {
        'avatar': user.picture,
        'name': user.friendlyName,
        'data': user,
      };
    });

    const filters = ui.box();
    filters.style.maxWidth = '250px';

    const reportsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'User reports',
      queryName: 'UserReports',
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t, {
          'showColumnLabels': false,
          'showRowHeader': false,
          'showColumnGridlines': false,
          'allowRowSelection': false,
          'allowColumnSelection': false,
          'allowBlockSelection': false,
          'showCurrentCellOutline': false,
          'defaultCellFont': '13px monospace',
        });
        filters.append(DG.Viewer.filters(t, filtersStyle).root);
        viewer.onBeforeDrawContent.subscribe(() => {
          viewer.columns.setOrder(['report_id', 'reporter', 'report_time', 'description', 'same_errors_count', 'error']);
          viewer.col('reporter')!.cellType = 'html';
          viewer.col('reporter')!.width = 30;
          viewer.col('report_id')!.cellType = 'html';
          viewer.col('report_id')!.width = 20;
          viewer.col('report_number')!.visible = false;
        });

        viewer.onCellPrepare(async function(gc) {
          if (gc.gridColumn.name === 'report_time') {
            gc.style.textColor = 0xFFB8BAC0;
            gc.style.font = '13px Roboto';
          }

          if (gc.gridColumn.name === 'reporter') {
            const user = users[gc.cell.value];
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

          if (gc.gridColumn.name === 'report_id') {
            const icon = ui.iconFA('arrow-to-bottom', () => {
              const indicator = DG.TaskBarProgressIndicator.create('Receiving user report...');
              //@ts-ignore
              grok.dapi.admin.getUserReport(gc.cell.value)
                .then((bytes: any) => DG.Utils.download(`report_${gc.cell.value}.zip`, bytes))
                .finally(() => indicator.close());
            });
            icon.style.marginTop = '7px';
            gc.style.element = ui.tooltip.bind(icon, 'Download report');
          }
        });

        return viewer;
      },
      processDataFrame: (t: DG.DataFrame) => {
        t.onSelectionChanged.subscribe(async () => {
          await this.showReportContextPanel(t);
        });
        t.onCurrentRowChanged.subscribe(async () => {
          t.selection.setAll(false);
          t.selection.set(t.currentRowIdx, true);
        });
        return t;
      },
    });

    reportsViewer.root.classList.add('ui-panel');
    this.viewers.push(reportsViewer);
    this.root.append(ui.splitH([
      filters,
      reportsViewer.root,
    ]));
  }

  async showReportContextPanel(table: DG.DataFrame): Promise<void> {
    if (!table.selection.anyTrue) return;
    let df = table.clone(table.selection);
    const reportId = df.getCol('report_id').get(0);
    if (!reportId) return
    grok.shell.o = ui.div([await DetailedLog.getAccordion(reportId)]);
  }
}
