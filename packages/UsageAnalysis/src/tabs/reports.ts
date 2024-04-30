import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import * as grok from 'datagrok-api/grok';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import * as DG from 'datagrok-api/dg';
import {DetailedLog, Grid} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {delay} from "rxjs/operators";
import {div} from "datagrok-api/ui";
import {ViewHandler} from "../view-handler";

const filtersStyle = {
  columnNames: ['error', 'reporter', 'assignee', 'time', 'report_number', 'is_acknowledged', 'label'],
};

export class ReportsView extends UaView {
  currentFilterGroup: DG.FilterGroup | null;
  private filters: HTMLDivElement = ui.box();

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Reports';
    this.currentFilterGroup = null;
    this.filters.style.maxWidth = '250px';
  }

  async initViewers(path?: string): Promise<void> {
    const users: {[_: string]: any} = {};
    (await grok.dapi.users.list()).forEach((user) => {
      users[user.friendlyName] = {
        'avatar': user.picture,
        'name': user.friendlyName,
        'data': user,
      };
    });

    const reportsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'User reports',
      queryName: 'UserReports',
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t);
        viewer.getInfo()
        viewer.sort(['is_acknowledged', 'errors', 'time'], [true, false, true]);
        this.reloadFilter(t);
        if (path != undefined) {
          const segments = path.split('/').filter((s) => s != '');
          if (segments.length > 1) {
            this.setReportFilter(segments[1]);
            this.updatePath(segments[1]);
          }
        }
        viewer.onBeforeDrawContent.subscribe(() => {
          viewer.columns.setOrder(['is_acknowledged', 'report_number', 'reporter', 'assignee', 'description', 'errors', 'jira', 'label', 'error', 'error_stack_trace', 'time', 'report_id']);
          viewer.col('reporter')!.cellType = 'html';
          viewer.col('reporter')!.width = 25;

          viewer.col('assignee')!.cellType = 'html';
          viewer.col('assignee')!.width = 25;
          viewer.col('errors')!.width = 35;
          viewer.col('errors')!.editable = false;

          viewer.col('jira')!.cellType = 'html';
          viewer.col('is_acknowledged')!.width = 30;

          viewer.col('report_id')!.visible = false;
          viewer.col('error_stack_trace_hash')!.visible = false;

          viewer.col('time')!.format = 'yyyy-MM-dd';
          viewer.col('time')!.width = 80;

          viewer.col('report_number')!.width = 35;

          viewer.col('description')!.width = 150;
          viewer.col('error')!.width = 200;
          viewer.col('error_stack_trace')!.width = 200;
        });

        viewer.onCellPrepare(async function(gc) {
          if ((gc.gridColumn.name === 'reporter' || gc.gridColumn.name === 'assignee') && gc.cell.value) {
            const user = users[gc.cell.value];
            const icon = DG.ObjectHandler.forEntity(user.data)?.renderIcon(user.data.dart)!;
            icon.style.top = 'calc(50% - 8px)';
            icon.style.left = 'calc(50% - 8px)';
            gc.style.element = ui.tooltip.bind(icon, () => {
              return DG.ObjectHandler.forEntity(user.data)?.renderTooltip(user.data.dart)!;
            });
          }
          if (gc.gridColumn.name === 'jira' && gc.cell.value) {
            const link = ui.link(gc.cell.value, `https://reddata.atlassian.net/jira/software/c/projects/GROK/issues/${gc.cell.value}`);
            link.addEventListener('click', (e) => {
              e.preventDefault();
              window.open(link.href, '_blank');
            });
            gc.style.element = ui.tooltip.bind(link, () => 'Link to JIRA ticket');
          }
        });
        const isAcknowledged = t.getCol('is_acknowledged');
        viewer.onCellRender.subscribe((gc) => {
          if (isAcknowledged.get(gc.cell.tableRowIndex ?? gc.cell.gridRow))
            gc.cell.style.textColor = 0xFFB8BAC0;
        });

        return viewer;
      },
      processDataFrame: (t: DG.DataFrame) => {
        t.getCol('label').setTag(DG.Tags.MultiValueSeparator, ',');
        t.onCurrentRowChanged.subscribe(async (_) => {
          const currentRow = t.currentRowIdx;
          if (currentRow === null || currentRow === undefined || currentRow === -1) return;
          const reportId = t.getCol('report_id').get(currentRow);
          if (!reportId) return;
          this.updatePath(t.getCol('report_number').get(currentRow));
          DetailedLog.showReportProperties(reportId, t, currentRow);
        });
        t.onValuesChanged.subscribe(async () => {
          await delay(500);
          (this.viewers[0].viewer as Grid)
            .sort(['is_acknowledged', 'errors', 'time'], [true, false, true]);
        });
        t.onCurrentCellChanged.subscribe(() => {
          if (t.currentCol.name === 'errors') {
            const errorHash = t.getCol('error_stack_trace_hash').get(t.currentRowIdx);
            const error = t.getCol('error').get(t.currentRowIdx);
            grok.shell.o = div([
              ui.actionLink('Get all errors...', async () => {
                const progress = DG.TaskBarProgressIndicator.create('Receiving errors...');
                try {
                  const result = await grok.functions
                    .call('UsageAnalysis:ReportSameErrors', {'stackTraceHash': errorHash, 'errorMessage': error});
                  grok.shell.addTableView(result);
                } catch (e: any) {
                  grok.shell.error(e.toString());
                } finally {
                  progress.close();
                }
              }),
            ]);
          }
        });
        return t;
      },
    });

    reportsViewer.root.classList.add('ui-panel');
    this.viewers.push(reportsViewer);
    this.root.append(ui.splitH([
      this.filters,
      reportsViewer.root,
    ]));
  }

  async reloadViewers(): Promise<void> {
    this.viewers = [];
    while (this.root.hasChildNodes())
      this.root.removeChild(this.root.lastChild!);
    await this.initViewers();
    for (const v of this.viewers)
      await v.reloadViewer();
  }

  async reloadFilter(table?: DG.DataFrame) {
    this.currentFilterGroup?.detach();
    while (this.filters.hasChildNodes())
      this.filters.removeChild(this.filters.lastChild!);
    table = table ?? await this.viewers[0].dataFrame;
    if (table) {
      const filters_ = DG.Viewer.filters(table, filtersStyle);
      this.currentFilterGroup = new DG.FilterGroup(filters_.dart);
      this.filters.append(filters_.root);
      //this.updateFilter();
    }
  }

  updatePath(reportNumber: string): void {
    const segments = window.location.href.split('/');
    const last = segments[segments.length - 1];
    if (last === 'reports')
      segments.push(reportNumber);
    else {
        if (/^-?\d+$/.test(last))
          segments[segments.length - 1] = reportNumber;
        else
          return;
      }

    window.history.pushState(
      null, 'Report ${detailedLog.reportNumber}', segments.join('/'));
  }

  setReportFilter(reportNumber: string): void {
    console.log(reportNumber);
    if (reportNumber != undefined && reportNumber != '') {
      this.currentFilterGroup?.updateOrAdd({
        type: DG.FILTER_TYPE.MULTI_VALUE,
        column: 'report_number',
        selected: [reportNumber],
      });
    }
  }
}
