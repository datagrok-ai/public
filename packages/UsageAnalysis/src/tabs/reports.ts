import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';


const filtersStyle = {
  columnNames: ['error', 'reporter', 'assignee', 'time', 'number', 'is_resolved', 'labels'],
};

export class ReportsView extends UaView {
  currentFilterGroup: DG.FilterGroup | null;
  private filters: HTMLDivElement = ui.box();
  users: { [_: string]: any; } = {};

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Reports';
    this.currentFilterGroup = null;
    this.filters.style.maxWidth = '250px';
  }

  async initViewers(path?: string): Promise<void> {
    (await grok.dapi.users.list()).forEach((user) => {
      this.users[user.id] = {
        'avatar': user.picture,
        'name': user.friendlyName,
        'data': user,
      };
    });

    const grid = ui.wait(async () => {
      let t: DG.DataFrame;
      let viewer: DG.Grid;
      if (path != undefined) {
        const segments = path.split('/').filter((s) => s != '');
        if (segments.length > 1) {
          const reportNumber = parseInt(segments[1]);
          t = await grok.dapi.reports.getReports(reportNumber);
          if (t && t.rowCount > 0) {
            t.currentRowIdx = 0;
            this.showPropertyPanel(t);
            grok.dapi.reports.getReports()
              .then(async (r: DG.DataFrame) => {
                r.rows.removeWhere((r) => r.get('number') === reportNumber);
                const df = t.append(r);
                const newViewer = DG.Viewer.grid(df);
                this.updateDf(df, newViewer, this.users);
                viewer.root.replaceWith(newViewer.root);
              });
          }
        }
      }
      // @ts-ignore
      if (!t)
        t = await grok.dapi.reports.getReports();
      viewer = DG.Viewer.grid(t);
      this.updateDf(t, viewer, this.users);
      return viewer.root;
    });

    this.root.append(ui.splitH([
      this.filters,
      grid,
    ]));
  }

  applyStyle(viewer: DG.Grid) {
    viewer.columns.setOrder(['is_resolved', 'number', 'reporter', 'assignee', 'description', 'jira', 'labels', 'error', 'error_stack_trace', 'time', 'id']);
    viewer.col('reporter')!.cellType = 'html';
    viewer.col('reporter')!.width = 25;

    viewer.col('assignee')!.cellType = 'html';
    viewer.col('assignee')!.width = 25;

    viewer.col('jira')!.cellType = 'html';
    viewer.col('is_resolved')!.width = 30;

    viewer.col('id')!.visible = false;
    viewer.col('error_stack_trace_hash')!.visible = false;

    viewer.col('time')!.format = 'yyyy-MM-dd';
    viewer.col('time')!.width = 80;

    viewer.col('number')!.width = 35;

    viewer.col('description')!.width = 150;
    viewer.col('error')!.width = 200;
    viewer.col('error_stack_trace')!.width = 200;
  }

  updateDf(t: DG.DataFrame, viewer: DG.Grid, users: { [_: string]: any; }) {
    if (t.rowCount === 0) return viewer.root
    viewer.sort(['is_resolved', 'time'], [true, false]);
    this._scroll(viewer);
    this.applyStyle(viewer);
    viewer.onCellPrepare(async function(gc) {
      if ((gc.gridColumn.name === 'reporter' || gc.gridColumn.name === 'assignee') && gc.cell.value) {
        const user = users[gc.cell.value];
        const icon = DG.ObjectHandler.forEntity(user.data)?.renderIcon(user.data.dart);
        if (icon) {
          icon.style.top = 'calc(50% - 8px)';
          icon.style.left = 'calc(50% - 8px)';
          gc.style.element = ui.tooltip.bind(icon, () => {
            return DG.ObjectHandler.forEntity(user.data)?.renderTooltip(user.data.dart)!;
          });
        }
      }
      if (gc.gridColumn.name === 'jira' && gc.cell.value) {
        const link = ui.link(gc.cell.value, `https://reddata.atlassian.net/jira/software/c/projects/GROK/issues/${gc.cell.value}`);
        link.addEventListener('click', (e) => {
          e.preventDefault();
          window.open(link.href, '_blank');
        });
        link.style.position = 'absolute';
        link.style.top = 'calc(50% - 8px)';
        link.style.left = '3px';
        gc.style.element = ui.tooltip.bind(link, () => 'Link to JIRA ticket');
      }
    });
    const isAcknowledged = t.getCol('is_resolved');
    viewer.onCellRender.subscribe((gc) => {
      if (isAcknowledged.get(gc.cell.tableRowIndex ?? gc.cell.gridRow))
        gc.cell.style.textColor = 0xFFB8BAC0;
    });

    t.getCol('labels').setTag(DG.Tags.MultiValueSeparator, ',');
    t.onCurrentRowChanged.subscribe(async (_: any) => this.showPropertyPanel(t));
    t.onValuesChanged.subscribe(async () => {
      viewer.sort(['is_resolved', 'time'], [true, false]);
      this._scroll(viewer);
    });
    grok.events.onEvent('d4-report-resolved').subscribe((r: DG.UserReport) => {
      const idCol = t.getCol('id');
      const length = idCol.length;
      for (let i = 0; i < length; i++) {
        if (idCol.get(i) === r.id) {
          t.cell(i, 'is_resolved').value = r.isResolved;
          t.fireValuesChanged();
        }
      }
    });

    grok.events.onEvent('d4-report-ticket_created').subscribe((r: DG.UserReport) => {
      const idCol = t.getCol('id');
      const length = idCol.length;
      for (let i = 0; i < length; i++) {
        if (idCol.get(i) === r.id) {
          t.cell(i, 'jira').value = r.jiraTicket;
          t.fireValuesChanged();
        }
      }
    });

    this.reloadFilter(t);
  }

  showPropertyPanel(t: DG.DataFrame) {
    const currentRow = t.currentRowIdx;
    if (currentRow === -1) return;
    const reportId = t.getCol('id').get(currentRow);
    if (!reportId) return;
    this.updatePath(t.getCol('number').get(currentRow));
    DG.DetailedLog.showReportProperties(reportId);
  }

  reloadFilter(table: DG.DataFrame) {
    this.currentFilterGroup?.detach();
    while (this.filters.hasChildNodes())
      this.filters.removeChild(this.filters.lastChild!);
    if (table) {
      const filters_ = DG.Viewer.filters(table, filtersStyle);
      this.currentFilterGroup = new DG.FilterGroup(filters_.dart);
      this.filters.append(filters_.root);
    }
  }

  _scroll(viewer: DG.Grid): void {
    const segments = window.location.href.split('/');
    const last = segments[segments.length - 1];
    if (last !== 'reports' && /^-?\d+$/.test(last)) {
      const num = parseInt(last);
      if (num) {
        setTimeout(() => {
          const df = viewer.dataFrame;
          const numbers = df.getCol('number');
          for (let i = 0; i < numbers.length; i++)
            if (numbers.get(i) == num) {
              viewer.scrollToCell('is_resolved', i);
              df.currentRowIdx = i;
              break;
            }
        }, 200);
      }
    }
  }

  updatePath(reportNumber: number): void {
    const segments = window.location.href.split('/');
    const last = segments[segments.length - 1];
    if (last === 'reports')
      segments.push(reportNumber.toString());
    else {
      if (/^-?\d+$/.test(last))
        segments[segments.length - 1] = reportNumber.toString();
      else
        return;
    }

    window.history.pushState(
      null, 'Report ${detailedLog.reportNumber}', segments.join('/'));
  }
}
