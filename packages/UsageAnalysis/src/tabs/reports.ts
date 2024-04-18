import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import * as grok from 'datagrok-api/grok';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import * as DG from 'datagrok-api/dg';
import {DetailedLog} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ViewHandler} from "../view-handler";

const filtersStyle = {
  columnNames: ['error', 'reporter', 'report_time', 'report_number', 'is_resolved'],
};

export class ReportsView extends UaView {
  currentFilterGroup: DG.FilterGroup | null;
  private filters: HTMLDivElement = ui.box();

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Reports';
    this.currentFilterGroup = null;
    this.filters.style.maxWidth = '250px';
    this.ribbonMenu = DG.Menu.create();
    this.ribbonMenu.item('Text', () => {
      console.log('From ribbon');
    });
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

    const reportsViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'User reports',
      queryName: 'UserReports',
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t, {
          'defaultCellFont': '13px monospace'
        });
        this.reloadFilter(t);
        viewer.onBeforeDrawContent.subscribe(() => {
          viewer.columns.setOrder(['is_resolved', 'report_number', 'reporter', 'assignee', 'description', 'same_errors_count', 'requests_count', 'error', 'error_stack_trace', 'report_time', 'report_id']);
          viewer.col('reporter')!.cellType = 'html';
          viewer.col('assignee')!.cellType = 'html';
          viewer.col('is_resolved')!.editable = false;
          viewer.col('requests_count')!.editable = false;
          viewer.col('same_errors_count')!.editable = false;
          viewer.col('error_stack_trace_hash')!.visible = false;
        });
        const isResolvedCol = t.getCol('is_resolved');
        viewer.onCellPrepare(async function(gc) {
          if (isResolvedCol.get(gc.gridRow))
            gc.style.textColor = 0xFFB8BAC0;
          if ((gc.gridColumn.name === 'reporter' || gc.gridColumn.name === 'assignee') && gc.cell.value) {
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
            const span = ui.span([ui.span([img, ui.label(user.name)], 'grok-markup-user')], 'd4-link-label');
            span.addEventListener('click', () => {
              grok.shell.o = user.data;
            });
            gc.style.element = ui.tooltip.bind(span, () => {
              return DG.ObjectHandler.forEntity(user.data)?.renderTooltip(user.data.dart)!;
            });
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
      this.updateFilter();
    }
  }

  updateFilter(): void {
    if (ViewHandler.getCurrentView().name === 'Reports' && ViewHandler.getInstance().getSearchParameters().has('report-number')) {
      let reportNumber = ViewHandler.getInstance().getSearchParameters().get('report-number');
      if (reportNumber) {
        this.currentFilterGroup?.updateOrAdd({
          type: DG.FILTER_TYPE.MULTI_VALUE,
          column: 'report_number',
          selected: [reportNumber]}
        );
      }
    }
  }

  async showReportContextPanel(table: DG.DataFrame): Promise<void> {
    if (!table.selection.anyTrue) return;
    const index = table.selection.getSelectedIndexes()[0];
    let df = table.clone(table.selection);
    const reportId = df.getCol('report_id').get(0);
    if (!reportId) return;
    grok.shell.o = ui.wait(async () => {
      const accordion = await DetailedLog.getAccordion(reportId);
      const pane = accordion.root.querySelectorAll('[d4-title="ACTIONS"]')[0];
      let reopen: HTMLLabelElement;
      let close: HTMLLabelElement;
      let isActive: boolean = false;

      async function resolve(e: MouseEvent, value: boolean, replacement: HTMLElement, action: Function) {
        if (isActive) return;
        isActive = true;
        (e.target as HTMLElement).replaceWith(ui.wait(async () => {
          (e.target as HTMLElement).style.color = '#286344';
          await action();
          table.getCol('is_resolved').set(index, value);
          isActive = false;
          return replacement;
        }));
      }
      reopen = ui.actionLink('Reopen issue...', async (e: MouseEvent) => {
        const dialog = DG.Dialog.create('Specify the reason');
        const description = ui.input.textArea('Description', {'size': {height: 150, width: 250}, 'nullable': false, 'tooltipText': 'Provide reason of reopening'});
        const notify = ui.input.bool('Notify assignee', {'tooltipText': 'Email with notification will be sent to the assignee'});
        dialog.add(description);
        dialog.add(notify);
        dialog.onOK(async () => await resolve(e, false, close, async () => await grok.dapi.reports.reopen(reportId, description.value, notify.value)));
        dialog.show();
      });
      close = ui.actionLink('Mark as resolved...', async (e: MouseEvent) => {
        await resolve(e, true, reopen, async () => await grok.dapi.reports.resolve(reportId));
      });
      pane.children[1].children[0].append(df.getCol('is_resolved').get(0) ? reopen : close);
      return ui.div([accordion.root]);
    });
  }
}
