import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {UaView} from "./ua";
import {UaToolbox} from "../ua-toolbox";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {UaFilter} from "../filter";
import {loadUsers, setupUserIconRenderer} from "../utils";
import {funcs, queries} from "../package-api";

import '../../css/usage_analysis.css';

const filtersStyle = {
  columnNames: ['event_time', 'user', 'error_message', 'is_reported'],
};

const ERROR_BAR_COLOR = 0xFFD9534F;
const SOURCE_BAR_COLOR = 0xFF5CB85C;

function createCountBarChart(t: DG.DataFrame, splitColumnName: string, title: string, barColor: number): DG.Viewer {
  return DG.Viewer.barChart(t, {
    'valueColumnName': 'count',
    'valueAggrType': 'sum',
    'barSortType': 'by value',
    'barSortOrder': 'desc',
    'showValueAxis': false,
    'showValueSelector': false,
    'splitColumnName': splitColumnName,
    'showCategoryValues': false,
    'showCategorySelector': false,
    'stackColumnName': '',
    'showStackSelector': false,
    'title': title,
    'barColor': barColor,
  });
}

export class ErrorsView extends UaView {
  private metrics?: {filter: UaFilter, errors: Promise<DG.ServerMetrics['errors']>};

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Errors';
  }

  /** The errors of `/admin/metrics` for the applied filter's date, by signature, fetched once for both charts. */
  errorMetrics(): Promise<DG.ServerMetrics['errors']> {
    const filter = this.uaToolbox.filterStream.value;
    if (this.metrics?.filter !== filter)
      this.metrics = {filter, errors: grok.dapi.admin.getMetrics({date: filter.date, limit: 15, errorsBy: 'signature'})
        .then((m) => m.errors)};
    return this.metrics.errors;
  }

  async initViewers(path?: string): Promise<void> {
    const users = await loadUsers();

    const filters = ui.box();
    filters.classList.add('ua-filters');

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

        viewer.col('id')!.visible = false;
        viewer.col('error_stack_trace_hash')!.visible = false;

        viewer.onCellPrepare((gc) => {
          if (gc.gridColumn.name === 'event_time') {
            gc.style.textColor = 0xFFB8BAC0;
            gc.style.font = '13px Roboto';
          }
        });
        setupUserIconRenderer(viewer, users, ['user']);

        return viewer;
      },
    });

    const topErrors = new UaFilterableQueryViewer(
      {
        filterSubscription: this.uaToolbox.filterStream,
        name: 'Top Errors',
        getDataFrame: async () => {
          const top = (await this.errorMetrics()).top;
          return DG.DataFrame.fromColumns([
            DG.Column.fromStrings('error', top.map((e) => e.message)),
            DG.Column.fromStrings('signature', top.map((e) => e.signature ?? '')),
            DG.Column.fromList(DG.COLUMN_TYPE.INT, 'count', top.map((e) => e.count)),
          ]);
        },
        createViewer: (t: DG.DataFrame) => {
          const viewer = createCountBarChart(t, 'error', 'Top errors', ERROR_BAR_COLOR);

          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            const df: DG.DataFrame | undefined = errorViewer.viewer?.dataFrame;
            const error = args.args.options.categories[0];
            const signatures = new Set([...Array(t.rowCount).keys()].filter((i) => t.get('error', i) === error)
              .map((i) => t.get('signature', i)));
            if (df)
              df.filter.handleClick((i) => signatures.has(df.get('error_stack_trace_hash', i)), new MouseEvent(''));
          });
          return viewer;
        }
      }
    );

    const topSources = new UaFilterableQueryViewer(
      {
        filterSubscription: this.uaToolbox.filterStream,
        name: 'Top Source',
        getDataFrame: async () => {
          const bySource = (await this.errorMetrics()).bySource;
          return DG.DataFrame.fromColumns([
            DG.Column.fromStrings('error_source', Object.keys(bySource)),
            DG.Column.fromList(DG.COLUMN_TYPE.INT, 'count', Object.values(bySource)),
          ]);
        },
        createViewer: (t: DG.DataFrame) =>
          createCountBarChart(t, 'error_source', 'Top source', SOURCE_BAR_COLOR),
      }
    );

    const errorsSummary = new UaFilterableQueryViewer(
      {
        filterSubscription: this.uaToolbox.filterStream,
        name: 'Errors Summary',
        queryName: 'EventsSources',
        processDataFrame: (t: DG.DataFrame) =>
          t.clone(DG.BitSet.create(t.rowCount, (i) => t.getCol('source').get(i) === 'error')),
        createViewer: (t: DG.DataFrame) => {
          return DG.Viewer.lineChart(t, {
            'xColumnName': 'time_start',
            'yColumnNames': ['count'],
            'showXSelector': false,
            'showYSelectors': false,
            'showAggrSelectors': false,
            'showSplitSelector': false,
            'chartTypes': ['Line Chart'],
            'lineColoringType': 'Custom',
            'lineColor': ERROR_BAR_COLOR,
            'markerColor': ERROR_BAR_COLOR,
            'title': 'Errors Summary'
          });
        }
      }
    );

    errorViewer.root.classList.add('ui-panel');
    this.viewers.push(errorViewer, errorsSummary, topErrors, topSources);
    this.root.append(ui.splitV([
      errorsSummary.root,
      ui.box(ui.splitH([topErrors.root, topSources.root]), {style: {maxHeight: '250px'}}),
      ui.splitH([
        filters,
        errorViewer.root
      ])
    ]));
  }

  showErrorContextPanel(table: DG.DataFrame): void {
    if (!table.selection.anyTrue) return;
    const rowIdx = table.selection.getSelectedIndexes()[0];
    const eventId = table.getCol('id').get(rowIdx);
    if (!eventId) return;
    const accordion = DG.Accordion.create();
    const properties = ui.div([accordion.root]);

    accordion.addPane('Details', () => ui.wait(async () => {
      const entity: DG.LogEvent = await grok.dapi.log.find(eventId);
      const users = await loadUsers();
      return ui.tableFromMap({
        'Error message': table.getCol('error_message').get(rowIdx),
        'Stack trace': table.getCol('error_stack_trace').get(rowIdx),
        'Handled': entity.parameters.find((p) => p.parameter?.name === 'handled')?.value,
        'Source': entity.parameters.find((p) => p.parameter?.name === 'source')?.value,
        'User': users[table.getCol('user').get(rowIdx)],
        'Reported': table.getCol('is_reported').get(rowIdx),
      });
    }), true);

    accordion.addPane('Statistics', () => ui.wait(async () => {
      const promises: Promise<any>[] = [
        grok.functions.call('UsageAnalysis:ReportsCount', {'event_id': eventId}),
        grok.functions.call('UsageAnalysis:SameErrors', {'event_id': eventId}),
      ];
      const results = await Promise.all(promises);
      const { count, report_number: reportNumber } = results[0];
      const detailsButton = ui.button('Details', async () => {
        grok.shell.addView(await funcs.reportsApp(`/${reportNumber}`));
      });
      detailsButton.classList.add('ua-details-button');
      const div = ui.divH([ui.span([count]), count > 0 ? detailsButton : null]);
      div.classList.add('ua-errors-reports');
      const map = {'Reports': div, 'Same errors': results[1]};
      return ui.tableFromMap(map);
    }));

    accordion.addPane('Alert', () => ui.wait(async () => {
      const t = await queries.errorAlerts(table.getCol('error_stack_trace_hash').get(rowIdx));
      if (t.rowCount === 0)
        return ui.divText('No error incident');
      return ui.tableFromMap({'Problem': t.get('problem', 0), 'Alert': t.get('alert', 0) ?? 'none',
        'Summary': t.get('summary', 0), 'Alerted': t.get('alerted_at', 0), 'Last seen': t.get('last_seen', 0)});
    }));
    grok.shell.o = properties;
  }
}
