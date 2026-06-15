import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {UaView} from "./ua";
import {UaToolbox} from "../ua-toolbox";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {loadUsers, setupUserIconRenderer} from "../utils";
import {funcs} from "../package-api";

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
  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Errors';
  }

  async initViewers(path?: string): Promise<void> {
    const users = await loadUsers();

    const filters = ui.box();
    filters.classList.add('ua-errors-filters');

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
        queryName: 'TopErrors',
        createViewer: (t: DG.DataFrame) => {
          const viewer = createCountBarChart(t, 'error', 'Top errors', ERROR_BAR_COLOR);

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

    const topSources = new UaFilterableQueryViewer(
      {
        filterSubscription: this.uaToolbox.filterStream,
        name: 'Top Source',
        queryName: 'TopErrorSources',
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
    grok.shell.o = properties;
  }
}
