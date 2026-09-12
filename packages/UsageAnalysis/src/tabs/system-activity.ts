import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';
import {loadUsers, setupUserIconRenderer, showEventDetails} from '../utils';
import '../../css/usage_analysis.css';

export class SystemActivityView extends UaView {
  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'System Activity';
  }

  async initViewers(path?: string): Promise<void> {
    const users = await loadUsers();
    const filters = ui.box();
    filters.classList.add('ua-filters');

    const summaryViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'System events over time',
      queryName: 'SystemActivitySummary',
      createViewer: (t: DG.DataFrame) => DG.Viewer.lineChart(t, {
        xColumnName: 'time_start',
        showXSelector: false,
        yColumnNames: ['count'],
        showYSelectors: false,
        showAggrSelectors: false,
        showSplitSelector: false,
        chartTypes: ['Line Chart'],
        title: 'System events',
        split: 'event',
      }),
    });

    const activityViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'System activity',
      queryName: 'SystemActivity',
      processDataFrame: (t: DG.DataFrame) => {
        t.onCurrentRowChanged.subscribe(() => showEventDetails(t));
        return t;
      },
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t, {
          showRowHeader: false,
          allowRowSelection: false,
          allowBlockSelection: false,
        });
        ui.empty(filters);
        filters.append(DG.Viewer.filters(t, {columnNames: ['event', 'user']}).root);
        viewer.columns.setOrder(['event_time', 'event', 'user', 'description', 'details']);
        viewer.col('ugid')!.visible = false;
        viewer.col('id')!.visible = false;
        setupUserIconRenderer(viewer, users, ['user']);
        return viewer;
      },
    });

    activityViewer.root.classList.add('ui-panel');
    this.viewers.push(summaryViewer, activityViewer);
    this.root.append(ui.splitV([
      summaryViewer.root,
      ui.splitH([filters, activityViewer.root]),
    ], {}, true));
  }
}
