import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';
import '../../css/usage_analysis.css';

export class AiUsageView extends UaView {
  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'AI Usage';
  }

  async initViewers(path?: string): Promise<void> {
    const grid = (name: string, queryName: string) => new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name,
      queryName,
      createViewer: (t: DG.DataFrame) => {
        const g = DG.Viewer.grid(t, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
        const ugid = g.col('ugid');
        if (ugid)
          ugid.visible = false;
        return g;
      },
    });

    const timeline = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Tokens over time',
      queryName: 'AiUsageTimeline',
      createViewer: (t: DG.DataFrame) => DG.Viewer.lineChart(t, {
        xColumnName: 'time_start',
        showXSelector: false,
        yColumnNames: ['input_tokens', 'output_tokens'],
        showYSelectors: false,
        showAggrSelectors: false,
        showSplitSelector: false,
        chartTypes: ['Line Chart'],
        title: 'Tokens over time',
      }),
    });
    const byUser = grid('By user', 'AiUsageByUser');
    const byGroup = grid('By group', 'AiUsageByGroup');
    const byFeature = grid('By feature', 'AiUsageByFeature');

    this.viewers.push(timeline, byUser, byGroup, byFeature);
    this.root.append(ui.splitV([
      ui.box(timeline.root, {style: {maxHeight: '320px'}}),
      ui.box(ui.tabControl({'By user': byUser.root, 'By group': byGroup.root, 'By feature': byFeature.root}).root),
    ]));
  }
}
