import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';


export class EventsView extends UaView {
  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Events';
  }

  async initViewers(path?: string): Promise<void> {
    const packagesViewer1 = new UaFilterableQueryViewer( {
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Sources',
      queryName: 'EventsSources',
      createViewer: (t: DG.DataFrame) => {
        return DG.Viewer.lineChart(t, {
          // 'overviewColumnName': 'date',
          xColumnName: 'time_start',
          showXSelector: false,
          yColumnNames: ['count'],
          showYSelectors: false,
          showAggrSelectors: false,
          showSplitSelector: false,
          // 'showMarkers': 'Never',
          chartTypes: ['Line Chart'],
          title: 'Sources',
          split: 'source',
        });
      }});

    const packagesViewer2 = new UaFilterableQueryViewer( {
      filterSubscription: this.uaToolbox.filterStream,
      name: 'User events',
      queryName: 'EventsUsersSources',
      createViewer: (t: DG.DataFrame) => {
        return DG.Viewer.scatterPlot(t, {
          x: 'time_start',
          y: 'source',
          size: 'count',
          color: 'user',
          jitterSize: 5,
          markerMinSize: 10,
          markerMaxSize: 30,
          showColorSelector: false,
          showSizeSelector: false,
          showXSelector: false,
          showYSelector: false,
          title: 'User events',
          invertYAxis: true,
        });
      }});

    this.viewers.push(packagesViewer1);
    this.viewers.push(packagesViewer2);
    const block1 = ui.block([packagesViewer1.root], {style: {height: '100%'}});
    const block2 = ui.block([packagesViewer2.root], {style: {height: '100%'}});
    (block1.firstChild as HTMLElement).style.height = '100%';
    (block2.firstChild as HTMLElement).style.height = '100%';
    this.root.append(ui.divV([block1, block2]));
  }
}
