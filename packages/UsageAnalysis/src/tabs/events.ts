import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
// import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';
// import {TopPackagesViewer} from '../drilldown_viewers/events/top-packages-viewer';
// import {TopPackageFunctionsViewer} from '../drilldown_viewers/events/top-package-functions-viewer';
// import {TopFunctionsViewer} from '../drilldown_viewers/events/top-functions-viewer';
// import {TopSourcesViewer} from '../drilldown_viewers/events/top-sources-viewer';


export class EventsView extends UaView {
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Events';
  }

  async initViewers(): Promise<void> {
    const packagesViewer1 = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Sources',
      'EventsSources',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.scatterPlot(t, {
          x: 'time_start',
          y: 'source',
          size: 'count',
          // color: 'user',
          jitterSize: 5,
          markerMinSize: 10,
          markerMaxSize: 30,
          showColorSelector: false,
          showSizeSelector: false,
          showXSelector: false,
          showYSelector: false,
        });
        return viewer;
      }, null, null);

    const packagesViewer2 = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Unique Errors',
      'UniqueErrors',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, {
          // x: 'time_start',
          // y: 'source',
          // size: 'count',
          // color: 'user',
          // jitterSize: 5,
          // markerMinSize: 10,
          // markerMaxSize: 30,
          // showColorSelector: false,
          showStackSelector: false,
          showValueSelector: false,
          showCategorySelector: false,
        });
        return viewer;
      }, null, null);

    this.viewers.push(packagesViewer1);
    this.viewers.push(packagesViewer2);
    const block1 = ui.block([packagesViewer1.root], {style: {height: '100%'}});
    const block2 = ui.block([packagesViewer2.root], {style: {height: '100%'}});
    (block1.firstChild as HTMLElement).style.height = '100%';
    (block2.firstChild as HTMLElement).style.height = '100%';
    this.root.append(ui.divV([block1, block2]));
  }
}
