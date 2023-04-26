import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';


export class EventsView extends UaView {
  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
    this.name = 'Events';
  }

  async initViewers(): Promise<void> {
    const packagesViewer1 = new UaFilterableQueryViewer( {
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Sources',
      queryName: 'EventsSources',
      viewerFunction: (t: DG.DataFrame) => {
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
          title: 'Sources',
        });
        return viewer;
      }});

    const packagesViewer2 = new UaFilterableQueryViewer( {
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Unique Errors',
      queryName: 'UniqueErrors',
      viewerFunction: (t: DG.DataFrame) => {
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
          title: 'Unique Errors',
        });
        return viewer;
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
