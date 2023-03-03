// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../css/usage_analysis.css';
import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua-view';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';


export class PackagesView extends UaView {
  static viewName = 'Packages';

  constructor(uaToolbox: UaToolbox) {
    super(uaToolbox);
  }

  async initViewers(): Promise<void> {
    const packagesViewer = new UaFilterableQueryViewer(
      this.uaToolbox.filterStream,
      'Packages',
      'PackagesUsage',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.scatterPlot(t, {
          x: 'time',
          y: 'package',
          size: 'count',
          color: 'user',
        }).root;
        // viewer.style.height = '1000px !important';
        return viewer;
      },
      // (host: HTMLElement) => {
      //   host.style.overflow='auto';
      //   host.style.height='1000px !important';
      // },
    );
    this.viewers.push(packagesViewer);

    this.root.append(
      packagesViewer.root,
    );
  }
}
