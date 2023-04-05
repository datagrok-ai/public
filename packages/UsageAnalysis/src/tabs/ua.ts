import * as DG from 'datagrok-api/dg';
// import * as ui from 'datagrok-api/ui';

import {UaToolbox} from '../ua-toolbox';
import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';


export interface Filter {
  time_start: number;
  time_end: number;
  groups?: string[];
  users: string[];
  packages: string[];
  functions?: string[];
}


export class UaView extends DG.ViewBase {
  uaToolbox: UaToolbox;
  viewers: UaQueryViewer[] = [];
  initialized: boolean = false;
  viewer: DG.Viewer;

  constructor(uaToolbox: UaToolbox, y?: string, ColorSelector: boolean = false) {
    super();
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;
    this.box = true;
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'time_start', []),
      DG.Column.fromList('string', y ?? '', []),
      DG.Column.fromList('string', 'count', []),
      DG.Column.fromList('string', 'user', [])]);
    this.viewer = DG.Viewer.scatterPlot(df, {
      x: 'time_start',
      y: y,
      size: 'count',
      // color: 'user',
      jitterSize: 5,
      markerMinSize: 10,
      markerMaxSize: 30,
      showColorSelector: ColorSelector,
      showSizeSelector: false,
      showXSelector: false,
      showYSelector: false,
    });
  }

  tryToinitViewers() {
    if (!this.initialized) {
      this.initialized = true;
      this.initViewers();
    }
  }

  getScatterPlot() {
    return this.viewers[0];
  }

  async initViewers(): Promise<void> {}
}
