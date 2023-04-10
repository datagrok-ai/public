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

  constructor(uaToolbox: UaToolbox) {
    super();
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;
    this.box = true;
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
