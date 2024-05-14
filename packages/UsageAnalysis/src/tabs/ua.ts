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
  systemId: string = '00000000-0000-0000-0000-000000000000';
  rout?: string;

  constructor(uaToolbox: UaToolbox) {
    super();
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;
    this.box = true;
  }

  async tryToinitViewers(path?: string): Promise<void> {
    if (!this.initialized) {
      this.initialized = true;
      await this.initViewers(path);
    }
  }

  async initViewers(path?: string): Promise<void> {}

  switchRout(): void {}
}
