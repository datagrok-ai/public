import * as DG from 'datagrok-api/dg';

import {UaToolbox} from '../ua-toolbox';
import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';

export class UaView extends DG.ViewBase {
  static viewName = 'Usage Analysis View';
  uaToolbox: UaToolbox;
  viewers: UaQueryViewer[] = [];
  initialized: boolean = false;

  constructor(uaToolbox: UaToolbox) {
    super();
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;
  }

  tryToinitViewers() {
    if (!this.initialized) {
      this.initialized = true;
      this.initViewers();
    }
  }

  async initViewers(): Promise<void> {

  }

  handleUrlParams(params: Map<string, string>) : void {
  }
}
