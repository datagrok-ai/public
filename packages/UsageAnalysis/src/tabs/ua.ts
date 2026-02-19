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
  uaToolbox!: UaToolbox;
  viewers: UaQueryViewer[] = [];
  initialized: boolean = false;
  systemId: string = '00000000-0000-0000-0000-000000000000';
  rout?: string;
  private _toolboxReady: Promise<void>;
  private _resolveToolbox!: () => void;

  constructor(uaToolbox?: UaToolbox) {
    super();
    this._toolboxReady = new Promise((resolve) => this._resolveToolbox = resolve);
    if (uaToolbox)
      this.setToolbox(uaToolbox);
    this.box = true;
  }

  setToolbox(uaToolbox: UaToolbox) {
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;
    this._resolveToolbox();
  }

  async tryToInitViewers(path?: string): Promise<void> {
    await this._toolboxReady;
    if (!this.initialized) {
      this.initialized = true;
      await this.initViewers(path);
      for (const viewer of this.viewers) {
        if (!viewer.activated) {
          viewer.activated = true;
          viewer.reloadViewer();
        }
      }
    }
  }

  async initViewers(path?: string): Promise<void> {}

  switchRout(): void {}
}
