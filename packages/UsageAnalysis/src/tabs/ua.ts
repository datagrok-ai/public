import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {zipSync} from 'fflate';

import {UaToolbox} from '../ua-toolbox';
import {UaQueryViewer} from '../viewers/abstract/ua-query-viewer';
import type {TimelineView} from './timeline';

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
    this.setRibbonPanels([[ui.iconFA('arrow-to-bottom', () => this.download(), 'Download data')]]);
  }

  static csvFiles(tables: {[name: string]: DG.DataFrame | null | undefined}): DG.FileInfo[] {
    return Object.entries(tables).filter(([_, df]) => df != null && df.rowCount > 0)
      .map(([name, df]) => DG.FileInfo.fromString(`${name}.csv`, df!.toCsv()));
  }

  async exportFiles(): Promise<DG.FileInfo[]> {
    return UaView.csvFiles(Object.fromEntries(this.viewers.filter((v) => v.errorDiv == null)
      .map((v) => [v.name, v.viewer?.dataFrame])));
  }

  async download(): Promise<void> {
    const progress = DG.TaskBarProgressIndicator.create('Preparing download...');
    const files = await this.exportFiles().finally(() => progress.close());
    const fileName = (s: string) => s.toLowerCase().replace(/[^a-z0-9.]+/g, '-');
    if (files.length === 0)
      grok.shell.warning('Nothing to download');
    else if (files.length === 1)
      DG.Utils.download(fileName(`${this.name}-${files[0].name}`), files[0].data as BlobPart);
    else
      DG.Utils.download(fileName(`${this.name}.zip`),
        zipSync(Object.fromEntries(files.map((f) => [fileName(f.name), f.data])), {level: 1}) as BlobPart);
  }

  setToolbox(uaToolbox: UaToolbox) {
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;
    this._resolveToolbox();
  }

  // Unblocks initViewers() for toolbox-independent hosts (e.g. the Release app), which reuse
  // toolbox-free tabs (Stress, Vulnerabilities, ...) without building the shared UaToolbox.
  markToolboxReady() {
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

  openTimeline(key: string, id: string): void {
    const handler = this.uaToolbox.viewHandler;
    handler.changeTab('Timeline');
    (handler.getCurrentView() as TimelineView).show(key, id);
  }

  switchRout(): void {}
}
