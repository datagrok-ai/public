/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MolTrackEms} from './ems/app';
import {DemoData} from './ems/demo-data';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.init()
  static async init(): Promise<void> {
    try {
      await grok.functions.call('Chem:initChemAutostart');
    } catch (e) {}
  }

  // loading the bundle registers MolTrack's table configuration, so a cold /domains link opens dressed
  @grok.decorators.autostart({meta: {autostartImmediate: 'true'}})
  static emsAutostart(): void {}

  @grok.decorators.app({
    name: 'MolTrack',
    description: 'Registers and searches compounds, batches and assay data',
    meta: {browsePath: 'Chem', icon: 'images/moltrack.png'},
  })
  static async molTrackApp(
    @grok.decorators.param({options: {metaUrl: true, optional: true, description: 'View to open: Compounds, Batches, Assays, AssayRuns, AssayResults or Additions'}}) path: string,
  ): Promise<DG.ViewBase> {
    return await MolTrackEms.open(path);
  }

  @grok.decorators.appTreeBrowser({app: 'MolTrack'})
  static async molTrackAppTreeBrowser(
    @grok.decorators.param({type: 'dynamic'}) appNode: DG.TreeViewGroup,
  ) {
    MolTrackEms.treeBrowser(appNode);
  }

  @grok.decorators.func({
    name: 'loadEmsDemoData',
    description: 'Loads the mol-track demo set into the moltrack domain schema where its tables are empty',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async loadEmsDemoData(): Promise<string> {
    return await DemoData.load();
  }
}
