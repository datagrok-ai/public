import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {PROPS as pdbPROPS} from '../viewers/ngl-viewer';

export class NglViewerApp {
  private readonly appFuncName: string;
  private pdb: string;
  private df: DG.DataFrame;

  constructor(appFuncName: string) {
    this.appFuncName = appFuncName;
  }

  async init(data?:{ligands: DG.DataFrame, macromolecule: string}): Promise<void> {
    if(data) {
      this.pdb = data.macromolecule;
      this.setData(data.ligands);
    }
    else 
      await this.loadData();
  }

  async loadData(): Promise<void> {
    const sdfBytes: Uint8Array = await _package.files.readAsBytes('samples/1bdq.sdf');
    const df: DG.DataFrame = (await grok.functions.call(
      'Chem:importSdf', {bytes: sdfBytes}))[0];
    this.pdb = await _package.files.readAsText('samples/protease.pdb');
    await this.setData(df);
  }

  async setData(df: DG.DataFrame): Promise<void> {
    this.df = df;

    await this.buildView();
  }

  // -- View --

  private view: DG.TableView;

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.appFuncName}`;

    const viewer: DG.JsViewer = (await this.view.dataFrame.plot.fromType('NglViewer', {
      [pdbPROPS.pdb]: this.pdb,
      [pdbPROPS.ligandColumnName]: 'molecule',
    })) as DG.JsViewer;
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.4);
  }
}
