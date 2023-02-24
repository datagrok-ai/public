import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {INglViewer, PROPS as pdbPROPS} from '../viewers/ngl-viewer';
import {Unsubscribable} from 'rxjs';

export class NglViewerApp {
  private readonly appFuncName: string;
  private pdb: string;
  private df: DG.DataFrame;

  constructor(appFuncName: string) {
    this.appFuncName = appFuncName;
  }

  async init(data?: { ligands: DG.DataFrame, macromolecule: string }): Promise<void> {
    if (data) {
      await this.setData(data.ligands, data.macromolecule);
    } else {
      const [ligandsDf, macromolecule] = await NglViewerApp.loadData();
      await this.setData(ligandsDf, macromolecule);
    }
  }

  /** Loads default data for {@link NglViewerApp} */
  static async loadData(): Promise<[DG.DataFrame, string]> {
    const sdfBytes: Uint8Array = await _package.files.readAsBytes('samples/1bdq.sdf');
    const ligandsDf: DG.DataFrame = (await grok.functions.call(
      'Chem:importSdf', {bytes: sdfBytes}))[0];
    const macromolecule: string = await _package.files.readAsText('samples/protease.pdb');
    return [ligandsDf, macromolecule];
  }

  async setData(ligandsDf: DG.DataFrame, macromolecule: string): Promise<void> {
    this.df = ligandsDf;
    this.pdb = macromolecule;

    await this.buildView();
  }

  // -- View --

  private view: DG.TableView;
  private viewSubs: Unsubscribable[] = [];

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.appFuncName}`;

    this.df.currentRowIdx = -1;

    const viewer: DG.JsViewer | INglViewer = (await this.df.plot.fromType('Ngl', {
      [pdbPROPS.pdb]: this.pdb,
      [pdbPROPS.ligandColumnName]: 'molecule',
    })) as DG.JsViewer | INglViewer;
    this.view.dockManager.dock(viewer as DG.JsViewer, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.4);

    this.viewSubs.push((viewer as INglViewer).onAfterBuildView.subscribe(() => {
      if (this.df.rowCount > 0)
        this.df.currentRowIdx = 0;
    }));
  }
}
