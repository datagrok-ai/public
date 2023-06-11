import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Observable, Subject, Unsubscribable} from 'rxjs';

import {INglViewer} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';

import {NglViewer, PROPS as pdbPROPS} from '../viewers/ngl-viewer';

import {_package} from '../package';

export class NglViewerApp {
  private readonly appFuncName: string;
  private pdb?: string;
  private df?: DG.DataFrame;

  private _onAfterBuildView: Subject<void> = new Subject<void>();

  public get onAfterBuildView(): Observable<void> { return this._onAfterBuildView; }

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

  private view?: DG.TableView;
  private viewSubs: Unsubscribable[] = [];

  async buildView(): Promise<void> {
    if (!this.df) throw new Error('df is not set');

    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.appFuncName}`;

    this.df.currentRowIdx = -1;

    const viewer: DG.Viewer & INglViewer = (await this.df.plot.fromType('NGL', {
      [pdbPROPS.pdb]: this.pdb,
      [pdbPROPS.ligandColumnName]: 'molecule',
    })) as DG.Viewer & INglViewer;
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.4);

    this.viewSubs.push(viewer.onAfterBuildView.subscribe(() => {
      this._onAfterBuildView.next();
    }));
  }
}
