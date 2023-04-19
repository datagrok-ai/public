import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import $ from 'cash-dom';
import wu from 'wu';
import {Observable, Subject} from 'rxjs';

import {TwinPviewer} from './twin-p-viewer';
import {Unsubscribable} from 'rxjs';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {LoaderParameters} from 'NGL';
import * as NGL from 'NGL';

export interface INglViewer {
  get pdb(): string;

  set pdb(value: string);

  get onAfterBuildView(): Observable<void>;
}

const enum PROPS_CATS {
  DATA = 'Data',
  STYLE = 'Style',
}

export const enum PROPS {
  // -- Data --
  pdb = 'pdb',
  pdbTag = 'pdbTag',
  ligandColumnName = 'ligandColumnName',

  // -- Style --
  representation = 'representation',
}

const pdbDefault: string = '';

enum RepresentationType {
  Cartoon = 'cartoon',
  Backbone = 'backbone',
  BallAndStick = 'ball+stick',
  Licorice = 'licorice',
  Hyperball = 'hyperball',
  Surface = 'surface'
}

/**
 * https://nglviewer.org/ngl/api/manual/example/snippets.html
 * https://nglviewer.org/ngl/api/manual/usage/file-formats.html
 */
export class NglViewer extends DG.JsViewer implements INglViewer {
  private viewed: boolean = false;
  private _onAfterBuildView = new Subject<void>();

  public get onAfterBuildView(): Observable<void> { return this._onAfterBuildView; }

  // -- Data --
  [PROPS.pdb]: string;
  [PROPS.pdbTag]: string;
  [PROPS.ligandColumnName]: string;

  // -- Style --
  [PROPS.representation]: string;

  private twinPviewer: TwinPviewer;

  constructor() {
    super();

    // -- Data --
    this.pdb = this.string(PROPS.pdb, pdbDefault,
      {category: PROPS_CATS.DATA, userEditable: false});
    this.pdbTag = this.string(PROPS.pdbTag, null,
      {category: PROPS_CATS.DATA, choices: []});
    this.ligandColumnName = this.string(PROPS.ligandColumnName, null,
      {category: PROPS_CATS.DATA, semType: DG.SEMTYPE.MOLECULE});

    // -- Style --
    this.representation = this.string(PROPS.representation, RepresentationType.Cartoon,
      {category: PROPS_CATS.STYLE, choices: Object.values(RepresentationType)});

    // --
    this.root.style.textAlign = 'center';
    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
  }

  override onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('BiostructureViewer: NglViewer.onPropertyChanged() property is null');
      return;
    }

    switch (property.name) {
    case PROPS.representation:
      this.updateView();
      break;
    }

    switch (property.name) {
    case PROPS.pdb:
    case PROPS.pdbTag:
      this.setData('onPropertyChanged');
      break;
    }
  }

  // effective PDB value (to plot)
  private pdbStr: string | null = null;

  override onTableAttached(): void {
    const superOnTableAttached = super.onTableAttached.bind(this);

    // -- Props editors --
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.pdbTag).choices = ['', ...dfTagNameList];

    this.viewPromise = this.viewPromise.then(async () => { // onTableAttached
      superOnTableAttached();
      await this.setData('onTableAttached');
    });
  }

  override detach(): void {
    const superDetach = super.detach.bind(this);
    this.viewPromise = this.viewPromise.then(async () => { // detach
      if (this.viewed) {
        await this.destroyView('detach');
        this.viewed = false;
      }
      superDetach();
    });
  }

  // -- Data --

  setData(purpose: string): void {
    _package.logger.debug(`NglViewer.setData(purpose='${purpose}') `);

    this.viewPromise = this.viewPromise.then(async () => { // setData
      if (this.viewed) {
        await this.destroyView('setData');
        this.viewed = false;
      }
    });

    // -- PDB data --
    let pdbTag: string = pdbTAGS.PDB;
    if (this.pdbTag) pdbTag = this.pdbTag;
    this.pdbStr = this.dataFrame.getTag(pdbTag);
    if (this.pdb && this.pdb != pdbDefault) this.pdbStr = this.pdb;

    // -- Ligand --
    if (!this.ligandColumnName) {
      const molCol: DG.Column | null = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      if (molCol)
        this.ligandColumnName = molCol.name;
    }

    this.viewPromise = this.viewPromise.then(async () => {
      if (!this.viewed) {
        await this.buildView('setData').then(() => { this._onAfterBuildView.next(); });
        this.viewed = true;
      }
    });
  }

  // -- View --

  private viewPromise: Promise<void> = Promise.resolve();
  private nglDiv?: HTMLDivElement;
  private stage?: NGL.Stage;

  private splashDiv?: HTMLDivElement;

  private viewSubs: Unsubscribable[] = [];

  private async destroyView(purpose: string): Promise<void> {
    _package.logger.debug(`NglViewer.destroyView(purpose='${purpose}') `);
    if (this.pdbStr) {
      if (this.nglDiv && this.stage)
        this.stage.removeAllComponents();
    }

    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }
  }

  private async buildView(purpose: string): Promise<void> {
    _package.logger.debug(`NglViewer.buildView(purpose='${purpose}') `);
    if (this.pdbStr)
      await this.buildViewWithPdb();
    else
      await this.buildViewWithoutPdb();
  }

  private async buildViewWithPdb() {
    if (!this.pdbStr) throw new Error('NglViewer.buildViewWithPdb() pdbStr is empty');

    if (!this.nglDiv) {
      this.nglDiv = ui.div([], {
        classes: 'd4-ngl-viewer',
        style: {width: '100%', height: '100%'}
      });
      this.root.appendChild(this.nglDiv);

      this.stage = new NGL.Stage(this.nglDiv);
    }

    const stage: NGL.Stage = this.stage!;
    const representation: string = this.representation;
    const pdbStr: string = this.pdbStr;
    const df: DG.DataFrame = this.dataFrame;

    const pdbBlob = new Blob([pdbStr], {type: 'text/plain'});
    await stage.loadFile(pdbBlob, {ext: 'pdb', compressed: false, binary: false, name: '<Name>'});

    //highlights in NGL
    // eslint-disable-next-line camelcase, prefer-const
    let scheme_buffer: string[][] = [];

    //TODO: remove - demo purpose only
    scheme_buffer.push(['#0069a7', `* and :A`]);
    scheme_buffer.push(['#f1532b', `* and :B`]);
    scheme_buffer.push(['green', `* and :R`]);
    scheme_buffer.push(['green', `* and :M`]);

    const schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
    const schemeObj = {color: schemeId};

    const repComp = stage.compList[0].addRepresentation(representation, {});
    stage.compList[0].autoView();

    // this.viewSubs.push(df.onCurrentRowChanged
    //   .subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
    this.viewSubs.push(df.onSelectionChanged
      .subscribe(this.dataFrameOnSelectionChanged.bind(this)));
  }

  private async buildViewWithoutPdb() {
    // preventing recreate nglDiv once again because of GL nature
    if (this.nglDiv) {
      this.stage!.dispose();
      delete this.stage;
      $(this.nglDiv).empty();
      delete this.nglDiv;
    }

    const fileEl: HTMLInputElement = ui.element('input');
    fileEl.type = 'file';
    fileEl.style.display = 'none';
    fileEl.addEventListener('change', async (event) => {
      const k = 11;
      if (fileEl.files != null && fileEl.files.length == 1) {
        const pdbStr: string = await fileEl.files[0]!.text();
        this.pdb = pdbStr;
        this.setData('onFileElChange');
      }
    });
    const fileLink = ui.link('Open...', '', '', {
      processNode: (node: HTMLElement) => {
        const k = 11;
      },
      // @ts-ignore // ui.link argument options.onClick: (node: HTMLElement) => void
      onClick: (event: PointerEvent) => {
        event.preventDefault();
        $(fileEl).trigger('click');
      }
    });
    this.splashDiv = ui.div([fileLink, fileEl],
      {style: {width: '100%', height: '100%', verticalAlign: 'middle', fontSize: 'larger'}});
    this.root.appendChild(this.splashDiv);
  }

  private updateView() {
    if (this.nglDiv && this.stage) {
      this.stage.compList[0].removeAllRepresentations();
      this.stage.compList[0].addRepresentation(this.representation, {});
    }
  }

  // -- Handle events --

  private rootOnSizeChanged(value: any): void {
    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;

    if (this.nglDiv && this.stage) {
      this.nglDiv.style.width = `${cw}px`;
      this.nglDiv.style.height = `${ch}px`;
      this.stage.handleResize();
    }

    if (this.splashDiv) {
      this.splashDiv.style.width = `${cw}px`;
      this.splashDiv.style.height = `${ch}px`;
      this.splashDiv.style.lineHeight = `${ch}px`;
    }
  }

  private dataFrameOnCurrentRowChanged(value: any): void {
    _package.logger.debug('NglViewer.dataFrameOnCurrentRowChanged() ');

    const dataFrame: DG.DataFrame = this.dataFrame;
    const ligandColumnName: string = this.ligandColumnName;
    const stage: NGL.Stage = this.stage!;
    this.viewPromise = this.viewPromise.then(async () => {
      await this.ligandsClear();

      if (!ligandColumnName || dataFrame.currentRowIdx == -1) return;

      // const ligandStr: string = await _package.files.readAsText('samples/1bdq.pdb');
      // const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
      // const ligandParams: LoaderParameters = {ext: 'sdf', compressed: false, binary: false, name: '<Ligand>'};

      const ligandBlob: Blob = this.getLigandBlob(dataFrame.currentRowIdx);
      const ligandParams: Partial<LoaderParameters> = {
        ext: 'sdf', compressed: false, binary: false, name: `<Ligand at row ${dataFrame.currentRowIdx}>`
      };

      const loadRes = await stage.loadFile(ligandBlob, ligandParams);
      const repComp = stage.compList[1].addRepresentation(RepresentationType.BallAndStick, {});
    });
  }

  private dataFrameOnSelectionChanged(value: any): void {
    _package.logger.debug('NglViewer.dataFrameOnCurrentRowChanged() ');

    const dataFrame: DG.DataFrame = this.dataFrame;
    const ligandColumnName: string = this.ligandColumnName;
    const stage: NGL.Stage = this.stage!;
    this.viewPromise = this.viewPromise.then(async () => {
      await this.ligandsClear();

      if (!ligandColumnName || !dataFrame.selection.anyTrue) return;

      const selIdxList: Int32Array = dataFrame.selection.getSelectedIndexes();
      for (let selI: number = 0; selI < selIdxList.length; selI++) {
        const selIdx: number = selIdxList[selI];
        const ligandBlob: Blob = this.getLigandBlob(selIdx);
        const ligandParams: Partial<LoaderParameters> = {
          ext: 'sdf',
          compressed: false,
          binary: false,
          name: `<Ligand at row ${selIdx}>`
        };

        const loadRes = await stage.loadFile(ligandBlob, ligandParams);
        const repComp = stage.compList[selI + 1].addRepresentation(RepresentationType.BallAndStick, {});
      }
    });
  }

  // -- Ligands routines --

  private getLigandBlob(rowIdx: number): Blob {
    const ligandMol: string = this.dataFrame.get(this.ligandColumnName, rowIdx);
    const ligandStr: string = ligandMol + '$$$$';
    const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
    return ligandBlob;
  }

  private async ligandsClear(): Promise<void> {
    if (!this.stage) return;

    // remove all components but the first (zero index)
    if (this.stage.compList.length > 1) {
      for (let compI = this.stage.compList.length - 1; compI > 0; compI--) {
        const comp = this.stage.compList[compI];
        this.stage.removeComponent(comp);
      }
    }
  }
}
