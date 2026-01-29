import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';
import {Observable, Subject, Unsubscribable} from 'rxjs';
import * as ngl from 'NGL';

import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {intToHtml} from '@datagrok-libraries/utils/src/color';
import {
  INglViewer,
  NglProps,
  NglPropsDefault,
  RepresentationType,
} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';
import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {testEvent} from '@datagrok-libraries/test/src/test';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';

import {awaitNgl} from './ngl-viewer-utils';
import {TwinPviewer} from './twin-p-viewer';

import {_package} from '../package';


const enum PROPS_CATS {
  DATA = 'Data',
  STYLE = 'Style',
  BEHAVIOUR = 'Behaviour',
}

export const enum PROPS {
  // -- Data --
  dataJson = 'dataJson',
  pdb = 'pdb',
  pdbTag = 'pdbTag',
  ligandColumnName = 'ligandColumnName',

  // -- Style --
  representation = 'representation',

  // -- Behaviour --
  showSelectedRowsLigands = 'showSelectedRowsLigands',
  showCurrentRowLigand = 'showCurrentRowLigand',
  showMouseOverRowLigand = 'showMouseOverRowLigand',
}

const pdbDefault: string = '';
const defaults: NglProps = NglPropsDefault;

export type LigandMapItem = { rowIdx: number, compIdx: number | null };

export type LigandMap = { selected: LigandMapItem[], current: LigandMapItem | null, hovered: LigandMapItem | null };

/**
 * https://nglviewer.org/ngl/api/manual/example/snippets.html
 * https://nglviewer.org/ngl/api/manual/usage/file-formats.html
 */
export class NglViewer extends DG.JsViewer implements INglViewer {
  private readonly logger: ILogger;
  private viewed: boolean = false;
  private _onAfterBuildView = new Subject<void>();

  public get onAfterBuildView(): Observable<void> { return this._onAfterBuildView; }

  // -- Data --
  [PROPS.dataJson]: string;
  [PROPS.pdb]: string;
  [PROPS.pdbTag]: string;
  [PROPS.ligandColumnName]: string;

  // -- Style --
  [PROPS.representation]: ngl.StructureRepresentationType;

  // -- Behavior --
  [PROPS.showSelectedRowsLigands]: boolean;
  [PROPS.showCurrentRowLigand]: boolean;
  [PROPS.showMouseOverRowLigand]: boolean;

  private twinPviewer: TwinPviewer;

  constructor() {
    super();

    this.helpUrl = '/help/visualize/viewers/ngl';

    // -- Data --
    this.dataJson = this.string(PROPS.dataJson, defaults.dataJson, {
      category: PROPS_CATS.DATA, userEditable: false,
      description: 'JSON encoded object of BiostructureData type with data value Base64 encoded data',
    });
    this.pdb = this.string(PROPS.pdb, defaults.pdb,
      {category: PROPS_CATS.DATA, userEditable: false});
    this.pdbTag = this.string(PROPS.pdbTag, defaults.pdbTag,
      {category: PROPS_CATS.DATA, choices: []});
    this.ligandColumnName = this.string(PROPS.ligandColumnName, defaults.ligandColumnName,
      {category: PROPS_CATS.DATA, semType: DG.SEMTYPE.MOLECULE});

    // -- Style --
    this.representation = this.string(PROPS.representation, defaults.representation,
      {category: PROPS_CATS.STYLE, choices: Object.values(RepresentationType)}) as ngl.StructureRepresentationType;

    // -- Behaviour --
    this.showSelectedRowsLigands = this.bool(PROPS.showSelectedRowsLigands, defaults.showSelectedRowsLigands,
      {category: PROPS_CATS.BEHAVIOUR});
    this.showCurrentRowLigand = this.bool(PROPS.showCurrentRowLigand, defaults.showCurrentRowLigand,
      {category: PROPS_CATS.BEHAVIOUR});
    this.showMouseOverRowLigand = this.bool(PROPS.showMouseOverRowLigand, defaults.showMouseOverRowLigand,
      {category: PROPS_CATS.BEHAVIOUR});

    // --
    this.root.style.textAlign = 'center';
    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

    this.logger = _package.logger;
    this.viewSyncer = new PromiseSyncer(this.logger);
  }

  private static viewerCounter: number = -1;
  private readonly viewerId: number = ++NglViewer.viewerCounter;

  private viewerToLog(): string { return `NglViewer<${this.viewerId}>`; }

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
      case PROPS.showCurrentRowLigand:
      case PROPS.showSelectedRowsLigands:
      case PROPS.showMouseOverRowLigand:
        this.rebuildViewLigands(`onPropertyChanged( name = ${property.name} )`);
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
  private dataEff: BiostructureData | null = null;

  override onTableAttached(): void {
    const superOnTableAttached = super.onTableAttached.bind(this);

    // -- Props editors --
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.pdbTag).choices = ['', ...dfTagNameList];

    superOnTableAttached();
    this.setData('onTableAttached');
  }

  override detach(): void {
    const callLog = `detach()`;
    const superDetach = super.detach.bind(this);
    this.viewSyncer.sync(callLog, async () => { // detach
      if (this.setDataInProgress) return; // check setDataInProgress synced
      if (this.viewed) {
        await this.destroyView(callLog);
        this.viewed = false;
      }
      superDetach();
    });
  }

  // -- Data --

  private _setDataCallCounter = -1;

  setData(caller: string): void {
    const callId = ++this._setDataCallCounter;
    const callLog = `setData( '${caller}', callId = ${callId} )`;
    this.viewSyncer.sync(callLog, async () => { // setData
      if (!this.setDataInProgress) this.setDataInProgress = true; else return; // check setDataInProgress synced
      try {
        if (this.viewed) {
          await this.destroyView(callLog);
          this.viewed = false;
        }
        // Wait whether this.dataFrame assigning has called detach() before continue set data and build view

        // -- PDB data --
        this.dataEff = null;
        let pdbTagName: string = pdbTAGS.PDB;
        let pdb: string | null = null;
        if (this.pdbTag) pdbTagName = this.pdbTag;
        if (this.dataFrame && this.dataFrame.tags.has(pdbTagName)) pdb = this.dataFrame.getTag(pdbTagName);
        if (this.pdb) pdb = this.pdb;
        if (pdb && pdb != pdbDefault)
          this.dataEff = {binary: false, ext: 'pdb', data: pdb!};
        if (this.dataJson && this.dataJson !== BiostructureDataJson.empty)
          this.dataEff = BiostructureDataJson.toData(this.dataJson);

        // -- Ligand --
        if (this.dataFrame && !this.ligandColumnName) {
          const molCol: DG.Column | null = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
          if (molCol)
            this.ligandColumnName = molCol.name;
        }

        if (!this.viewed) {
          await this.buildView(callLog);
          this._onAfterBuildView.next();
          this.viewed = true;
        }
      } catch (err: any) {
        const errMsg = err instanceof Error ? err.message : err.toString();
        const stack = err instanceof Error ? err.stack : undefined;
        grok.shell.error(errMsg);
        this.logger.error(errMsg, undefined, stack);
      } finally {
        this.setDataInProgress = false;
      }
    });
  }

  // -- View --

  private viewSyncer: PromiseSyncer;
  private setDataInProgress: boolean = false;

  private nglDiv?: HTMLDivElement;
  public stage?: ngl.Stage;

  private splashDiv?: HTMLDivElement;

  private viewSubs: Unsubscribable[] = [];

  private async destroyView(caller: string): Promise<void> {
    const callLog = `destroyView( ${caller} )`;
    const callLogPrefix = `${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${callLogPrefix}, start`);
    if (this.nglDiv && this.stage) {
      await this.destroyViewLigands();
      this.stage.removeAllComponents();
    }

    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }

    if (this.nglDiv) {
      this.logger.debug(`${callLogPrefix}, stage removing`);
      this.stage!.dispose();
      delete this.stage;
      $(this.nglDiv).empty();
      this.nglDiv.remove();
      delete this.nglDiv;
      this.logger.debug(`${callLogPrefix}, stage removed`);
    }
    this.logger.debug(`${callLogPrefix}, end`);
  }

  private async buildView(caller: string): Promise<void> {
    const callLog = `buildView( ${caller} )`;
    const callLogPrefix = `${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${callLogPrefix}, start`);
    if (this.dataEff)
      await this.buildViewWithPdb(callLog);
    else
      await this.buildViewWithoutPdb(callLog);
    this.logger.debug(`${callLogPrefix}, end`);
  }

  private async buildViewWithPdb(caller: string) {
    const callLog = `buildViewWithPdb( ${caller} )`;
    const callLogPrefix = `${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${callLogPrefix}, start`);
    if (!this.dataEff) throw new Error(`${callLogPrefix}, ` + 'pdbStr is empty');

    if (!this.nglDiv) {
      this.nglDiv = ui.div([], {
        classes: 'd4-ngl-viewer',
        style: {width: '100%', height: '100%'},
      });
      this.root.appendChild(this.nglDiv);

      this.logger.debug(`${callLogPrefix}, stage creating`);
      this.stage = new ngl.Stage(this.nglDiv);
      await awaitNgl(this.stage, callLogPrefix);
      this.logger.debug(`${callLogPrefix}, stage created`);
    }

    const stage: ngl.Stage = this.stage!;
    const representation: ngl.StructureRepresentationType = this.representation;
    const df: DG.DataFrame = this.dataFrame;

    const dataVal: Blob = this.dataEff.binary ? new Blob([this.dataEff.data as BlobPart]) :
      new Blob([this.dataEff.data as string], {type: 'text/plain'});
    await stage.loadFile(dataVal, {
      ext: this.dataEff.ext, compressed: false, binary: this.dataEff.binary, name: '<Name>',
      defaultRepresentation: true
    });

    //highlights in NGL
    /* eslint-disable camelcase, prefer-const */
    let scheme_buffer: ngl.SelectionSchemeData[] = [];

    //TODO: remove - demo purpose only
    scheme_buffer.push(['#0069a7', `* and :A`, undefined]);
    scheme_buffer.push(['#f1532b', `* and :B`, undefined]);
    scheme_buffer.push(['green', `* and :R`, undefined]);
    scheme_buffer.push(['green', `* and :M`, undefined]);
    /* eslint-enable camelcase, prefer-const */

    const schemeId = ngl.ColormakerRegistry.addSelectionScheme(scheme_buffer);
    const _schemeObj = {color: schemeId};

    // const _repComp = stage.compList[0].addRepresentation(representation, {});
    stage.compList[0].autoView();

    if (df && this.ligandColumnName) {
      this.viewSubs.push(df.onSelectionChanged.subscribe(this.dataFrameOnSelectionChanged.bind(this)));
      this.viewSubs.push(df.onCurrentRowChanged.subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
      this.viewSubs.push(df.onMouseOverRowChanged.subscribe(this.dataFrameOnMouseOverRowChanged.bind(this)));
      await this.buildViewLigands(callLog);
    }
    this.logger.debug(`${callLogPrefix}, end`);
  }

  private async buildViewWithoutPdb(caller: string) {
    const callLog = `buildViewWithoutPdb( ${caller} )`;
    const callLogPrefix = `${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${callLogPrefix}, start`);
    // preventing recreate nglDiv once again because of GL nature
    if (this.nglDiv) {
      this.logger.debug(`${callLogPrefix}, stage removing`);
      this.stage!.dispose();
      delete this.stage;
      $(this.nglDiv).empty();
      delete this.nglDiv;
      this.logger.debug(`${callLogPrefix}, stage removed`);
    }

    const fileEl: HTMLInputElement = ui.element('input');
    fileEl.type = 'file';
    fileEl.style.display = 'none';
    fileEl.addEventListener('change', async (_event) => {
      if (fileEl.files != null && fileEl.files.length == 1) {
        const pdbStr: string = await fileEl.files[0]!.text();
        this.pdb = pdbStr;
        this.setData('onFileElChange'); // isolates async, errors get lost
      }
    });
    const fileLink = ui.link('Open...', '', '', {
      processNode: (_node: HTMLElement) => {
        const _k = 11;
      },
      // @ts-ignore // ui.link argument options.onClick: (node: HTMLElement) => void
      onClick: (event: PointerEvent) => {
        event.preventDefault();
        $(fileEl).trigger('click');
      },
    });
    this.splashDiv = ui.div([fileLink, fileEl],
      {style: {width: '100%', height: '100%', verticalAlign: 'middle', fontSize: 'larger'}});
    this.root.appendChild(this.splashDiv);

    this.logger.debug(`${callLogPrefix}, end`);
  }

  private updateView() {
    if (this.nglDiv && this.stage) {
      this.stage.compList[0].removeAllRepresentations();
      this.stage.compList[0].addRepresentation(this.representation, {});
    }
  }

  // -- Handle events --

  private handleError(err: any): void {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.error(errMsg);
    this.logger.error(errMsg, undefined, errStack);
  }

  private rootOnSizeChanged(_value: any): void {
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

  private dataFrameOnSelectionChanged(_value: any): void {
    this.logger.debug(`${this.viewerToLog()}.dataFrameOnCurrentRowChanged() `);
    if (this.showSelectedRowsLigands) this.rebuildViewLigands('dataFrameSelectionChanged()');
  }

  private dataFrameOnCurrentRowChanged(_value: any): void {
    this.logger.debug(`${this.viewerToLog()}.dataFrameOnCurrentRowChanged() `);
    if (this.showCurrentRowLigand) this.rebuildViewLigands('dataFrameOnCurrentRowChanged()');
  }

  private dataFrameOnMouseOverRowChanged(_value: any) {
    this.logger.debug(`${this.viewerToLog()}.dataFrameOnMouseOverRowChanged() `);

    if (this.showMouseOverRowLigand) this.rebuildViewLigands('dataFrameOnMouseOverRowChanged()');
  }

  // -- Ligands routines --

  private ligands: LigandMap = {selected: [], current: null, hovered: null};

  private getLigandBlobOfRow(rowIdx: number): Blob {
    const ligandMol: string = this.dataFrame.get(this.ligandColumnName, rowIdx);
    const ligandStr: string = ligandMol + '$$$$';
    const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
    return ligandBlob;
  }

  private rebuildViewLigands(caller: string): void {
    const callLog = `rebuildViewLigands( ${caller} )`;
    this.viewSyncer.sync(callLog, async () => {
      if (this.viewed) {
        await this.destroyViewLigands();
        await this.buildViewLigands(callLog);
      }
    });
  }

  private async destroyViewLigands(): Promise<void> {
    if (!this.stage) throw new Error('The stage is not created'); // return; // There is not PDB data
    if (!this.ligandColumnName) return;

    const allLigands: LigandMapItem[] = [
      ...this.ligands.selected,
      ...(this.ligands.current ? [this.ligands.current] : []),
      ...(this.ligands.hovered ? [this.ligands.hovered] : []),
    ];
    const desc = <T>(a: T, b: T): number => (a > b ? -1 : 1);
    for (const compIdx of allLigands.map((l) => l.compIdx!).sort(desc)) {
      const comp = this.stage.compList[compIdx];
      this.stage.removeComponent(comp);
    }
    for (const ligand of allLigands) ligand.compIdx = null; // unbind with this.stage.compList
  }

  /** Builds up ligands on the stage view */
  private async buildViewLigands(caller: string): Promise<void> {
    const callLog = `buildViewLigands( ${caller} )`;
    const callLogPrefix = `${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${callLogPrefix}, start`);
    try {
      if (!this.stage) throw new Error('The stage is not created'); // return; // There is not PDB data
      if (!this.dataFrame || !this.ligandColumnName) return;

      this.ligands.selected = !this.showSelectedRowsLigands ? [] :
        wu(this.dataFrame.selection.getSelectedIndexes())
          .take(25) // Limit selected ligands to display
          .map((selRowIdx) => { return {rowIdx: selRowIdx, compIdx: null}; })
          .toArray();
      this.ligands.current = !this.showCurrentRowLigand ? null :
        this.dataFrame.currentRowIdx >= 0 ? {rowIdx: this.dataFrame.currentRowIdx, compIdx: null} : null;
      this.ligands.hovered = !this.showMouseOverRowLigand ? null :
        this.dataFrame.mouseOverRowIdx >= 0 ? {rowIdx: this.dataFrame.mouseOverRowIdx, compIdx: null} : null;

      const addLigandOnStage = async (rowIdx: number, color: DG.Color | null): Promise<number> => {
        const ligandBlob = this.getLigandBlobOfRow(rowIdx);
        const ligandParams: Partial<ngl.LoaderParameters> =
          {ext: 'sdf', compressed: false, binary: false, name: `<Ligand at row ${rowIdx}`};
        await this.stage!.loadFile(ligandBlob, ligandParams); // assume this all adds last compList
        const compIdx: number = this.stage!.compList.length - 1;

        const params: Partial<ngl.RepresentationParameters> = {
          ...(color ? {color: intToHtml(color as number)} : {}),
        };

        this.stage!.compList[compIdx].addRepresentation(RepresentationType.BallAndStick, params);
        return compIdx;
      };

      const selCount = this.ligands.selected.length;
      for (const [selectedLigand, selI] of wu.enumerate(this.ligands.selected)) {
        const color =
          this.showCurrentRowLigand || this.showMouseOverRowLigand ?
            (selCount > 1 ? DG.Color.selectedRows : null) :
            (selCount > 1 ? DG.Color.scaleColor(selI, 0, selCount, 0.5) : null);

        selectedLigand.compIdx = await addLigandOnStage(selectedLigand.rowIdx, color);
      }
      if (this.ligands.current) {
        const color = this.showSelectedRowsLigands ? DG.Color.currentRow : null;

        this.ligands.current.compIdx = await addLigandOnStage(this.ligands.current.rowIdx, color);
      }
      if (this.ligands.hovered) {
        // TODO: color hovered ligand
        const color =
          this.showSelectedRowsLigands || this.showCurrentRowLigand ?
            DG.Color.mouseOverRows : null;

        this.ligands.hovered.compIdx = await addLigandOnStage(this.ligands.hovered.rowIdx, color);
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      throw new Error(`NglViewer.buildViewLigands() error: ${errMsg}`, {cause: err});
    } finally {
      this.logger.debug(`${callLogPrefix}, end`);
    }
  }

  // -- IRenderer --

  private _onRendered: Subject<void> = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    this.viewSyncer.sync('invalidate(${caller ? ` <- ${caller} ` : \'\'})', async () => {
      // update view / render
      this._onRendered.next();
    });
  }

  async awaitRendered(timeout: number | undefined = 5000): Promise<void> {
    await testEvent(this.onRendered, () => {}, () => {
      this.invalidate();
    }, timeout);

    // Rethrow stored syncer error (for test purposes)
    const viewErrors = this.viewSyncer.resetErrors();
    if (viewErrors.length > 0) throw viewErrors[0];
  }
}
