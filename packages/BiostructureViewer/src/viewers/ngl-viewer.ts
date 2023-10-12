import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';
import {Observable, Subject} from 'rxjs';

import {TwinPviewer} from './twin-p-viewer';
import {Unsubscribable} from 'rxjs';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {LoaderParameters, RepresentationParameters} from 'NGL';
import * as NGL from 'NGL';
import {intToHtml} from '@datagrok-libraries/utils/src/color';
import {
  INglViewer,
  NglProps,
  NglPropsDefault,
  RepresentationType,
} from '@datagrok-libraries/bio/src/viewers/ngl-gl-viewer';

import {_package} from '../package';

const enum PROPS_CATS {
  DATA = 'Data',
  STYLE = 'Style',
  BEHAVIOUR = 'Behaviour',
}

export const enum PROPS {
  // -- Data --
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
  private viewed: boolean = false;
  private _onAfterBuildView = new Subject<void>();

  public get onAfterBuildView(): Observable<void> { return this._onAfterBuildView; }

  // -- Data --
  [PROPS.pdb]: string;
  [PROPS.pdbTag]: string;
  [PROPS.ligandColumnName]: string;

  // -- Style --
  [PROPS.representation]: NGL.StructureRepresentationType;

  // -- Behavior --
  [PROPS.showSelectedRowsLigands]: boolean;
  [PROPS.showCurrentRowLigand]: boolean;
  [PROPS.showMouseOverRowLigand]: boolean;

  private twinPviewer: TwinPviewer;

  constructor() {
    super();
    this.helpUrl = '/help/visualize/viewers/ngl';

    // -- Data --
    this.pdb = this.string(PROPS.pdb, defaults.pdb,
      {category: PROPS_CATS.DATA, userEditable: false});
    this.pdbTag = this.string(PROPS.pdbTag, defaults.pdbTag,
      {category: PROPS_CATS.DATA, choices: []});
    this.ligandColumnName = this.string(PROPS.ligandColumnName, defaults.ligandColumnName,
      {category: PROPS_CATS.DATA, semType: DG.SEMTYPE.MOLECULE});

    // -- Style --
    this.representation = this.string(PROPS.representation, defaults.representation,
      {category: PROPS_CATS.STYLE, choices: Object.values(RepresentationType)}) as NGL.StructureRepresentationType;

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
      case PROPS.showCurrentRowLigand:
      case PROPS.showSelectedRowsLigands:
      case PROPS.showMouseOverRowLigand:
        this.rebuildViewLigands();
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
    this.detachPromise = this.detachPromise.then(async () => { // detach
      await this.viewPromise;
      if (this.setDataInProgress) return; // check setDataInProgress synced
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
      if (!this.setDataInProgress) this.setDataInProgress = true; else return; // check setDataInProgress synced
      try {
        if (this.viewed) {
          await this.destroyView('setData');
          this.viewed = false;
        }

        await this.detachPromise;
        // Wait whether this.dataFrame assigning has called detach() before continue set data and build view

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

        if (!this.viewed) {
          await this.buildView('setData').then(() => { this._onAfterBuildView.next(); });
          this.viewed = true;
        }
      } catch (err: any) {
        const errMsg = err instanceof Error ? err.message : err.toString();
        const stack = err instanceof Error ? err.stack : undefined;
        grok.shell.error(errMsg);
        _package.logger.error(errMsg, undefined, stack);
      } finally {
        this.setDataInProgress = false;
      }
    });
  }

  // -- View --

  private viewPromise: Promise<void> = Promise.resolve();
  private detachPromise: Promise<void> = Promise.resolve();
  private setDataInProgress: boolean = false;

  private nglDiv?: HTMLDivElement;
  private stage?: NGL.Stage;

  private splashDiv?: HTMLDivElement;

  private viewSubs: Unsubscribable[] = [];

  private async destroyView(purpose: string): Promise<void> {
    _package.logger.debug(`NglViewer.destroyView(purpose='${purpose}') `);
    if (this.pdbStr) {
      if (this.nglDiv && this.stage) {
        await this.destroyViewLigands();
        this.stage.removeAllComponents();
      }
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
        style: {width: '100%', height: '100%'},
      });
      this.root.appendChild(this.nglDiv);

      this.stage = new NGL.Stage(this.nglDiv);
    }

    const stage: NGL.Stage = this.stage!;
    const representation: NGL.StructureRepresentationType = this.representation;
    const pdbStr: string = this.pdbStr;
    const df: DG.DataFrame = this.dataFrame;

    const pdbBlob = new Blob([pdbStr], {type: 'text/plain'});
    await stage.loadFile(pdbBlob, {ext: 'pdb', compressed: false, binary: false, name: '<Name>'});

    //highlights in NGL
    /* eslint-disable camelcase, prefer-const */
    let scheme_buffer: [string, string][] = [];

    //TODO: remove - demo purpose only
    scheme_buffer.push(['#0069a7', `* and :A`]);
    scheme_buffer.push(['#f1532b', `* and :B`]);
    scheme_buffer.push(['green', `* and :R`]);
    scheme_buffer.push(['green', `* and :M`]);
    /* eslint-enable camelcase, prefer-const */

    const schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
    const _schemeObj = {color: schemeId};

    const _repComp = stage.compList[0].addRepresentation(representation, {});
    stage.compList[0].autoView();

    this.viewSubs.push(df.onSelectionChanged.subscribe(this.dataFrameOnSelectionChanged.bind(this)));
    this.viewSubs.push(df.onCurrentRowChanged.subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
    this.viewSubs.push(df.onMouseOverRowChanged.subscribe(this.dataFrameOnMouseOverRowChanged.bind(this)));

    await this.buildViewLigands();
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
  }

  private updateView() {
    if (this.nglDiv && this.stage) {
      this.stage.compList[0].removeAllRepresentations();
      this.stage.compList[0].addRepresentation(this.representation, {});
    }
  }

  // -- Handle events --

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
    _package.logger.debug('NglViewer.dataFrameOnCurrentRowChanged() ');
    if (this.showSelectedRowsLigands) this.rebuildViewLigands();
  }

  private dataFrameOnCurrentRowChanged(_value: any): void {
    _package.logger.debug('NglViewer.dataFrameOnCurrentRowChanged() ');
    if (this.showCurrentRowLigand) this.rebuildViewLigands();
  }

  private dataFrameOnMouseOverRowChanged(_value: any) {
    _package.logger.debug('NglViewer.dataFrameOnMouseOverRowChanged() ');

    if (this.showMouseOverRowLigand) this.rebuildViewLigands();
  }

  // -- Ligands routines --

  private ligands: LigandMap = {selected: [], current: null, hovered: null};

  private getLigandBlobOfRow(rowIdx: number): Blob {
    const ligandMol: string = this.dataFrame.get(this.ligandColumnName, rowIdx);
    const ligandStr: string = ligandMol + '$$$$';
    const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
    return ligandBlob;
  }

  private rebuildViewLigands(): void {
    this.viewPromise = this.viewPromise.then(async () => {
      if (this.viewed) {
        await this.destroyViewLigands();
        await this.buildViewLigands();
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
  private async buildViewLigands(): Promise<void> {
    if (!this.stage) throw new Error('The stage is not created'); // return; // There is not PDB data
    if (!this.ligandColumnName) return;

    this.ligands.selected = !this.showSelectedRowsLigands ? [] :
      wu(this.dataFrame.selection.getSelectedIndexes())
        .map((selRowIdx) => { return {rowIdx: selRowIdx, compIdx: null}; })
        .toArray();
    this.ligands.current = !this.showCurrentRowLigand ? null :
      this.dataFrame.currentRowIdx >= 0 ? {rowIdx: this.dataFrame.currentRowIdx, compIdx: null} : null;
    this.ligands.hovered = !this.showMouseOverRowLigand ? null :
      this.dataFrame.mouseOverRowIdx >= 0 ? {rowIdx: this.dataFrame.mouseOverRowIdx, compIdx: null} : null;

    const addLigandOnStage = async (rowIdx: number, color: DG.Color | null): Promise<number> => {
      const ligandBlob = this.getLigandBlobOfRow(rowIdx);
      const ligandParams: Partial<LoaderParameters> =
        {ext: 'sdf', compressed: false, binary: false, name: `<Ligand at row ${rowIdx}`};
      await this.stage!.loadFile(ligandBlob, ligandParams); // assume this all adds last compList
      const compIdx: number = this.stage!.compList.length - 1;

      const params: RepresentationParameters = {
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
  }
}
