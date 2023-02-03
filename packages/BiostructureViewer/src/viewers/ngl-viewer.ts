import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';

import {TwinPviewer} from './twin-p-viewer';
import {Unsubscribable} from 'rxjs';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {_package} from '../package';
import {LoaderParameters} from 'NGL';

import * as NGL from 'NGL';

export interface INglViewer {
  get pdb(): string;

  set pdb(value: string);
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
      this.setData();
      break;
    }
  }

  // effective PDB value (to plot)
  private pdbStr: string | null = null;

  override onTableAttached(): void {
    super.onTableAttached();

    // -- Editors --
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.pdbTag).choices = ['', ...dfTagNameList];

    this.setData();
  }

  override detach(): void {
    super.detach();

    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }
  }

  // -- Data --

  setData(): void {
    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }

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
      this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private nglDiv?: HTMLDivElement;
  private stage?: NGL.Stage;

  private splashDiv?: HTMLDivElement;

  private viewSubs: Unsubscribable[] = [];

  private destroyView(): void {
    console.debug('BiostructureViewer: NglViewer.destroyView() ');
    if (this.pdbStr) {
      if (this.nglDiv && this.stage)
        this.stage.removeAllComponents();
    }

    for (const sub of this.viewSubs) sub.unsubscribe();

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }
  }

  private buildView(): void {
    console.debug('BiostructureViewer: NglViewer.buildView() ');
    if (this.pdbStr) {
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
      window.setTimeout(async () => {
        if (pdbStr) {
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

          this.viewSubs.push(this.dataFrame.onCurrentRowChanged
            .subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
        }
      }, 0 /* next event cycle */);
    } else {
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
          this.setData();
        }
      });
      const fileLink = ui.link('Open...', '', '', {
        onClick: (node) => {
          const k = 11;
          $(fileEl).trigger('click');
        }
      });
      this.splashDiv = ui.div([fileLink, fileEl],
        {style: {width: '100%', height: '100%', verticalAlign: 'middle', fontSize: 'larger'}});
      this.root.appendChild(this.splashDiv);
    }
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
    console.debug('BiostructureViewer: NglViewer.dataFrameOnCurrentRowChanged() ');

    const dataFrame: DG.DataFrame = this.dataFrame;
    const ligandColumnName: string = this.ligandColumnName;
    const stage: NGL.Stage = this.stage!;
    window.setTimeout(async () => {
      if (!ligandColumnName || dataFrame.currentRowIdx == -1) return;

      // remove all components but the first
      if (stage.compList.length > 1) {
        for (let compI = stage.compList.length - 1; compI > 0; compI--) {
          const comp = stage.compList[compI];
          stage.removeComponent(comp);
        }
      }

      // const ligandStr: string = await _package.files.readAsText('samples/1bdq.pdb');
      // const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
      // const ligandParams: LoaderParameters = {ext: 'sdf', compressed: false, binary: false, name: '<Ligand>'};

      const ligandMol: string = dataFrame.get(ligandColumnName, dataFrame.currentRowIdx);
      const ligandStr: string = ligandMol + '$$$$';
      const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
      const ligandParams: LoaderParameters = {ext: 'sdf', compressed: false, binary: false, name: '<Ligand>'};

      const loadRes = await stage.loadFile(ligandBlob, ligandParams);
      const repComp = stage.compList[1].addRepresentation(RepresentationType.BallAndStick, {});
    }, 0);
  }
}
