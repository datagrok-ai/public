import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {TwinPviewer} from './twin-p-viewer';
import {Unsubscribable} from 'rxjs';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';

//@ts-ignore
//import * as NGL from 'NGL';

export interface INglViewer {
  get pdb(): string;

  set pdb(value: string);
}

const enum PROPS_CATS {
  DATA = 'Data',
  APPEARANCE = 'Appearance',
}

export const enum PROPS {
  // -- Data --
  pdb = 'pdb',
  pdbTag = 'pdbTag',
  ligandColumnName = 'ligandColumnName',

  // -- Appearance --
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

export class NglViewer extends DG.JsViewer implements INglViewer {
  private viewed: boolean = false;

  // -- Data --
  [PROPS.pdb]: string;
  [PROPS.pdbTag]: string;
  [PROPS.ligandColumnName]: string;

  // -- Appearance --
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

    // -- Appearance --
    this.representation = this.string(PROPS.representation, RepresentationType.Cartoon,
      {category: PROPS_CATS.APPEARANCE, choices: Object.values(RepresentationType)});

    // --
    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
  }

  override onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);

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

  async onTableAttached() {
    super.onTableAttached();

    // -- Editors --
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.pdbTag).choices = ['', ...dfTagNameList];

    this.setData();
  }

  // -- Data --

  private setData(): void {
    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }

    // -- PDB data --
    let pdbTag: string = pdbTAGS.PDB;
    if (this.pdbTag) pdbTag = this.pdbTag;
    this.pdbStr = this.dataFrame.getTag(pdbTag);
    if (this.pdb && this.pdb != pdbDefault) this.pdbStr = this.pdb;

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

  private nglDiv: HTMLDivElement;

  //@ts-ignore
  private stage: NGL.Stage;

  private viewSubs: Unsubscribable[] = [];

  private destroyView(): void {
    console.debug('BiostructureViewer: NglViewer.destroyView() ');
    this.stage.removeAllComponents();

    for (const sub of this.viewSubs) sub.unsubscribe();
  }

  private buildView(): void {
    if (!this.nglDiv) {
      this.nglDiv = ui.div([], {
        classes: 'd4-ngl-viewer',
        style: {width: '100%', height: '100%'}
      });
      this.root.appendChild(this.nglDiv);

      //@ts-ignore
      this.stage = new NGL.Stage(this.nglDiv);
    }

    //@ts-ignore
    const stage: NGL.Stage = this.stage;
    const representation: string = this.representation;
    const pdbStr: string = this.pdbStr;
    window.setTimeout(async () => {
      if (pdbStr) {
        const pdbBlob = new Blob([pdbStr], {type: 'text/plain'});
        await stage.loadFile(pdbBlob, {ext: 'pdb', compressed: false, binary: false, name: '<Name>'});

        //highlights in NGL
        let scheme_buffer: string[][] = [];

        //TODO: remove - demo purpose only
        scheme_buffer.push(['#0069a7', `* and :A`]);
        scheme_buffer.push(['#f1532b', `* and :B`]);
        scheme_buffer.push(['green', `* and :R`]);
        scheme_buffer.push(['green', `* and :M`]);

        //@ts-ignore
        const schemeId = NGL.ColormakerRegistry.addSelectionScheme(scheme_buffer);
        const schemeObj = {color: schemeId};

        const repComp = stage.compList[0].addRepresentation(representation, schemeObj);
        stage.compList[0].autoView();

        this.viewSubs.push(this.dataFrame.onCurrentRowChanged.subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
      }
    }, 0 /* next event cycle */);
  }

  private updateView() {
    this.stage.compList[0].representation = this.representation;
  }

  // -- Handle events --

  private rootOnSizeChanged(value: any): void {
    if (!this.nglDiv || !this.stage) return;

    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;

    this.nglDiv.style.width = `${cw}px`;
    this.nglDiv.style.height = `${ch}px`;
    this.stage.handleResize();
  }

  private dataFrameOnCurrentRowChanged(value: any): void {
    console.debug('BiostructureViewer: NglViewer.dataFrameOnCurrentRowChanged() ');

    if (!this.ligandColumnName || this.dataFrame.currentRowIdx == -1) return;
    const ligandMol: string = this.dataFrame.get(this.ligandColumnName, this.dataFrame.currentRowIdx);
    const k = 11;
  }
}


