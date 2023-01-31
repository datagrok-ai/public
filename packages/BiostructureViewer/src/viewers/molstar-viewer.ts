import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';

//@ts-ignore
import {Viewer} from 'molstar/build/viewer/molstar';
import {Unsubscribable} from 'rxjs';

// TODO: find out which extensions are needed.
/*const Extensions = {
    'cellpack': PluginSpec.Behavior(CellPack),
    'dnatco-confal-pyramids': PluginSpec.Behavior(DnatcoConfalPyramids),
    'pdbe-structure-quality-report': PluginSpec.Behavior(PDBeStructureQualityReport),
    'rcsb-assembly-symmetry': PluginSpec.Behavior(RCSBAssemblySymmetry),
    'rcsb-validation-report': PluginSpec.Behavior(RCSBValidationReport),
    'anvil-membrane-orientation': PluginSpec.Behavior(ANVILMembraneOrientation),
    'g3d': PluginSpec.Behavior(G3DFormat),
    'model-export': PluginSpec.Behavior(ModelExport),
    'mp4-export': PluginSpec.Behavior(Mp4Export),
    'geo-export': PluginSpec.Behavior(GeometryExport),
    'ma-quality-assessment': PluginSpec.Behavior(MAQualityAssessment),
};*/

const DefaultViewerOptions = {
  extensions: [],
  layoutIsExpanded: false,
  layoutShowControls: false,
  layoutShowRemoteState: false,
  layoutShowSequence: false,
  layoutShowLog: false,
  layoutShowLeftPanel: false,
  collapseLeftPanel: true,
  collapseRightPanel: true,

  viewportShowExpand: false,
  viewportShowControls: false,

  pdbProvider: 'rcsb',
  emdbProvider: 'rcsb',
};

const enum PROPS_CATS {
  DATA = 'Data',
  STYLE = 'Style',
  LAYOUT = 'Layout'
}

export enum PROPS {
  // -- Data --
  pdb = 'pdb',
  pdbTag = 'pdbTag',
  ligandColumnName = 'ligandColumnName',
  pdbProvider = 'pdbProvider',
  emdbProvider = 'emdbProvider',

  // -- Style --
  representation = 'representation',

  // --Layout --
  layoutIsExpanded = 'layoutIsExpanded',
  layoutShowControls = 'layoutShowControls',
  layoutShowRemoteState = 'layoutShowRemoteState',
  layoutShowSequence = 'layoutShowSequence',
  layoutShowLog = 'layoutShowLog',
  layoutShowLeftPanel = 'layoutShowLeftPanel',
  collapseLeftPanel = 'collapseLeftPanel',
  collapseRightPanel = 'collapseRightPanel',
  viewportShowExpand = 'viewportShowExpand',
  viewportShowControls = 'viewportShowControls',
}

export interface IMolstarViewer {
  get pdb(): string;

  set pdb(value: string);
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

export class MolstarViewer extends DG.JsViewer implements IMolstarViewer {
  private viewed: boolean = false;

  // -- Data --
  [PROPS.pdb]: string;
  [PROPS.pdbTag]: string;
  [PROPS.ligandColumnName]: string;
  [PROPS.pdbProvider]: string;
  [PROPS.emdbProvider]: string;

  // --Style --
  [PROPS.representation]: string;

  // -- Layout --
  [PROPS.layoutIsExpanded]: boolean;
  [PROPS.layoutShowControls]: boolean;
  [PROPS.layoutShowRemoteState]: boolean;
  [PROPS.layoutShowSequence]: boolean;
  [PROPS.layoutShowLog]: boolean;
  [PROPS.layoutShowLeftPanel]: boolean;
  [PROPS.collapseLeftPanel]: boolean;
  [PROPS.collapseRightPanel]: boolean;
  [PROPS.viewportShowExpand]: boolean;
  [PROPS.viewportShowControls]: boolean;

  constructor() {
    super();

    // -- Data --
    this.pdb = this.string(PROPS.pdb, pdbDefault,
      {category: PROPS_CATS.DATA, userEditable: false});
    this.pdbTag = this.string(PROPS.pdbTag, null,
      {category: PROPS_CATS.DATA, choices: []});
    this.ligandColumnName = this.string(PROPS.ligandColumnName, null,
      {category: PROPS_CATS.DATA, semType: DG.SEMTYPE.MOLECULE});
    this.pdbProvider = this.string(PROPS.pdbProvider, 'rcsb',
      {category: PROPS_CATS.DATA});
    this.emdbProvider = this.string(PROPS.emdbProvider, 'rcsb',
      {category: PROPS_CATS.DATA});

    // -- Style --
    this.representation = this.string(PROPS.representation, RepresentationType.Cartoon,
      {category: PROPS_CATS.STYLE, choices: Object.values(RepresentationType)});

    // -- Layout --
    this.layoutIsExpanded = this.bool(PROPS.layoutIsExpanded, false,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowControls = this.bool(PROPS.layoutShowControls, false,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowRemoteState = this.bool(PROPS.layoutShowRemoteState, false,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowSequence = this.bool(PROPS.layoutShowSequence, false,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowLog = this.bool(PROPS.layoutShowLog, false,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowLeftPanel = this.bool(PROPS.layoutShowLeftPanel, false,
      {category: PROPS_CATS.LAYOUT});
    this.collapseLeftPanel = this.bool(PROPS.collapseLeftPanel, true,
      {category: PROPS_CATS.LAYOUT});
    this.collapseRightPanel = this.bool(PROPS.collapseRightPanel, true,
      {category: PROPS_CATS.LAYOUT});
    this.viewportShowExpand = this.bool(PROPS.viewportShowExpand, false,
      {category: PROPS_CATS.LAYOUT});
    this.viewportShowControls = this.bool(PROPS.viewportShowControls, false,
      {category: PROPS_CATS.LAYOUT});

    // --
  }

  override onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('BiostructureViewer: MolstarViewer.onPropertyChanged() property is null');
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

    // -- Ligand --
    // TODO: Ligand

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private splashDiv?: HTMLDivElement;

  private viewSubs: Unsubscribable[] = [];

  private destroyView(): void {
    console.debug('BiostrctureViewer: MolstarViewer.destroyView() ');
    if (this.pdbStr) {
      // Clear viewer
    }

    for (const sub of this.viewSubs) sub.unsubscribe();

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }
  }

  private buildView(): void {
    console.debug('BiostructureViewer: MolstarViewer.buildView() ');
    if (this.pdbStr) {
      // Fill in viewer
    }
  }

  private updateView(): void {

  }

  // -- Handle events --

  private dataFrameOnCurrentRowChanged(value: any): void {
    console.debug('BiostructureViewer: MolstarViewer.dataFrameOnCurrentRowChanged() ');
  }
}

export function initViewer(viewName: string = 'Mol*') {
  const view = grok.shell.newView(viewName);
  const viewerContainer = view.root;
  const viewer = Viewer.create(viewerContainer, DefaultViewerOptions);
  return viewer;
}

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} pdbID ID of an entity in PDB.
 */
export async function byId(pdbID: string) {
  initViewer()
    .then((v: any) => {
      v.loadPdb(pdbID);
    });
  //v.handleResize();
}

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} pdbData Data in PDB
 */
export async function byData(pdbData: string, name: string = 'Mol*') {
  initViewer(name)
    .then((v: any) => {
      v.loadStructureFromData(pdbData, 'pdb');
    });
  //v.handleResize();
}

