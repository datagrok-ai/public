import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

// import {Viewer as rcsbViewer, ViewerProps as rcsbViewerProps} from '@rcsb/rcsb-molstar/src/viewer';
import {Viewer as RcsbViewer, ViewerProps as RcsbViewerProps} from '@rcsb/rcsb-molstar/build/src/viewer';

import $ from 'cash-dom';
import wu from 'wu';

//@ts-ignore
import {Viewer} from 'molstar/build/viewer/molstar';
import {Unsubscribable} from 'rxjs';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {ColorNames} from 'molstar/lib/mol-util/color/names';
import {ModelUrlProvider} from '@rcsb/rcsb-molstar/build/src/viewer/types';
import {_package} from '../package';
import {PluginContext} from 'molstar/lib/mol-plugin/context';

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
  LAYOUT = 'Layout',
  CONTROLS = 'Controls',
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

  // -- Controls --
  showImportControls = 'showImportControls',
}

const DefaultViewerProps: Partial<RcsbViewerProps> = {
  showImportControls: false,
  showSessionControls: false,
  showStructureSourceControls: true,
  // showMeasurementsControls: true,
  // showStrucmotifSubmitControls: true,
  showSuperpositionControls: true,
  // showQuickStylesControls: false,
  // showStructureComponentControls: true,
  // showVolumeStreamingControls: true,
  // showAssemblySymmetryControls: true,
  showValidationReportControls: true,

  showMembraneOrientationPreset: false,
  // showNakbColorTheme: false,
  /**
   * Needed when running outside of sierra. If set to true, the strucmotif UI will use an absolute URL to sierra-prod.
   * Otherwise, the link will be relative on the current host.
   */
  detachedFromSierra: false,
  modelUrlProviders: [
    (pdbId: string) => ({
      url: `https://models.rcsb.org/${pdbId.toLowerCase()}.bcif`,
      format: 'mmcif',
      isBinary: true
    }),
    (pdbId: string) => ({
      url: `https://files.rcsb.org/download/${pdbId.toLowerCase()}.cif`,
      format: 'mmcif',
      isBinary: false
    })
  ] as ModelUrlProvider[],

  layoutIsExpanded: false,
  layoutShowControls: true,
  layoutShowSequence: true,
  layoutShowLog: false,

  viewportShowExpand: true,
  viewportShowSelectionMode: true,
  volumeStreamingServer: 'https://maps.rcsb.org/',

  backgroundColor: ColorNames.white,
  showWelcomeToast: true
};

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

  // -- Controls --
  [PROPS.showImportControls]: boolean;

  // propsEngine = new class {
  //
  // }(this);

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

    // -- Controls --
    this.showImportControls = this.bool(PROPS.showImportControls, DefaultViewerProps.showImportControls,
      {category: PROPS_CATS.CONTROLS});

    this.subs.push(
      ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
  }

  override onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('BiostructureViewer: MolstarViewer.onPropertyChanged() property is null');
      return;
    }

    // const setPropsLayout = (props: Partial<RcsbViewerProps>) => {
    //   Object.assign(this.viewerLayoutState, props);
    //   if (this.viewer) {
    //     const viewerPlugin: PluginContext = this.viewer.getPlugin();
    //     viewerPlugin.layout.setProps();
    //
    //     // for (const propName in props) {
    //     //   const propValue: any = (props as { [propName: string]: any })[propName];
    //     //   this.viewer.customState[propName] = propValue;
    //     // }
    //   }
    // };

    switch (property.name) {
    case PROPS.representation:
      this.updateView();
      break;
    case PROPS.showImportControls:
    case PROPS.layoutIsExpanded:
    // case PROPS.layoutShowControls:
    //   const propValue: any = this[property.name];
    //   this.propsEngine.apply(property.name, propValue);
    //   _package.logger.debug('PhylocanvasGlViewer.onPropertyChanged() ' +
    //     `${property.name} = ${propValue.toString()} .`);
    //   break;
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
    // TODO: Ligand

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private viewerDiv?: HTMLDivElement;
  private viewer?: RcsbViewer;

  /** Container to store prop values while {@link viewer} is not created yet */
  private viewerProps: Partial<RcsbViewerProps> = {};

  private splashDiv?: HTMLDivElement;

  private viewSubs: Unsubscribable[] = [];

  private destroyView(): void {
    _package.logger.debug('MolstarViewer.destroyView() ');
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
    _package.logger.debug('MolstarViewer.buildView() ');
    if (this.pdbStr) {
      // Fill in viewer
      if (!this.viewerDiv) {
        this.viewerDiv = ui.div([], {
          classes: 'd4-molstar-viewer',
          style: {width: '100%', height: '100%'},
        });
        this.root.appendChild(this.viewerDiv);

        const props: Partial<RcsbViewerProps> = {};

        Object.assign(props, this.viewerProps);
        this.viewer = new RcsbViewer(this.viewerDiv, props);
      }
    }
  }

  private updateView(): void {

  }

  private calcSize(): void {
    if (!this.viewer) return;

    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    _package.logger.debug('MolstarViewer.calcSize( ${cw.toString()} x ${ch.toString()} )');

    if (this.viewerDiv) {
      this.viewerDiv.style.width = `${cw}px`;
      this.viewerDiv.style.height = `${ch}px`;

      this.viewer.handleResize();
    }
  }

  // -- Handle events --

  private rootOnSizeChanged(value: any): void {
    _package.logger.debug('MolstarViewer.rootOnSizeChanged() ');
    this.calcSize();
  }

  private dataFrameOnCurrentRowChanged(value: any): void {
    _package.logger.debug('MolstarViewer.dataFrameOnCurrentRowChanged() ');
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

