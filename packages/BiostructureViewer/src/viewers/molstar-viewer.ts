import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

// import {Viewer as rcsbViewer, ViewerProps as rcsbViewerProps} from '@rcsb/rcsb-molstar/src/viewer';
import {Viewer as RcsbViewer, ViewerProps as RcsbViewerProps} from '@rcsb/rcsb-molstar/build/src/viewer';
import {PROPS as pdbPROPS} from './../viewers/ngl-viewer';


import {_package} from '../package';
import $ from 'cash-dom';
import wu from 'wu';

import {Unsubscribable} from 'rxjs';

import {
  BiostructureProps,
  BiostructurePropsDefault,
  IBiostructureViewer,
  PluginLayoutControlsDisplayType,
  RegionStateOptionsType,
  RepresentationType,
  SimpleRegionStateOptionsType
} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb';
import {ColorNames} from 'molstar/lib/mol-util/color/names';
import {PluginCommands} from 'molstar/lib/mol-plugin/commands';
import {ModelUrlProvider} from '@rcsb/rcsb-molstar/build/src/viewer/types';
import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {PluginLayoutControlsDisplay} from 'molstar/lib/mol-plugin/layout';
import {PdbResDataFrame} from '../utils/pdb-helper';
import {BuiltInTrajectoryFormat} from 'molstar/lib/mol-plugin-state/formats/trajectory';
import {intToHtmlA} from '@datagrok-libraries/utils/src/color';
import {Color as msColor} from 'molstar/lib/mol-util/color';
import {getPdbHelper, IPdbHelper, PdbResDataFrameType} from '@datagrok-libraries/bio/src/pdb/pdb-helper';


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
  layoutRegionStateLeft = 'layoutRegionStateLeft',
  layoutRegionStateTop = 'layoutRegionStateTop',
  layoutRegionStateRight = 'layoutRegionStateRight',
  layoutRegionStateBottom = 'layoutRegionStateBottom',
  layoutControlsDisplay = 'layoutControlsDisplay',

  layoutShowRemoteState = 'layoutShowRemoteState',
  layoutShowSequence = 'layoutShowSequence',
  layoutShowLog = 'layoutShowLog',
  layoutShowLeftPanel = 'layoutShowLeftPanel',
  collapseLeftPanel = 'collapseLeftPanel',
  collapseRightPanel = 'collapseRightPanel',
  viewportShowExpand = 'viewportShowExpand',
  viewportShowControls = 'viewportShowControls',

  // -- Controls --
  showWelcomeToast = 'showWelcomeToast',
  showImportControls = 'showImportControls',
}

// const MolstarViewerDefaultOptions: Partial<RcsbViewerProps> = {
//   // extensions: [],
//   layoutIsExpanded: false,
//   layoutShowControls: false,
//   // layoutControlsDisplay: 'outside',
//   // layoutShowRemoteState: false,
//   layoutShowSequence: false,
//   layoutShowLog: false,
//   // layoutShowLeftPanel: false,
//   // collapseLeftPanel: true,
//   // collapseRightPanel: true,
//   //
//   // viewportShowExpand: false,
//   // viewportShowControls: false,
//   //
//   // pdbProvider: 'rcsb',
//   // emdbProvider: 'rcsb',
//   showWelcomeToast: false
// };

const pdbDefault: string = '';
const defaults: BiostructureProps = BiostructurePropsDefault;

export class MolstarViewer extends DG.JsViewer implements IBiostructureViewer {
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
  [PROPS.layoutRegionStateLeft]: RegionStateOptionsType;
  [PROPS.layoutRegionStateTop]: SimpleRegionStateOptionsType;
  [PROPS.layoutRegionStateRight]: SimpleRegionStateOptionsType;
  [PROPS.layoutRegionStateBottom]: SimpleRegionStateOptionsType;
  [PROPS.layoutControlsDisplay]: string;

  [PROPS.layoutShowRemoteState]: boolean;
  [PROPS.layoutShowSequence]: boolean;
  [PROPS.layoutShowLog]: boolean;
  [PROPS.layoutShowLeftPanel]: boolean;
  [PROPS.collapseLeftPanel]: boolean;
  [PROPS.collapseRightPanel]: boolean;
  [PROPS.viewportShowExpand]: boolean;
  [PROPS.viewportShowControls]: boolean;

  // -- Controls --
  [PROPS.showWelcomeToast]: boolean;
  [PROPS.showImportControls]: boolean;

  // propsEngine = new class {
  //
  // }(this);

  constructor() {
    super();

    // -- Data --
    this.pdb = this.string(PROPS.pdb, pdbDefault,
      {category: PROPS_CATS.DATA, userEditable: false});
    this.pdbTag = this.string(PROPS.pdbTag, defaults.pdbTag,
      {category: PROPS_CATS.DATA, choices: []});
    this.ligandColumnName = this.string(PROPS.ligandColumnName, defaults.ligandColumnName,
      {category: PROPS_CATS.DATA, semType: DG.SEMTYPE.MOLECULE});
    this.pdbProvider = this.string(PROPS.pdbProvider, defaults.pdbProvider,
      {category: PROPS_CATS.DATA});
    this.emdbProvider = this.string(PROPS.emdbProvider, defaults.emdbProvider,
      {category: PROPS_CATS.DATA});

    // -- Style --
    this.representation = this.string(PROPS.representation, defaults.representation,
      {category: PROPS_CATS.STYLE, choices: Object.values(RepresentationType)});

    // -- Layout --
    this.layoutIsExpanded = this.bool(PROPS.layoutIsExpanded, defaults.layoutIsExpanded,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowControls = this.bool(PROPS.layoutShowControls, defaults.layoutShowControls,
      {category: PROPS_CATS.LAYOUT});
    this.layoutRegionStateLeft = this.string(PROPS.layoutRegionStateLeft, defaults.layoutRegionStateLeft,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(RegionStateOptionsType)}) as RegionStateOptionsType;
    this.layoutRegionStateTop = this.string(PROPS.layoutRegionStateTop, defaults.layoutRegionStateTop,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SimpleRegionStateOptionsType)}
    ) as SimpleRegionStateOptionsType;
    this.layoutRegionStateRight = this.string(PROPS.layoutRegionStateRight, defaults.layoutRegionStateRight,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SimpleRegionStateOptionsType)}
    ) as SimpleRegionStateOptionsType;
    this.layoutRegionStateBottom = this.string(PROPS.layoutRegionStateBottom, defaults.layoutRegionStateBottom,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SimpleRegionStateOptionsType)}
    ) as SimpleRegionStateOptionsType;
    this.layoutControlsDisplay = this.string(PROPS.layoutControlsDisplay, defaults.layoutControlsDisplay,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(PluginLayoutControlsDisplayType)});

    this.layoutShowRemoteState = this.bool(PROPS.layoutShowRemoteState, defaults.layoutShowRemoteState,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowSequence = this.bool(PROPS.layoutShowSequence, defaults.layoutShowSequence,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowLog = this.bool(PROPS.layoutShowLog, defaults.layoutShowLog,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowLeftPanel = this.bool(PROPS.layoutShowLeftPanel, defaults.layoutShowLeftPanel,
      {category: PROPS_CATS.LAYOUT});
    this.collapseLeftPanel = this.bool(PROPS.collapseLeftPanel, defaults.collapseLeftPanel,
      {category: PROPS_CATS.LAYOUT});
    this.collapseRightPanel = this.bool(PROPS.collapseRightPanel, defaults.collapseRightPanel,
      {category: PROPS_CATS.LAYOUT});
    this.viewportShowExpand = this.bool(PROPS.viewportShowExpand, defaults.viewportShowExpand,
      {category: PROPS_CATS.LAYOUT});
    this.viewportShowControls = this.bool(PROPS.viewportShowControls, defaults.viewportShowControls,
      {category: PROPS_CATS.LAYOUT});

    // -- Controls --
    this.showWelcomeToast = this.bool(PROPS.showWelcomeToast, defaults.showWelcomeToast,
      {category: PROPS_CATS.CONTROLS});
    this.showImportControls = this.bool(PROPS.showImportControls, defaults.showImportControls,
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
    const applyProperty = async (propName: string, value: any) => {
      if (!this.viewer) throw new Error('viewer does not exists');

      const plugin: PluginContext = this.viewer.plugin;
      switch (property.name) {
      case PROPS.layoutShowLog: {
        //plugin.layout.setProps({layoutShowLog: value;});

        //PluginCommands.Layout.Update(plugin, );
        const k = 11;
      }
        break;
      case PROPS.layoutShowControls: {
        // eslint-disable-next-line new-cap
        await PluginCommands.Layout.Update(plugin, {state: {showControls: this.layoutShowControls}});
        const k = 11;
      }
        break;
      case PROPS.layoutIsExpanded: {
        //eslint-disable-next-line new-cap
        await PluginCommands.Layout.Update(plugin, {state: {isExpanded: this.layoutIsExpanded}});
      }
        break;
      case PROPS.layoutRegionStateLeft:
      case PROPS.layoutRegionStateTop:
      case PROPS.layoutRegionStateRight:
      case PROPS.layoutRegionStateBottom: {
        //eslint-disable-next-line new-cap
        await PluginCommands.Layout.Update(plugin, {
          state: {
            regionState: {
              left: this.layoutRegionStateLeft,
              top: this.layoutRegionStateTop,
              right: this.layoutRegionStateRight,
              bottom: this.layoutRegionStateBottom
            }
          }
        });
      }
        break;
      case PROPS.layoutControlsDisplay: {
        //eslint-disable-next-line new-cap
        await PluginCommands.Layout.Update(plugin,
          {state: {controlsDisplay: this.layoutControlsDisplay as PluginLayoutControlsDisplay}});
      }
        break;
        // case PROPS.viewportShowExpand: {
        //   await PluginCommands.State.ToggleExpanded(plugin,
        //     {state: {isExpanded: this.viewportShowExpand}})
        // }
        //   break;
      }
    };

    switch (property.name) {
    case PROPS.representation:
      break;
    case PROPS.showImportControls:
      break;
    case PROPS.layoutIsExpanded:
      this.viewerProps[PROPS.layoutIsExpanded] = this.layoutIsExpanded;
      break;
    }

    const propName: string = property.name;
    const propValue: any = this.props.get(propName);
    if (this.viewer) {
      this.viewPromise = this.viewPromise.then(async () => {
        await applyProperty(propName, propValue);
      });
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
    _package.logger.debug('MolstarViewer.onTableAttached(), ');
    const superOnTableAttached = super.onTableAttached.bind(this);

    // -- Props editors --
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.pdbTag).choices = ['', ...dfTagNameList];

    superOnTableAttached();
    this.setData();
  }

  override detach(): void {
    if (this.setDataInProgress) return;
    _package.logger.debug('MolstarViewer.detach(), ');

    const superDetach = super.detach.bind(this);
    this.detachPromise = this.detachPromise.then(async () => { // detach
      if (this.viewed) {
        await this.destroyView('detach');
        this.viewed = false;
      }
      superDetach();
    }).catch((reason: any) => {
      grok.shell.error(reason.toString());
    });
  }

  // -- Data --

  setData(): void {
    if (!this.setDataInProgress) this.setDataInProgress = true; else return;
    _package.logger.debug(`MolstarViewer.setData() `);

    this.viewPromise = this.viewPromise.then(async () => { // setData
      if (this.viewed) {
        await this.destroyView('setData');
        this.viewed = false;
      }
    }).then(async () => {
      await this.detachPromise;
      // Wait whether this.dataFrame assigning has called detach() before continue set data and build view

      // -- PDB data --
      let pdbTag: string = pdbTAGS.PDB;
      if (this.pdbTag) pdbTag = this.pdbTag;
      this.pdbStr = this.dataFrame.getTag(pdbTag);
      if (this.pdb && this.pdb != pdbDefault) this.pdbStr = this.pdb;

      // -- Ligand --
      // TODO: Ligand
    }).then(async () => {
      if (!this.viewed) {
        await this.buildView('setData');
        this.viewed = true;
      }
    }).catch((reason: any) => {
      grok.shell.error(reason.toString());
    }).finally(() => {
      this.setDataInProgress = false;
    });
  }

  // -- View --

  private viewPromise: Promise<void> = Promise.resolve();
  private detachPromise: Promise<void> = Promise.resolve();
  private setDataInProgress: boolean = false;

  private viewerDiv?: HTMLDivElement;
  private viewer?: RcsbViewer;

  /** Container to store prop values while {@link this.viewer} is not created yet */
  private viewerProps: Partial<RcsbViewerProps> = {
    [PROPS.layoutShowControls]: false,
    [PROPS.showWelcomeToast]: false,
  };

  private splashDiv?: HTMLDivElement;
  private viewSubs: Unsubscribable[] = [];


  private async destroyView(purpose: string): Promise<void> {
    _package.logger.debug(`MolstarViewer.destroyView( purpose='${purpose}' ) `);
    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }

    if (this.viewer) {
      // Clear viewer
      await this.viewer.clear();
    }
  }

  private async buildView(purpose: string): Promise<void> {
    _package.logger.debug(`MolstarViewer.buildView( purpose='${purpose}' ) `);
    if (this.pdbStr)
      await this.buildViewWithPdb();
    else
      await this.buildViewWithoutPdb();

    this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(
      this.rootOnSizeChanged.bind(this)));
    this.viewSubs.push(this.dataFrame.onSelectionChanged.subscribe(
      this.dataFrameOnSelectionChanged.bind(this)));
  }

  private async buildViewWithPdb() {
    if (!this.pdbStr) throw new Error('MolstarViewer.buildViewWithPdb() pdbStr is empty');

    // Fill in viewer
    if (!this.viewerDiv) {
      this.viewerDiv = ui.div([], {
        classes: 'd4-molstar-viewer',
        style: {width: '100%', height: '100%'},
      });
      this.root.style.overflow = 'hidden'; /* Prevent blinking viewer size changed */
      this.root.appendChild(this.viewerDiv);

      const props: Partial<RcsbViewerProps> = {};

      Object.assign(props, this.viewerProps);
      this.viewer = new RcsbViewer(this.viewerDiv, props);
    }
    if (!this.viewer) throw new Error(`The 'viewer' is not created.`);

    await this.viewer.loadStructureFromData(this.pdbStr, 'pdb', false);
    const plugin: PluginContext = this.viewer.plugin;
    this.viewSubs.push(plugin.commands.subscribe(PluginCommands.Layout.Update, () => {

    }));
  }

  private async buildViewWithoutPdb() {
    if (this.viewerDiv) {
      await this.viewer!.clear();
      delete this.viewer;

      $(this.viewerDiv).empty();
      this.viewerDiv.remove();
      delete this.viewerDiv;
    }

    const fileEl: HTMLInputElement = ui.element('input');
    fileEl.type = 'file';
    fileEl.style.display = 'none';
    fileEl.addEventListener('change', async (event) => {
      const k = 11;
      if (fileEl.files != null && fileEl.files.length == 1) {
        const [pdbStr, pdbHelper]: [string, IPdbHelper] = await Promise.all([
          await fileEl.files[0]!.text(), getPdbHelper()]);
        this.setOptions({pdb: pdbStr});
      }
    });
    const fileLink = ui.link('Open...', '', '', {
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

  private updateView(): void {

  }

  private calcSize(): void {
    if (!this.viewer || !this.viewerDiv) return;

    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    _package.logger.debug(`MolstarViewer.calcSize( ${cw.toString()} x ${ch.toString()} )`);

    this.viewerDiv.style.width = `${cw}px`;
    this.viewerDiv.style.height = `${ch}px`;

    this.viewer.handleResize();
  }

  // -- Handle events --

  private rootOnSizeChanged(value: any): void {
    _package.logger.debug('MolstarViewer.rootOnSizeChanged() ');
    this.calcSize();
  }

  private dataFrameOnCurrentRowChanged(value: any): void {
    _package.logger.debug('MolstarViewer.dataFrameOnCurrentRowChanged() ');
  }

  private dataFrameOnSelectionChanged(value: any): void {
    this.viewPromise = this.viewPromise.then(() => {
      const k = 11;
    });
  }
}

export async function initViewer(viewName: string = 'Mol*'): Promise<RcsbViewer> {
  const view = grok.shell.newView(viewName);
  const viewerContainer = view.root;
  const viewer = new RcsbViewer(viewerContainer, castProps(defaults));
  return viewer;
}

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} pdbId ID of an entity in PDB.
 */
export async function byId(pdbId: string) {
  initViewer()
    .then((v: any) => {
      v.loadPdb(pdbId);
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
    .then(async (viewer: RcsbViewer) => {
      await viewer.loadStructureFromData(pdbData, 'pdb', false);
    });
  //v.handleResize();
}

export async function viewMolstarUI(content: string, name?: string): Promise<void> {
  await byData(content, name);
}

/** Creates view with Molstar viewer to preview Biostructure (PDB) */
export function previewMolstarUI(file: DG.FileInfo): DG.View {
  const view = DG.View.create();
  const viewer = new RcsbViewer(view.root, castProps(defaults));
  const subs: Unsubscribable[] = [];
  subs.push(ui.onSizeChanged(view.root).subscribe((value: any) => {
    viewer.handleResize();
  }));
  subs.push(grok.events.onViewRemoved.subscribe((evtView) => {
    if (evtView.id === view.id)
      for (const sub of subs) sub.unsubscribe();
  }));

  function loadString(data: string) {
    const binary: boolean = false;
    viewer.loadStructureFromData(data, file.extension as BuiltInTrajectoryFormat, binary)
      .then(() => {}); // Ignoring Promise returned
  }

  function loadBytes(bytes: any) {
    const binary: boolean = false;
    viewer.loadStructureFromData(bytes, file.extension as BuiltInTrajectoryFormat, binary)
      .then(() => {}); // Ignoring Promise returned
  }

  if (['pdb', 'sdf'].includes(file.extension))
    file.readAsString().then(loadString);
  else
    file.readAsBytes().then(loadBytes);

  return view;
}

function castProps(src: BiostructureProps): Partial<RcsbViewerProps> {
  const res: Partial<RcsbViewerProps> = {
    showImportControls: src.showImportControls,
    showWelcomeToast: src.showWelcomeToast,
    layoutIsExpanded: src.layoutIsExpanded,
    layoutShowLog: src.layoutShowLog,
    layoutControlsDisplay: src.layoutControlsDisplay,
    layoutShowControls: src.layoutShowControls,
    layoutShowSequence: src.layoutShowSequence,
    backgroundColor: src.backgroundColor as msColor,
    detachedFromSierra: src.detachedFromSierra,
    showExportControls: src.showExportControls,
    showMembraneOrientationPreset: src.showMembraneOrientationPreset,
    showSessionControls: src.showSessionControls,
    showStructureSourceControls: src.showStructureSourceControls,
    viewportShowExpand: src.viewportShowExpand,
    showSuperpositionControls: src.showSuperpositionControls,
    viewportShowSelectionMode: src.viewportShowSelectionMode,
    volumeStreamingServer: src.volumeStreamingServer,
    modelUrlProviders: [],
    extensions: []
  };
  return res;
}
