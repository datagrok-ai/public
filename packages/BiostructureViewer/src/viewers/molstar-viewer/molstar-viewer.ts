import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

// import {Viewer as rcsbViewer, ViewerProps as rcsbViewerProps} from '@rcsb/rcsb-molstar/src/viewer';
import {Viewer as RcsbViewer, ViewerProps as RcsbViewerProps} from '@rcsb/rcsb-molstar/build/src/viewer';

import $ from 'cash-dom';
import wu from 'wu';
import {Observable, Subject, Unsubscribable} from 'rxjs';

import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {PluginLayoutControlsDisplay, PluginLayoutStateProps} from 'molstar/lib/mol-plugin/layout';
import {StructureComponentRef} from 'molstar/lib/mol-plugin-state/manager/structure/hierarchy-state';
import {BuiltInTrajectoryFormat, TrajectoryFormatProvider} from 'molstar/lib/mol-plugin-state/formats/trajectory';
import {PluginCommands} from 'molstar/lib/mol-plugin/commands';
import {to_mmCIF} from 'molstar/lib/mol-model/structure/export/mmcif';
import {utf8ByteCount, utf8Write} from 'molstar/lib/mol-io/common/utf8';
import {zip} from 'molstar/lib/mol-util/zip/zip';
import {SyncRuntimeContext} from 'molstar/lib/mol-task/execution/synchronous';
import {AssetManager} from 'molstar/lib/mol-util/assets';

import {testEvent} from '@datagrok-libraries/utils/src/test';
import {
  BiostructureData, BiostructureDataJson, BiostructureDataProviderFunc
} from '@datagrok-libraries/bio/src/pdb/types';
import {
  IBiostructureViewer,
  PluginLayoutControlsDisplayType,
  RegionStateOptionsType,
  SimpleRegionStateOptionsType,
} from '@datagrok-libraries/bio/src/viewers/molstar-viewer';
import {TAGS as pdbTAGS} from '@datagrok-libraries/bio/src/pdb/index';
import {Molecule3DUnits} from '@datagrok-libraries/bio/src/molecule-3d/molecule-3d-units-handler';
import {CHEM_ATOM_PICKER_LINKED_COL, IMolecule3DBrowser, Molecule3DData} from '@datagrok-libraries/bio/src/viewers/molecule3d';
import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {getDataProviderList} from '@datagrok-libraries/bio/src/utils/data-provider';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';

import {addLigandOnStage, buildSplash, LigandData, parseAndVisualsData, removeVisualsData} from './molstar-viewer-open';
import {
  AtomMapping3D, CHEM_MOL3D_HOVER_EVENT, CHEM_SELECTION_EVENT,
  ChemSelectionEventArgs, Mol3DHoverEventArgs,
  computeSerials, getSelectionCache, selectionCacheKey,
} from './molstar-highlight-utils';
import {defaults, molecule3dFileExtensions} from './consts';
import {createRcsbViewer, disposeRcsbViewer} from './utils';

import {_package} from '../../package';
import {convertWasm} from '../../conversion/wasm/converterWasm';
import {StateObjectRef} from 'molstar/lib/mol-state';
import {createStructureRepresentationParams} from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params';
import {StructureRepresentationRegistry} from 'molstar/lib/mol-repr/structure/registry';
import {StateTransforms} from 'molstar/lib/mol-plugin-state/transforms';
import {StateElements} from 'molstar/lib/examples/proteopedia-wrapper/helpers';
import {MolScriptBuilder} from 'molstar/lib/mol-script/language/builder';
import {Script} from 'molstar/lib/mol-script/script';
import {Bond, Structure, StructureSelection, StructureElement, StructureProperties} from 'molstar/lib/mol-model/structure';
import {InteractivityManager} from 'molstar/lib/mol-plugin-state/manager/interactivity';
import {Color} from 'molstar/lib/mol-util/color';
import {
  setStructureOverpaint, clearStructureOverpaint,
} from 'molstar/lib/mol-plugin-state/helpers/structure-overpaint';

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
  BEHAVIOUR = 'Behaviour',
}

export enum PROPS {
  // -- Data --
  dataJson = 'dataJson',
  pdb = 'pdb',
  pdbTag = 'pdbTag',
  ligandColumnName = 'ligandColumnName',
  // pdbProvider = 'pdbProvider',
  // emdbProvider = 'emdbProvider',
  biostructureIdColumnName = 'biostructureIdColumnName',
  /** DG.Func nqName */
  biostructureDataProvider = 'biostructureDataProvider',
  ligandValue = 'ligandValue',

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
  // -- Behaviour --
  showSelectedRowsLigands = 'showSelectedRowsLigands',
  showCurrentRowLigand = 'showCurrentRowLigand',
  showMouseOverRowLigand = 'showMouseOverRowLigand',
  zoom = 'zoom',
}

const pdbDefault: string = '';

export type LigandMapItem = {
  rowIdx: number,
  structureRefs: Array<string> | null
};
export type LigandMap = {
  selected: LigandMapItem[],
  current: LigandMapItem | null,
  hovered: LigandMapItem | null,
};

export const DebounceIntervals = {
  setData: 20,
  currentRow: 20,
  ligands: 20,
};

export type LigandValueData = {
  value: string,
  semType: string,
  units: string,
}

export class MolstarViewer extends DG.JsViewer implements IBiostructureViewer, IMolecule3DBrowser {
  private viewed: boolean = false;

  // -- Data --
  dataJson: string;
  pdb: string;
  pdbTag: string;

  /** Column name */
  biostructureIdColumnName: string | null;
  /** DG.Func nqName */
  biostructureDataProvider: string | null;

  [PROPS.ligandColumnName]: string;
  // [PROPS.pdbProvider]: string;
  // [PROPS.emdbProvider]: string;
  [PROPS.ligandValue]: string;

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

  // -- Behaviour --
  [PROPS.showSelectedRowsLigands]: boolean;
  [PROPS.showCurrentRowLigand]: boolean;
  [PROPS.showMouseOverRowLigand]: boolean;
  // propsEngine = new class {
  //
  // }(this);
  [PROPS.zoom]: boolean;

  constructor() {
    super();

    this.helpUrl = '/help/visualize/viewers/biostructure';

    // -- Data --
    this.dataJson = this.string(PROPS.dataJson, defaults.dataJson, {
      category: PROPS_CATS.DATA, userEditable: false,
      description: 'JSON encoded object of BiostructureData type with data value Base64 encoded data',
    });
    this.pdb = this.string(PROPS.pdb, pdbDefault,
      {category: PROPS_CATS.DATA, userEditable: false});
    this.pdbTag = this.string(PROPS.pdbTag, defaults.pdbTag,
      {category: PROPS_CATS.DATA, choices: []});

    this.biostructureIdColumnName = this.string(PROPS.biostructureIdColumnName, defaults.biostructureIdColumnName,
      {category: PROPS_CATS.DATA});
    this.biostructureDataProvider = this.string(PROPS.biostructureDataProvider, defaults.biostructureDataProvider,
      {category: PROPS_CATS.DATA, /* fill choices in setData() */});

    this.ligandColumnName = this.string(PROPS.ligandColumnName, defaults.ligandColumnName,
      {category: PROPS_CATS.DATA});
    // this.pdbProvider = this.string(PROPS.pdbProvider, defaults.pdbProvider,
    //   {category: PROPS_CATS.DATA});
    // this.emdbProvider = this.string(PROPS.emdbProvider, defaults.emdbProvider,
    //   {category: PROPS_CATS.DATA});
    this.ligandValue = this.string(PROPS.ligandValue, null, {category: PROPS_CATS.DATA, userEditable: false});

    // -- Layout --
    this.layoutIsExpanded = this.bool(PROPS.layoutIsExpanded, defaults.layoutIsExpanded,
      {category: PROPS_CATS.LAYOUT});
    this.layoutShowControls = this.bool(PROPS.layoutShowControls, defaults.layoutShowControls,
      {category: PROPS_CATS.LAYOUT});
    this.layoutRegionStateLeft = this.string(PROPS.layoutRegionStateLeft, defaults.layoutRegionStateLeft,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(RegionStateOptionsType)}) as RegionStateOptionsType;
    this.layoutRegionStateTop = this.string(PROPS.layoutRegionStateTop, defaults.layoutRegionStateTop,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SimpleRegionStateOptionsType)},
    ) as SimpleRegionStateOptionsType;
    this.layoutRegionStateRight = this.string(PROPS.layoutRegionStateRight, defaults.layoutRegionStateRight,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SimpleRegionStateOptionsType)},
    ) as SimpleRegionStateOptionsType;
    this.layoutRegionStateBottom = this.string(PROPS.layoutRegionStateBottom, defaults.layoutRegionStateBottom,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SimpleRegionStateOptionsType)},
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

    // -- Behaviour --
    this.showSelectedRowsLigands = this.bool(PROPS.showSelectedRowsLigands, false,
      {category: PROPS_CATS.BEHAVIOUR});
    this.showCurrentRowLigand = this.bool(PROPS.showCurrentRowLigand, true,
      {category: PROPS_CATS.BEHAVIOUR});
    this.showMouseOverRowLigand = this.bool(PROPS.showMouseOverRowLigand, true,
      {category: PROPS_CATS.BEHAVIOUR});

    // -- Style --
    this.representation = this.string(PROPS.representation, defaults.representation,
      {category: PROPS_CATS.STYLE, choices: Object.keys(StructureRepresentationRegistry.BuiltIn)});

    // --
    this.zoom = this.bool(PROPS.zoom, false, {userEditable: false});
    this.root.style.textAlign = 'center';

    this.logger = _package.logger;
    this.viewSyncer = new PromiseSyncer(this.logger);

    this.setDataRequest = new Subject<void>();
    this.subs.push(DG.debounce(this.setDataRequest, DebounceIntervals.setData)
      .subscribe(() => { this.onSetDataRequestDebounced(); }));
    this.viewSubs.push(this.onContextMenu.subscribe(this.onContextMenuHandler.bind(this)));
    this._initButtonExpand();
    this._initMenu();
  }

  private static viewerCounter: number = -1;
  private readonly viewerId: number = ++MolstarViewer.viewerCounter;

  private _initButtonExpand() {
    const button = $('.msp-btn.msp-btn-icon.msp-btn-link-toggle-off');
    button.on('click', () => this.root.requestFullscreen());
  }

  private async _initProps() {
    if (!this.dataFrame || this.dataJson) return;

    // -- Pdb or Pdb Id --
    if (!this.biostructureIdColumnName) {
      const pdbCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE3D) || this.dataFrame.columns.bySemType(DG.SEMTYPE.PDB_ID);
      if (pdbCol) {
        this.biostructureIdColumnName = pdbCol.name;

        if (pdbCol.semType === DG.SEMTYPE.PDB_ID) {
          const funcs = await getDataProviderList(DG.SEMTYPE.MOLECULE3D);
          this.biostructureDataProvider = funcs[0]?.nqName;
        }
      }
    }

    // -- Ligand --
    if (!this.ligandColumnName) {
      const molCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE3D) ?? this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
      if (molCol)
        this.ligandColumnName = molCol.name;
    }

    // Set the two-way link tag on the SMILES column so the Chem cell
    // renderer knows which Molecule3D column activates the atom picker.
    // Persist via `col.tags[...]` so the link survives save/reload; use
    // the shared constant from `@datagrok-libraries/bio` to stay in sync
    // with the BSV widget (embedded in Docking's AutoDock panel), the
    // Docking pipeline post-run write, and Chem's reader.
    if (this.ligandColumnName) {
      const ligandCol = this.dataFrame.col(this.ligandColumnName);
      if (ligandCol) {
        const smilesCol = this.dataFrame.columns.toList().find(
          (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE && c.name !== this.ligandColumnName);
        if (smilesCol)
          smilesCol.tags[CHEM_ATOM_PICKER_LINKED_COL] = this.ligandColumnName;
      }
    }
  }

  private viewerToLog(): string { return `MolstarViewer<${this.viewerId}>`; }

  override onPropertyChanged(property: DG.Property | null): void {
    const logIndent: number = 0;
    const callLog: string = `onPropertyChanged( '${property?.name}' )`;
    const logPrefix: string = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);
    super.onPropertyChanged(property);

    if (!property) {
      this.logger.debug(`${logPrefix}, property is null`);
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
    const applyProperty = async (_propName: string, _value: any) => {
      this.logger.debug(`${logPrefix}.applyProperty(), start`);

      if (!this.viewer) throw new Error('viewer does not exists');

      const plugin: PluginContext = this.viewer.plugin;
      const state: Partial<PluginLayoutStateProps> = {};
      switch (property.name) {
      case PROPS.layoutShowLog: {
        // PluginCommands.Layout.Update(plugin, );
        const _k = 11;
        break;
      }
      case PROPS.layoutShowControls: {
        state.showControls = this.layoutShowControls;
        break;
      }
      case PROPS.layoutIsExpanded: {
        state.isExpanded = this.layoutIsExpanded;
        break;
      }
      case PROPS.layoutRegionStateLeft:
      case PROPS.layoutRegionStateTop:
      case PROPS.layoutRegionStateRight:
      case PROPS.layoutRegionStateBottom: {
        state.regionState = {
          left: this.layoutRegionStateLeft,
          top: this.layoutRegionStateTop,
          right: this.layoutRegionStateRight,
          bottom: this.layoutRegionStateBottom,
        };
        break;
      }
      case PROPS.layoutControlsDisplay: {
        state.controlsDisplay = this.layoutControlsDisplay as PluginLayoutControlsDisplay;
        break;
      }
        // case PROPS.viewportShowExpand: {
        //   await PluginCommands.State.ToggleExpanded(plugin,
        //     {state: {isExpanded: this.viewportShowExpand}})
        //   break;
        // }
      }
      //eslint-disable-next-line new-cap
      await PluginCommands.Layout.Update(plugin, {state: state});
      this.logger.debug(`${logPrefix}.applyProperty(), end`);
    };

    // Type checks
    switch (property.name) {
    case PROPS.dataJson: {
      if (typeof this.dataJson !== 'string')
        throw new Error(`Property '${property.name}' value is of invalid type '${typeof this.dataJson}'.`);
      break;
    }
    }

    switch (property.name) {
    case PROPS.biostructureDataProvider: {
      // reset for setData
      this.biostructureDataProviderFunc = null;
      break;
    }
    case PROPS.pdb:
    case PROPS.dataJson: {
      // reset for setData
      this.dataEff = null;
      break;
    }

    case PROPS.representation: {
      const checkAndUpdate = () => {
        if (!this.setDataInProgress && this.viewer)
          this.updateView(this.representation);
        else
          requestAnimationFrame(checkAndUpdate);
      };
      checkAndUpdate();
      break;
    }

    case PROPS.showImportControls:
      break;
    case PROPS.layoutIsExpanded:
      this.viewerProps[PROPS.layoutIsExpanded] = this.layoutIsExpanded;
      break;
    case PROPS.showSelectedRowsLigands:
    case PROPS.showCurrentRowLigand:
    case PROPS.showMouseOverRowLigand:
      this.onRebuildViewLigandsRequest.next();
      break;
    }

    const propName: string = property.name;
    const propValue: any = this.props.get(propName);
    this.viewSyncer.sync(logPrefix, async () => {
      if (this.viewer)
        await applyProperty(propName, propValue);
    });

    switch (property.name) {
    case PROPS.dataJson:
    case PROPS.pdb:
    case PROPS.pdbTag:
    case PROPS.biostructureDataProvider:
    case PROPS.biostructureIdColumnName:
    case PROPS.ligandValue:
    case PROPS.ligandColumnName: {
      this.setData(logIndent + 1, callLog);
      break;
    }
    }
    this.logger.debug(`${logPrefix}, end`);
  }

  // effective PDB value (to plot)
  /** {@link IBiostructureViewer.dataEff} */
  public dataEff: BiostructureData | null = null;
  private dataEffStructureRefs: string[] | null = null;

  private biostructureDataProviderList: DG.Func[] = [];
  private biostructureDataProviderFunc: DG.Func | null = null;

  override onTableAttached(): void {
    const logIndent: number = 0;
    const callLog = 'onTableAttached()';
    const logPrefix: string = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    const superOnTableAttached = super.onTableAttached.bind(this);

    // -- Props editors --
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.pdbTag).choices = ['', ...dfTagNameList];

    // Update biostructureDataProviderList in queue before setData, sync handles error also
    // TODO: Optimize
    this.viewSyncer.sync(logPrefix, async () => {
      this.biostructureDataProviderList = await getDataProviderList(DG.SEMTYPE.MOLECULE3D);
      this.props.getProperty(PROPS.biostructureDataProvider).choices =
        ['', ...this.biostructureDataProviderList.map((f) => f.nqName)];
    });
    this._initProps();

    //this.subs.push(this.onContextMenu.subscribe(this.onContextMenuHandler.bind(this)));

    superOnTableAttached();
    this.setData(logIndent + 1, callLog);
    ui.tools.waitForElementInDom(this.root).then(() => this.initializeResizeHandling());
    this.logger.debug(`${logPrefix}, end`);
  }

  initializeResizeHandling(): void {
    const dialogPanel = this.root.closest('.dialog-floating') as HTMLElement;
    const accPanel = this.root.closest('.d4-accordion-pane-content') as HTMLElement;
    const resizeTarget = dialogPanel || accPanel;

    if (!resizeTarget) return;

    const setDimensions = (element: HTMLElement | null, width: number, height: number) => {
      if (element) {
        element.style.width = `${width}px`;
        element.style.height = `${height}px`;
      }
    };

    const resizeCanvasAndViewer = (width: number, height: number) => {
      const canvas = this.viewer?.plugin.canvas3dContext?.canvas;
      setDimensions(this.viewerDiv!, width, height);
      if (canvas)
        setDimensions(canvas, width, height);
      this.viewer?.plugin.canvas3d?.handleResize();
      this.viewer?.plugin.handleResize();
    };

    this.subs.push(ui.onSizeChanged(resizeTarget).subscribe(() => {
      const width = resizeTarget.clientWidth;
      const height = resizeTarget.clientHeight;

      if (dialogPanel) {
        setDimensions(this.root, width, height);
        resizeCanvasAndViewer(width, height);

        const waitParentEl = this.root.parentElement;
        if (waitParentEl?.classList.contains('grok-wait'))
          setDimensions(waitParentEl, width, height);
      }

      if (accPanel)
        resizeCanvasAndViewer(width, height);
    }));
  }

  override detach(): void {
    const logIndent: number = 0;
    const callLog: string = `detach()`;
    const logPrefix: string = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    const superDetach = super.detach.bind(this);
    this.viewSyncer.sync(logPrefix, async () => { // detach
      if (this.setDataInProgress) return; // check setDataInProgress synced
      if (this.viewed) {
        await this.destroyView(true, logIndent + 1, callLog);
        this.viewed = false;
      }
      superDetach();
    });
    this.logger.debug(`${logPrefix}, end`);
  }

  _initMenu(): void {
    this.root.oncontextmenu = (e) => {
      e.preventDefault();
      e.stopImmediatePropagation();
      const menu = DG.Menu.popup();
      this.onContextMenuHandler(menu);
      menu.show();
    };
  }

  onContextMenuHandler(menu: DG.Menu): void {
    menu
      .group('Download')
      .item('As CIF', async () => await this.writeToFormat('cif'))
      .item('As PDB', async () => await this.writeToFormat('pdb'));
  }

  async writeToFormat(format: 'cif' | 'pdb') {
    const fileName = 'file';
    const {structures} = this.viewer!.plugin.managers.structure.hierarchy.current;
    const _s = structures[0];
    const s = _s.transform?.cell.obj?.data ?? _s.cell.obj?.data;
    const cif = to_mmCIF(fileName, s!) as string;

    if (format === 'cif')
      DG.Utils.download(`${fileName}.cif`, cif);
    else if (format === 'pdb') {
      const pdb = await convertWasm(cif);
      DG.Utils.download(`${fileName}.pdb`, pdb);
    }
  }

  async _exportHierarchy(plugin: PluginContext, options?: { format?: 'cif' }) {
    const format = options?.format ?? 'cif';
    const {structures} = plugin.managers.structure.hierarchy.current;
    const files: [name: string, data: string | Uint8Array][] = [];
    const entryMap = new Map<string, number>();
    for (const _s of structures) {
      const s = _s.transform?.cell.obj?.data ?? _s.cell.obj?.data;
      if (!s) continue;
      const name = s.model.entryId || 'unnamed';

      const fileName = entryMap.has(name) ?
        `${name}_${entryMap.get(name)! + 1}.${format}` :
        `${name}.${format}`;
      entryMap.set(name, (entryMap.get(name) ?? 0) + 1);

      if (s.elementCount > 100000)
        await new Promise((res) => setTimeout(res, 50));

      try {
        files.push([fileName, to_mmCIF(name, s)]);
      } catch (e) {
        this.logger.error(e);
      }
    }

    let blob;
    if (files.length === 1)
      blob = new Blob([files[0][1] as BlobPart]), files[0][0];
    else if (files.length > 1) {
      const zipData: Record<string, Uint8Array> = {};
      for (const [fn, data] of files) {
        if (data instanceof Uint8Array)
          zipData[fn] = data;
        else {
          const bytes = new Uint8Array(utf8ByteCount(data));
          utf8Write(bytes, 0, data);
          zipData[fn] = bytes;
        }
      }
      const assetManager = new AssetManager();
      const ctx = {runtime: SyncRuntimeContext, assetManager};
      const buffer = await zip(ctx.runtime, zipData);
      //@ts-ignore
      blob = new Blob([new Uint8Array(buffer, 0, buffer.byteLength)]);
    }
    return await blob?.text();
  }

  // -- Data --

  private _setDataCallCounter = -1;

  private setDataRequest: Subject<void>;

  setData(logIndent: number, caller: string): void {
    const callLog = `setData( <- ${caller})`;
    this.logger.debug(`${margin(logIndent)}${this.viewerToLog()}.${callLog}`);
    this.setDataRequest.next();
  }

  onSetDataRequestDebounced(): void {
    const logIndent: number = 0;
    const callId = ++this._setDataCallCounter;
    const callLog = `onSetDataRequestDebounced( callId = ${callId} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    this.viewSyncer.sync(logPrefix, async () => { // setData
      if (!this.setDataInProgress) this.setDataInProgress = true; else return; // check setDataInProgress synced
      try {
        if (this.viewed) {
          await this.destroyView(false, 1, callLog);
          this.viewed = false;
        }

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

        if (this.biostructureDataProvider) {
          if (!this.biostructureDataProviderFunc) {
            this.biostructureDataProviderFunc = this.biostructureDataProviderList
              .find((f) => f.nqName === this.biostructureDataProvider) ?? null;
          }
        }

        if (!this.viewed) {
          await this.buildView(1, callLog);
          this.viewed = true;
        }
      } finally {
        this.setDataInProgress = false;
      }
    });
    this.logger.debug(`${logPrefix}, end`);
  }

  // -- View --

  private readonly logger: ILogger;
  private readonly viewSyncer: PromiseSyncer;
  private setDataInProgress: boolean = false;

  private viewerDiv?: HTMLDivElement;
  private viewer?: RcsbViewer;

  /** Container to store prop values while {@link viewer} is not created yet */
  private viewerProps: Partial<RcsbViewerProps> = {
    [PROPS.layoutShowControls]: false,
    [PROPS.showWelcomeToast]: false,
  };

  private splashDiv?: HTMLDivElement;
  private viewSubs: Unsubscribable[] = [];
  private onRebuildViewLigandsRequest: Subject<void> = new Subject<void>();

  private async destroyView(free: boolean, logIndent: number, caller: string): Promise<void> {
    const callLog: string = `destroyView( <- ${caller} )`;
    const logPrefix: string = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start `);

    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    if (this.splashDiv) {
      $(this.splashDiv).empty();
      this.splashDiv.remove();
      delete this.splashDiv;
    }

    if (this.dataFrame && this.ligandColumnName)
      await this.destroyViewLigands(logIndent + 1, callLog);

    if (free /* detach */) {
      if (this.viewerDiv) {
        // Clear viewer
        // await this.viewer.clear();
        await disposeRcsbViewer(this.viewer!, this.viewerDiv);
        delete this.viewer;
        this.viewerDiv.remove();
        delete this.viewerDiv;
      }
    } else {
      if (this.dataEffStructureRefs && this.viewer?.plugin)
        await removeVisualsData(this.viewer!.plugin, this.dataEffStructureRefs, callLog);
    }
    this.logger.debug(`${logPrefix}, end `);
  }

  private async buildView(logIndent: number, caller: string): Promise<void> {
    const callLog: string = `buildView( <- ${caller} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}}`;
    this.logger.debug(`${logPrefix}, start `);

    if (!this.dataEff && (!this.dataFrame || !this.biostructureDataProvider && !this.biostructureIdColumnName))
      await this.buildViewWithoutData(logIndent + 1, callLog);
    else
      await this.buildViewWithData(logIndent + 1, callLog);

    this.logger.debug(`${logPrefix}, end `);
  }

  /** Gets structure for the current row from the data provider and displays it. */
  private async rebuildViewCurrentRow(
    oldStructureRefs: string[] | null, newCurrentRowIdx: number, logIndent: number, caller: string
  ): Promise<[BiostructureData, string[]]> {
    const callLog = `rebuildViewCurrentRow( <- ${caller}, currentRowIdx = ${newCurrentRowIdx} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;

    if (!this.viewer)
      this.logger.warning(`${logPrefix}, no viewer`);

    const col = this.dataFrame.getCol(this.biostructureIdColumnName!);
    const id = col.get(newCurrentRowIdx);
    const splash = buildSplash(this.root, `Loading data of '${id}'.`);
    try {
      // while loading the next, the previous structure is covered by splash
      const dataString = col.semType === DG.SEMTYPE.MOLECULE3D ?
        BiostructureDataJson.fromData({binary: false, data: id, ext: 'pdb'}) :
        await (async () => {
          const fetchData = await this.biostructureDataProviderFunc!.prepare({id}).call();
          return fetchData.getOutputParamValue() as string;
        })();

      const dataEff = BiostructureDataJson.toData(dataString);
      const plugin = this.viewer!.plugin;
      await this.destroyViewLigands(0, callLog);
      await removeVisualsData(plugin, oldStructureRefs, callLog);
      const newStructureRefs: string[] = await parseAndVisualsData(plugin, dataEff, callLog);
      await this.buildViewLigands(0, callLog);
      return [dataEff, newStructureRefs];
    } finally {
      splash.close();
    }
  }

  private async buildViewWithData(logIndent: number, caller: string): Promise<void> {
    const callLog: string = `buildViewWithData( <- ${caller} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    // Fill in viewer
    if (!this.viewerDiv) {
      this.viewerDiv = ui.div([], {
        classes: 'd4-molstar-viewer',
        style: {width: '100%', height: '100%'},
      });
      this.root.style.overflow = 'hidden'; /* Prevent blinking viewer size changed */
      this.root.appendChild(this.viewerDiv);
      this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(
        this.rootOnSizeChanged.bind(this)));

      const props: Partial<RcsbViewerProps> = {};

      Object.assign(props, this.viewerProps);
      // this.viewer = new RcsbViewer(this.viewerDiv, props);
      this.viewer = await createRcsbViewer(this.viewerDiv, props);
    }
    if (!this.viewer) throw new Error(`The 'viewer' is not created.`);

    const plugin = this.viewer.plugin;
    this.viewSubs.push(plugin.commands.subscribe(PluginCommands.Layout.Update, () => {

    }));

    const isPdbColumn = this.biostructureIdColumnName ?
      this.dataFrame.getCol(this.biostructureIdColumnName).semType === DG.SEMTYPE.MOLECULE3D :
      false;

    const hasValidDataProvider = this.biostructureDataProviderFunc && !isPdbColumn;
    if (this.dataFrame && (hasValidDataProvider || isPdbColumn)) {
      this.viewSubs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, DebounceIntervals.currentRow).subscribe(
        this.dataFrameOnCurrentRowChangedDebounced.bind(this)));
      [this.dataEff, this.dataEffStructureRefs] = await this.rebuildViewCurrentRow(
        null, this.dataFrame.currentRowIdx, logIndent + 1, callLog);
    } else if (this.dataEff) {
      // display a structure of this.dataEff
      this.dataEffStructureRefs = await parseAndVisualsData(plugin, this.dataEff, callLog);
    }

    if (this.dataEff) {
      if (!(this.dataEff.ext in molecule3dFileExtensions))
        throw new Error(`Unsupported file extension '${this.dataEff.ext}'.`);

      // Subscribe at instance level so we can apply highlights in real-time
      // when the viewer IS alive. Works for BOTH the AutoDock Molstar panel
      // and the 3D Structure Molstar panel. The module-level listener
      // (registered at module load) handles caching for replay.
      this.viewSubs.push(grok.events.onCustomEvent(
        'chem-interactive-selection-changed').subscribe((_args: unknown) => {
        const args = _args as ChemSelectionEventArgs;
        try {
          if (args?.atoms?.length >= 0) {
            this.logger.debug(`[molstar-picker] live highlight atomsLen=${args.atoms.length} rowIdx=${args?.rowIdx} persistent=${args?.persistent}`);
            // Pass the event's atoms + mapping directly so transient
            // (preview) highlights work even though they're not cached.
            this.highlightAllLigandAtoms({
              rowIdx: args.rowIdx,
              atoms: args.atoms,
              mapping3D: args.mapping3D ?? null,
            });
          }
        } catch (err: any) {
          this.logger.error(
            `${CHEM_SELECTION_EVENT} handler failed: ${err?.message ?? err}`);
        }
      }));

      // Reapply highlights on any click/interaction in the Molstar canvas.
      // This catches zoom, focus, right-click etc. that may trigger
      // structure rebuilds which clear overpaint.
      const canvas = plugin.canvas3d?.webgl?.gl?.canvas as HTMLElement | undefined;
      if (canvas) {
        const reapplyOnClick = () => {
          // Queue through syncer so it runs AFTER any triggered rebuilds
          // (which also go through the syncer). This ensures fresh
          // components are used, not stale references from before rebuild.
          if (this._stateChangeTimer) clearTimeout(this._stateChangeTimer);
          this._stateChangeTimer = setTimeout(() => {
            this._stateChangeTimer = null;
            if (getSelectionCache())
              this._replayHighlightIfCached();
          }, 500) as any;
        };
        canvas.addEventListener('click', reapplyOnClick);
        canvas.addEventListener('contextmenu', reapplyOnClick);
        canvas.addEventListener('mouseup', reapplyOnClick);
        this.viewSubs.push({unsubscribe: () => {
          canvas.removeEventListener('click', reapplyOnClick);
          canvas.removeEventListener('contextmenu', reapplyOnClick);
          canvas.removeEventListener('mouseup', reapplyOnClick);
        }} as any);
      }

      // --- Reverse 3D→2D hover bridge ------------------------------------
      // Subscribe to Molstar's hover observable. When the user hovers a
      // LIGAND (non-polymer) atom in the 3D viewer, fire
      // CHEM_MOL3D_HOVER_EVENT so Chem's rdkit-cell-renderer can highlight
      // the corresponding 2D atom. Receptor atoms are filtered out; H atoms
      // are ignored (2D renderings have only heavy atoms). The `mode` field
      // mirrors the 2D modifier scheme (shift = paint, ctrl+shift = erase).
      //
      // Dedup: we track the last (rowIdx, atomSerial, mode) tuple and skip
      // repeat fires so the downstream 2D renderer isn't flooded with
      // redundant work on every mousemove pixel within the same atom radius.
      this._lastHoverFired = null;
      try {
        this.viewSubs.push(plugin.behaviors.interaction.hover.subscribe(
          (event: InteractivityManager.HoverEvent) => {
            try {
              if (!this.dataFrame || !this.ligandColumnName) return;
              const rowIdx = this.dataFrame.currentRowIdx;
              if (rowIdx == null || rowIdx < 0) return;

              const loci = event.current.loci;
              if (loci.kind === 'empty-loci') {
                this._fireMol3DHover(rowIdx, null, 'preview');
                return;
              }

              // Resolve the hovered atom from whichever loci kind Molstar
              // produced. Element-loci (direct atom pick) is the ideal case.
              // Bond-loci is what Molstar's default ball-and-stick picker
              // emits when hovering between atoms — pick the `a` side so the
              // user still gets a single-atom highlight in 2D. Other kinds
              // (structure-loci, volume-loci, etc.) are ignored.
              let loc: StructureElement.Location | null = null;
              if (loci.kind === 'element-loci') {
                if (StructureElement.Loci.isEmpty(loci)) {
                  this._fireMol3DHover(rowIdx, null, 'preview');
                  return;
                }
                loc = StructureElement.Loci.getFirstLocation(loci) ?? null;
              } else if (loci.kind === 'bond-loci' && loci.bonds.length > 0) {
                const bond = loci.bonds[0] as Bond.Location;
                loc = StructureElement.Location.create(loci.structure, bond.aUnit, bond.aUnit.elements[bond.aIndex]);
              }
              if (!loc) return;
              // Filter: only ligand (non-polymer) atoms. Receptor residues ignored.
              // When the ligand is loaded as a standalone SDF (no receptor), every
              // atom is 'non-polymer' anyway, so this is a no-op in that case.
              try {
                if (StructureProperties.entity.type(loc) !== 'non-polymer') return;
              } catch { /* entity lookup failed — skip rather than risk false-pos */ return; }
              // Skip explicit H — 2D SMILES have only heavy atoms.
              if (StructureProperties.atom.type_symbol(loc) === 'H') return;

              const atomSerial = StructureProperties.atom.id(loc);
              const {modifiers} = event;
              const mode: 'preview' | 'paint' | 'erase' =
                modifiers.shift && modifiers.control ? 'erase' : (modifiers.shift ? 'paint' : 'preview');
              this._fireMol3DHover(rowIdx, atomSerial, mode);
            } catch (err: unknown) {
              this.logger.error(
                `${CHEM_MOL3D_HOVER_EVENT} emit failed: ${err instanceof Error ? err.message : String(err)}`);
            }
          },
        ));
      } catch (err: unknown) {
        // Older Molstar builds without `behaviors.interaction.hover` — reverse
        // bridge just won't activate; forward direction remains unaffected.
        this.logger.debug(
          `Molstar hover subscription unavailable: ${err instanceof Error ? err.message : String(err)}`);
      }

      if (this.dataFrame && this.ligandColumnName) {
        this.viewSubs.push(this.dataFrame.onSelectionChanged.subscribe(
          this.dataFrameOnSelectionChanged.bind(this)));
        this.viewSubs.push(this.dataFrame.onCurrentRowChanged.subscribe(
          this.dataFrameOnCurrentRowChanged.bind(this)));
        this.viewSubs.push(this.dataFrame.onMouseOverRowChanged.subscribe(
          this.dataFrameOnMouseOverRowChanged.bind(this)));
        this.viewSubs.push(DG.debounce(this.onRebuildViewLigandsRequest, DebounceIntervals.ligands).subscribe(
          this.onRebuildViewLigandsDebounced.bind(this)));

        await this.buildViewLigands(logIndent, callLog);
      } else {
        // Non-ligand viewer (e.g. "3D Structure" panel showing a single
        // molecule). Replay any cached selection after the structure loads.
        this._replayHighlightIfCached();
      }
    }

    this.logger.debug(`${logPrefix}, end`);
  }

  private async buildViewWithoutData(logIndent: number, caller: string) {
    const callLog: string = `buildViewWithoutData( <- ${caller} )`;
    const logPrefix: string = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    if (this.viewerDiv) {
      await disposeRcsbViewer(this.viewer!, this.viewerDiv);
      delete this.viewer;
      this.viewerDiv.remove();
      delete this.viewerDiv;
    }

    const dataFileProp = DG.Property.fromOptions({name: 'dataFile', caption: 'Data file', type: 'file'});
    const dataFileInput = DG.InputBase.forProperty(dataFileProp);
    dataFileInput.captionLabel.innerText = 'Data file';
    this.viewSubs.push(dataFileInput.onChanged.subscribe(async () => {
      const dataFi: DG.FileInfo = dataFileInput.value;
      const dataA: Uint8Array = dataFi.data ? dataFi.data /* User's file*/ : await dataFi.readAsBytes()/* Shares */;
      const data: BiostructureData = {binary: true, ext: dataFi.extension, data: dataA};
      this.setOptions({
        [PROPS.dataJson]: BiostructureDataJson.fromData(data),
      });
      // openBtn.disabled = !(dataFi && (dataFi.data || dataFi.url));
    }));
    this.splashDiv = ui.div(
      ui.divV([dataFileInput]/*, openBtn*/),
      {classes: 'bsv-viewer-splash'} /* splash */);
    this.root.appendChild(this.splashDiv);
    this.logger.debug(`${logPrefix}, end`);
  }

  private async updateView(type: any) {
    const entries = this.viewer!.plugin.managers.structure.selection.entries;
    const state = this.viewer!.plugin.state;
    entries.forEach(async ({selection}, ref) => {
      const cell = StateObjectRef.resolveAndCheck(state!.data, ref);
      if (cell) {
        const components = this.viewer!.plugin.build().to(cell);
        const repr = createStructureRepresentationParams(this.viewer!.plugin, void 0, {
          type: type
        });
        components.applyOrUpdate(StateElements.SequenceVisual, StateTransforms.Representation.StructureRepresentation3D, repr);
        await components.commit();
      }
    });
  }

  private calcSize(logIndent: number, caller: string): void {
    const callLog = `calcSize( <- ${caller} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;

    if (this.viewerDiv && this.viewer) {
      this.viewerDiv.style.width = `${cw}px`;
      this.viewerDiv.style.height = `${ch}px`;
      this.viewer.handleResize();
    }

    if (this.splashDiv) {
      this.splashDiv.style.width = `${cw}px`;
      this.splashDiv.style.height = `${ch}px`;
      this.splashDiv.style.lineHeight = `${ch}px`;
    }
    this.logger.debug(`${logPrefix}, end`);
  }

  // -- Handle events --

  private rootOnSizeChanged(_value: any): void {
    const logIndent: number = 0;
    const callLog = `rootOnSizeChanged()`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);
    this.calcSize(logIndent + 1, callLog);
    this.logger.debug(`${logPrefix}, end`);
  }

  private dataFrameOnCurrentRowChangedDebounced(_value: any): void {
    const logIndent: number = 0;
    const oldStructureRefs: string[] | null = this.dataEffStructureRefs;
    const newCurrentRowIdx = this.dataFrame.currentRowIdx;
    const callLog = `dataFrameOnCurrentRowChangedDebounced( newCurrentRowIdx = ${newCurrentRowIdx} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.viewSyncer.sync(logPrefix, async () => {
      [this.dataEff, this.dataEffStructureRefs] = await this.rebuildViewCurrentRow(
        oldStructureRefs, newCurrentRowIdx, 0, callLog);
      // Protein structure was rebuilt — reapply highlights on ligands.
      this._replayHighlightIfCached();
    });
  }

  private dataFrameOnSelectionChanged(_value: any): void {
    this.logger.debug(`${this.viewerToLog()}.dataFrameOnSelectionChanged() `);
    this.onRebuildViewLigandsRequest.next();
  }

  private dataFrameOnCurrentRowChanged(_value: any): void {
    this.logger.debug(`${this.viewerToLog()}.dataFrameOnCurrentRowChanged() `);
    this.onRebuildViewLigandsRequest.next();
  }

  private dataFrameOnMouseOverRowChanged(_value: any): void {
    this.logger.debug(`${this.viewerToLog()}.dataFrameOnMouseOverRowChanged() `);
    this.onRebuildViewLigandsRequest.next();
  }

  private onRebuildViewLigandsDebounced(): void {
    const callLog = `onRebuildViewLigandsDebounced()`;
    this.logger.debug(`${this.viewerToLog()}.${callLog} `);
    this.rebuildViewLigands(0, callLog);
  }

  // -- Ligands routines --

  public ligands: LigandMap = {selected: [], current: null, hovered: null};

  /** Unify get mol* component key/ref, not static for performance
   * @param {StructureComponentRef} comp
   * @return {string | null} */
  getCompKey(comp: StructureComponentRef): string | null {
    return comp.cell.sourceRef ?? null; // comp.version
  }

  private getLigandStr(rowIdx: number): LigandData {
    if (!this.dataFrame || !this.ligandColumnName)
      throw new Error(`${this.viewerToLog()}.getLigandStr(), no dataFrame or ligandColumnName`);

    const ligandCol: DG.Column = this.dataFrame.getCol(this.ligandColumnName);
    const ligandObject: LigandValueData = this.ligandValue ? JSON.parse(this.ligandValue) : {};
    const ligandUnits: string = ligandObject.units ?? ligandCol.meta.units;
    const ligandCellValue: string = ligandObject.value ?? ligandCol.get(rowIdx);
    const ligandSemType: string = ligandObject.semType ?? ligandCol.semType;
    let ligandValue: string;
    let ligandFormat: BuiltInTrajectoryFormat | TrajectoryFormatProvider;
    switch (ligandSemType) {
    case DG.SEMTYPE.MOLECULE: {
      switch (ligandUnits) {
      default: {
        ligandFormat = 'sdf';
        ligandValue = ligandCellValue + '$$$$';
      }
      }
      break;
    }

    case DG.SEMTYPE.MOLECULE3D: {
      switch (ligandUnits) {
      case Molecule3DUnits.pdb: {
        ligandFormat = 'pdb';
        ligandValue = ligandCellValue;
        break;
      }
      case Molecule3DUnits.pdbqt: {
        ligandFormat = 'pdbqt';
        ligandValue = ligandCellValue;
        break;
      }
      default:
        throw new Error(`Unsupported units '${ligandUnits}' of '${DG.SEMTYPE.MOLECULE3D}' ligand.`);
      }
      break;
    }

    default:
      throw new Error(`Unsupported ligand semantic type '${ligandSemType}'.`);
    }
    // const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
    // return ligandBlob;
    return {data: ligandValue, format: ligandFormat, rowIdx: rowIdx};
  }

  private rebuildViewLigands(logIndent: number, caller: string): void {
    const callLog = `rebuildViewLigands( <- ${caller} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.viewSyncer.sync(`${logPrefix}`, async () => {
      await this.destroyViewLigands(0, callLog);
      await this.buildViewLigands(0, callLog);
    });
  }

  private async destroyViewLigands(logIndent: number, caller: string): Promise<void> {
    const callLog = `destroyViewLigands( <- ${caller} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start`);

    if (!this.viewer) throw new Error('The viewer is not created'); // return; // There is not PDB data
    if (!this.ligandColumnName) return;

    const allLigands: LigandMapItem[] = [
      ...this.ligands.selected,
      ...(this.ligands.current ? [this.ligands.current] : []),
      ...(this.ligands.hovered ? [this.ligands.hovered] : []),
    ];

    for (const ligand of allLigands) {
      if (!ligand.structureRefs) continue;
      await removeVisualsData(this.viewer!.plugin, ligand.structureRefs, callLog);
    }
    for (const ligand of allLigands) ligand.structureRefs = null; // unbind with this.stage.compList

    this.logger.debug(`${logPrefix}, end`);
  }

  /** Builds up ligands on the stage view with postponed sync */
  private async buildViewLigands(logIndent: number, caller: string): Promise<void> {
    const callLog = `buildViewLigands( <- ${caller} )`;
    const logPrefix = `${margin(logIndent)}${this.viewerToLog()}.${callLog}`;
    this.logger.debug(`${logPrefix}, start `);

    if (!this.viewer) throw new Error('The mol* viewer is not created'); // return; // There is not PDB data
    if (!this.dataFrame || !this.ligandColumnName || (this.dataFrame.col(this.ligandColumnName)?.semType !== DG.SEMTYPE.MOLECULE && this.dataFrame.col(this.ligandColumnName)?.semType !== DG.SEMTYPE.MOLECULE3D)) return;
    const newLigands: LigandMap = {selected: [], current: null, hovered: null};
    newLigands.selected = !this.showSelectedRowsLigands ? [] :
      wu(this.dataFrame.selection.getSelectedIndexes())
        .take(25)
        .map((selRowIdx) => { return {rowIdx: selRowIdx, structureRefs: null}; })
        .toArray();
    newLigands.current = !this.showCurrentRowLigand ? null :
      this.dataFrame.currentRowIdx >= 0 ? {rowIdx: this.dataFrame.currentRowIdx, structureRefs: null} : null;
    newLigands.hovered = !this.showMouseOverRowLigand ? null :
      this.dataFrame.mouseOverRowIdx >= 0 ? {rowIdx: this.dataFrame.mouseOverRowIdx, structureRefs: null} : null;

    const plugin = this.viewer!.plugin;
    const ligandTaskList: (() => Promise<void>)[] = [];

    const selCount = newLigands.selected.length;
    for (const [selectedLigand, selI] of wu.enumerate(newLigands.selected)) {
      const color =
        this.showCurrentRowLigand || this.showMouseOverRowLigand ?
          (selCount > 1 ? DG.Color.selectedRows : null) :
          (selCount > 1 ? DG.Color.scaleColor(selI, 0, selCount, 0.5) : null);

      const selectedLigandData = this.getLigandStr(selectedLigand.rowIdx);
      ligandTaskList.push(async () => {
        selectedLigand.structureRefs = await addLigandOnStage(plugin, selectedLigandData, color, this.zoom);
      });
    }
    // When both current and hovered are loaded, always use distinct
    // colors so the user can tell which pose belongs to which row.
    const hasBothPoses = !!newLigands.current && !!newLigands.hovered;
    if (newLigands.current) {
      const color = hasBothPoses || this.showSelectedRowsLigands ?
        (DG.Color.currentRow || 0xFF8C00) : // dark orange fallback
        null;

      const currentLigandData = this.getLigandStr(newLigands.current.rowIdx);
      const currentLigand = newLigands.current;
      ligandTaskList.push(async () => {
        currentLigand.structureRefs = await addLigandOnStage(plugin, currentLigandData, color, this.zoom);
      });
    }
    if (newLigands.hovered) {
      const color = hasBothPoses || this.showSelectedRowsLigands || this.showCurrentRowLigand ?
        (DG.Color.mouseOverRows || 0x2196F3) : // blue fallback
        null;
      const hoveredLigandData = this.getLigandStr(newLigands.hovered.rowIdx);
      const hoveredLigand = newLigands.hovered;
      ligandTaskList.push(async () => {
        hoveredLigand.structureRefs = await addLigandOnStage(plugin, hoveredLigandData, color, this.zoom);
      });
    }
    this.logger.debug(`${logPrefix},\nnewLigands = ${JSON.stringify(newLigands)}`);

    // Because of the async nature of loading structures to .viewer, the .dataFrame property can be changed (to null).
    // So collect data from the .dataFrame synchronously and then add ligands to the .viewer with postponed sync.
    await Promise.all(ligandTaskList.map(async (task) => {
      try {
        await task();
      } catch (e) {
        _package.logger.error(e);
      }
    })).then(() => {
      this.ligands = newLigands;
    });

    // Replay highlights after ligands are built (via syncer queue).
    this._replayHighlightIfCached();

    this.logger.debug(`${logPrefix}, end`);
  }

  // -- IRenderer--

  private _onRendered: Subject<void> = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    const logPrefix = `${this.viewerToLog()}.invalidate(${caller ? ` <- ${caller} ` : ''})`;
    this.logger.debug(`${logPrefix}, ` + '---> put to syncer');
    // Put the event trigger in the tail of the synced calls queue.
    this.viewSyncer.sync(logPrefix, async () => {
      // update view / render
      this.logger.debug(`${logPrefix}, ` + '---> triggered _onRendered');
      this._onRendered.next();
    });
  }

  async awaitRendered(timeout: number = 10000): Promise<void> {
    const callLog = `awaitRendered( ${timeout} )`;
    const logPrefix = `${this.viewerToLog()}.${callLog}`;
    await testEvent(this._onRendered, (args) => {
      this.logger.debug(`${logPrefix}, ` + '_onRendered event caught');
    }, () => {
      this.invalidate(callLog);
    }, timeout);

    // Rethrow stored syncer error (for test purposes)
    const viewErrors = this.viewSyncer.resetErrors();
    if (viewErrors.length > 0) throw viewErrors[0];
  }

  // -- IMolecule3DBrowser --

  showStructure(data: Molecule3DData) {
    throw new Error('Not implemented');
  }

  /** Guard against concurrent async highlight updates. */
  private _highlightInProgress = false;
  private _stateChangeTimer: ReturnType<typeof setTimeout> | null = null;

  /** Last-fired 3D→2D hover tuple — used to dedup Molstar's
   *  per-mousemove-pixel hover firing. null = no atom currently hovered. */
  private _lastHoverFired: {rowIdx: number; atomSerial: number | null; mode: string} | null = null;

  /** Emits CHEM_MOL3D_HOVER_EVENT for the 3D→2D hover bridge, with dedup.
   *  `atomSerial: null` means the cursor left all ligand atoms — the 2D
   *  side clears its preview. */
  private _fireMol3DHover(
    rowIdx: number, atomSerial: number | null, mode: 'preview' | 'paint' | 'erase',
  ): void {
    const last = this._lastHoverFired;
    if (last && last.rowIdx === rowIdx && last.atomSerial === atomSerial && last.mode === mode)
      return;
    this._lastHoverFired = {rowIdx, atomSerial, mode};
    const args: Mol3DHoverEventArgs = {
      mol3DColumnName: this.ligandColumnName!, rowIdx, atom3DSerial: atomSerial, mode,
    };
    grok.events.fireCustomEvent(CHEM_MOL3D_HOVER_EVENT, args);
  }

  /** Returns all Molstar structure components, or null if unavailable. */
  private _getAllComponents(): any[] | null {
    const plugin = this.viewer?.plugin;
    if (!plugin) return null;
    const structures = plugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) return null;
    const comps = structures.flatMap((s: any) => s.components ?? []);
    return comps.length > 0 ? comps : null;
  }
  /** Last received atom selection event args. Replayed after the viewer
   *  rebuilds its ligands so the highlight persists across cell switches
   *  (switching away destroys the viewer; switching back rebuilds it). */

  // -- Interactive atom highlighting (2D ↔ 3D bridge) ----------------------

  /** Replays highlights for all loaded ligands. Reads directly from the
   *  SMILES column's temp providers (the source of truth) so highlights
   *  appear automatically whenever Molstar rebuilds, without needing a
   *  fresh selection event. */
  private _replayHighlightIfCached(): void {
    this.logger.debug('[molstar-picker] scheduling replay after buildViewLigands');
    // Queue the replay through the viewSyncer so it runs AFTER all
    // pending structure rebuilds (rebuildViewCurrentRow, buildViewLigands)
    // have completed. This prevents the common issue where our overpaint
    // is applied to structures that are then immediately replaced.
    this.viewSyncer.sync('replayHighlight', async () => {
      this.logger.debug('[molstar-picker] replay firing now (synced)');
      try {
        const comps = this._getAllComponents();
        if (comps && this.viewer?.plugin)
          await clearStructureOverpaint(this.viewer.plugin, comps);
      } catch {
        /* best-effort */
      }
      await this._applyBaseColors();
      await this.highlightAllLigandAtoms();
    });
  }

  /** Applies distinct base colors to loaded ligand structures so the user
   *  can tell them apart. Called after ligands are built and on each
   *  highlight update. Does NOT require any atom selection — runs
   *  unconditionally when 2+ poses are loaded. */
  private async _applyBaseColors(): Promise<void> {
    const plugin = this.viewer?.plugin;
    if (!plugin) return;

    const loadedLigands: { rowIdx: number, structureRefs: string[] | null }[] = [];
    const seenRows = new Set<number>();
    const lig = this.ligands;
    if (lig) {
      const currentRow = lig.current?.rowIdx;
      if (currentRow != null && currentRow >= 0 && !seenRows.has(currentRow)) {
        loadedLigands.push(lig.current!);
        seenRows.add(currentRow);
      }
      const hoveredRow = lig.hovered?.rowIdx;
      if (hoveredRow != null && hoveredRow >= 0 && !seenRows.has(hoveredRow)) {
        loadedLigands.push(lig.hovered!);
        seenRows.add(hoveredRow);
      }
      if (lig.selected) {
        for (const sel of lig.selected) {
          if (sel.rowIdx >= 0 && !seenRows.has(sel.rowIdx)) {
            loadedLigands.push(sel);
            seenRows.add(sel.rowIdx);
          }
        }
      }
    }
    if (loadedLigands.length < 2) return;

    const allComponents = this._getAllComponents();
    if (!allComponents) return;

    // Determine current (selected) row — leave its pose in default colors.
    // Only color the OTHER (comparison) pose in teal.
    const currentRowIdx = this.dataFrame?.currentRowIdx ?? -1;

    for (const ligand of loadedLigands) {
      if (ligand.rowIdx === currentRowIdx) continue; // skip current row

      let structure: Structure | undefined;
      if (ligand.structureRefs && ligand.structureRefs.length >= 4) {
        const cell = plugin.state.data.cells.get(ligand.structureRefs[3]);
        if (cell?.obj?.data) structure = cell.obj.data;
      }
      if (!structure) continue;

      const targetStructure = structure;
      await setStructureOverpaint(
        plugin, allComponents, Color(0x616161),
        async (structureData: Structure) => {
          if (structureData !== targetStructure)
            return StructureElement.Loci.none(structureData);
          // Only color carbon atoms — leave N, O, S, etc. in their
          // default element colors so heteroatoms stay recognizable.
          const query = MolScriptBuilder.struct.generator.atomGroups({
            'atom-test': MolScriptBuilder.core.rel.eq([
              MolScriptBuilder.acp('elementSymbol'),
              MolScriptBuilder.es('C'),
            ]),
          });
          const sel = Script.getStructureSelection(query, structureData);
          return StructureSelection.toLociWithSourceUnits(sel);
        },
      );
    }
  }

  /** Reads the highlighted atom indices for a given row from the SMILES
   *  column's temp substruct providers. Returns empty array if none. */
  /** Reads highlighted atom indices for a given row from a specific
   *  molecule column's temp substruct providers. The column is identified
   *  by name (not searched by semType) to avoid ambiguity when a DataFrame
   *  has multiple molecule columns. */
  private _getHighlightedAtomsFromColumn(rowIdx: number, molColName?: string): number[] {
    const df = this.dataFrame;
    if (!df) return [];
    // Use the explicitly provided column name, or fall back to finding
    // the molecule column linked to this viewer's ligand column.
    let molCol: DG.Column | null = null;
    if (molColName)
      molCol = df.col(molColName);
    if (!molCol) {
      // Fall back: find the molecule column that has atom-picker providers.
      molCol = df.columns.toList().find(
        (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE &&
          ((c.temp as Record<string, unknown>)?.['chem-substruct-providers'] as unknown[] ?? [])
            .some((p: any) => p.__atomPicker)) ?? null;
    }
    if (!molCol) return [];
    // Read the atom-picker provider for this row from col.temp.
    const providers = ((molCol.temp as Record<string, unknown>)?.['chem-substruct-providers'] ?? []) as
      Array<{__atomPicker?: boolean; __rowIdx?: number; __atoms?: Set<number>}>;
    const picker = providers.find(
      (p) => p.__atomPicker && p.__rowIdx === rowIdx);
    if (!picker?.__atoms || picker.__atoms.size === 0) return [];
    return [...picker.__atoms];
  }

  /** Applies overpaint for ALL loaded ligands. Reads atom highlights
   *  directly from the SMILES column's providers (source of truth),
   *  falling back to the event cache. */
  public async highlightAllLigandAtoms(
    liveEvent?: { rowIdx: number, atoms: number[], mapping3D: AtomMapping3D | null } | null,
  ): Promise<void> {
    const plugin = this.viewer?.plugin;
    if (!plugin) return;
    // NOTE: no _highlightInProgress guard here — this method is the
    // single entry point for applying highlights and must always run,
    // even if a previous overpaint is mid-flight.

    this.logger.debug('[molstar-picker] highlightAllLigandAtoms called');

    // Collect all loaded ligands with their rows and structure refs.
    // Deduplicate by rowIdx — if current and hovered point to the same
    // row, only process it once (avoids applying the same mapping twice
    // and prevents stale highlights when the hovered row changes).
    type LigandInfo = { rowIdx: number, structureRefs: string[] | null };
    const loadedLigands: LigandInfo[] = [];
    const seenRows = new Set<number>();
    const lig = this.ligands;
    if (lig) {
      if (lig.current && lig.current.rowIdx >= 0 && !seenRows.has(lig.current.rowIdx)) {
        loadedLigands.push(lig.current);
        seenRows.add(lig.current.rowIdx);
      }
      if (lig.hovered && lig.hovered.rowIdx >= 0 && !seenRows.has(lig.hovered.rowIdx)) {
        loadedLigands.push(lig.hovered);
        seenRows.add(lig.hovered.rowIdx);
      }
      if (lig.selected) {
        for (const sel of lig.selected) {
          if (sel.rowIdx >= 0 && !seenRows.has(sel.rowIdx)) {
            loadedLigands.push(sel);
            seenRows.add(sel.rowIdx);
          }
        }
      }
    }

    this.logger.debug('[molstar-picker] loaded ligands:', () =>
      loadedLigands.map((l) => ({row: l.rowIdx, hasRefs: !!l.structureRefs})));

    // Collect serials per ligand row.
    const currentRowIdx = this.dataFrame?.currentRowIdx ?? -1;
    type LigandHighlight = {
      serials: number[], isCurrent: boolean, dataRef: string | null
    };
    const highlights: LigandHighlight[] = [];

    // Build cache key prefix from the current dataframe + molecule column.
    const df = this.dataFrame;
    const molCol = df?.columns.toList().find(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE);
    const dfId = df?.id ?? '';
    const dfName = df?.name ?? '';
    const colName = molCol?.name ?? '';
    const cache = getSelectionCache();

    for (const ligand of loadedLigands) {
      let atoms: number[] = [];
      let mapping3D: AtomMapping3D | null = null;

      if (liveEvent && liveEvent.rowIdx === ligand.rowIdx && liveEvent.atoms.length > 0) {
        atoms = liveEvent.atoms;
        mapping3D = liveEvent.mapping3D;
      }
      if (atoms.length === 0) {
        const key = selectionCacheKey(dfId, dfName, colName, ligand.rowIdx);
        const cached = cache.get(key);
        if (cached && cached.atoms.length > 0) {
          atoms = cached.atoms;
          mapping3D = cached.mapping3D ?? null;
        }
      }
      if (atoms.length === 0) {
        atoms = this._getHighlightedAtomsFromColumn(ligand.rowIdx, colName);
        const key = selectionCacheKey(dfId, dfName, colName, ligand.rowIdx);
        const cached = cache.get(key);
        mapping3D = cached?.mapping3D ?? null;
      }
      if (atoms.length === 0) continue;

      let structure: Structure | undefined;
      if (ligand.structureRefs && ligand.structureRefs.length >= 4) {
        const cell = plugin.state.data.cells.get(ligand.structureRefs[3]);
        if (cell?.obj?.data) structure = cell.obj.data;
      }

      const serials = computeSerials(atoms, mapping3D, structure);
      if (serials.length === 0) continue;

      // Store the data ref (first ref from addLigandOnStage) to scope
      // the overpaint to this ligand's subtree only.
      const dataRef = ligand.structureRefs?.[0] ?? null;
      highlights.push({
        serials,
        isCurrent: ligand.rowIdx === currentRowIdx,
        dataRef,
      });
    }

    const allComponents = this._getAllComponents();
    if (!allComponents) return;

    await clearStructureOverpaint(plugin, allComponents);
    await this._applyBaseColors();

    for (const h of highlights) {
      const color = h.isCurrent ? Color(0xFFFF00) : Color(0x00BCD4);
      // Apply overpaint to ALL components but scope the loci getter
      // to only match atoms in structures descended from this ligand's
      // data ref. This avoids the stale Structure identity issue while
      // still preventing cross-molecule contamination.
      await setStructureOverpaint(
        plugin, allComponents, color,
        async (structureData: Structure) => {
          // If we have a data ref, check if this structure descends
          // from it by walking up the state tree.
          if (h.dataRef) {
            let found = false;
            try {
              // Find the state cell for this structure and walk parents.
              for (const [, cell] of plugin.state.data.cells) {
                if (cell.obj?.data === structureData) {
                  let cur: any = cell;
                  for (let d = 0; cur && d < 6; d++) {
                    if (cur.transform?.ref === h.dataRef) {
                      found = true;
                      break;
                    }
                    const pRef = cur.transform?.parent;
                    cur = pRef ? plugin.state.data.cells.get(pRef) : null;
                  }
                  break;
                }
              }
            } catch {
              // skip check, allow fallback
              found = true;
            }
            if (!found) return StructureElement.Loci.none(structureData);
          }

          const atomSet = MolScriptBuilder.set(...h.serials);
          const query = MolScriptBuilder.struct.generator.atomGroups({
            'atom-test': MolScriptBuilder.core.set.has([
              atomSet,
              MolScriptBuilder.ammp('id'),
            ]),
          });
          const sel = Script.getStructureSelection(query, structureData);
          return StructureSelection.toLociWithSourceUnits(sel);
        },
      );
    }
  }

  /**
   * Highlights specific atoms in the 3D Molstar view using overpaint.
   *
   * @param atomIndices  0-based atom indices from the 2D picker.
   * @param mapping3D    Optional pre-computed mapping.
   * @param precomputedSerials  If true, atomIndices are already 1-based
   *                            serials (skip the 2D→3D mapping step).
   */
  public async highlightLigandAtoms(
    atomIndices: number[],
    mapping3D?: AtomMapping3D | null,
    precomputedSerials?: boolean,
  ): Promise<void> {
    const plugin = this.viewer?.plugin;
    if (!plugin) return;
    if (this._highlightInProgress) return; // skip concurrent calls
    this._highlightInProgress = true;

    try {
      // Find the current structure's components via the hierarchy manager.
      // These are the StructureComponentRef objects that setStructureOverpaint
      // needs to iterate and apply overpaint layers to.
      const structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures || structures.length === 0) return;

      // Collect ALL components from all structures so overpaint covers everything.
      const allComponents = structures.flatMap((s: any) => s.components ?? []);
      if (allComponents.length === 0) return;

      // Clear any previous overpaint first.
      await clearStructureOverpaint(plugin, allComponents);

      if (atomIndices.length === 0) return;

      let mol3DSerials: number[];
      if (precomputedSerials) {
        // Already 1-based serials — use directly.
        mol3DSerials = atomIndices;
      } else {
        // Compute serials from 2D atom indices.
        let structure: Structure | undefined;
        const currentLigand = this.ligands?.current;
        if (currentLigand?.structureRefs && currentLigand.structureRefs.length >= 4) {
          const structureRef = currentLigand.structureRefs[3];
          const structureCell = plugin.state.data.cells.get(structureRef);
          if (structureCell?.obj?.data)
            structure = structureCell.obj.data;
        }
        if (!structure && structures.length > 0)
          structure = structures[0]?.cell?.obj?.data;

        mol3DSerials = computeSerials(atomIndices, mapping3D, structure);
      }
      if (mol3DSerials.length === 0) return;

      _package.logger.debug(`[molstar-picker] overpaint serials [${mol3DSerials.slice(0, 10)}]`);

      const serialSet = mol3DSerials;

      // Apply yellow overpaint to matching atoms across all components.
      await setStructureOverpaint(
        plugin,
        allComponents,
        Color(0xFFFF00), // yellow
        async (structureData: Structure) => {
          // Build the query against the provided structure
          const atomSet = MolScriptBuilder.set(...serialSet);
          const query = MolScriptBuilder.struct.generator.atomGroups({
            'atom-test': MolScriptBuilder.core.set.has([
              atomSet,
              MolScriptBuilder.ammp('id'),
            ]),
          });
          const sel = Script.getStructureSelection(query, structureData);
          return StructureSelection.toLociWithSourceUnits(sel);
        },
      );
      this.logger.debug('[molstar-picker] overpaint applied');
    } catch (err: any) {
      this.logger.error(
        `highlightLigandAtoms failed: ${err?.message ?? err}`);
    } finally { this._highlightInProgress = false; }
  }

  /**
   * Clears any atom-level overpaint applied by {@link highlightLigandAtoms}.
   */
  public async clearLigandAtomHighlight(): Promise<void> {
    const plugin = this.viewer?.plugin;
    if (!plugin) return;
    try {
      const structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures || structures.length === 0) return;
      const allComponents = structures.flatMap((s: any) => s.components ?? []);
      if (allComponents.length === 0) return;
      await clearStructureOverpaint(plugin, allComponents);
    } catch { /* viewer may be disposed */ }
  }
}

/** margin indent */
function margin(indent: number): string {
  return '  '.repeat(indent);
}

