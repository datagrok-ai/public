import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

// import {Viewer as rcsbViewer, ViewerProps as rcsbViewerProps} from '@rcsb/rcsb-molstar/src/viewer';
import {Viewer as RcsbViewer, ViewerProps as RcsbViewerProps} from '@rcsb/rcsb-molstar/build/src/viewer';
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
import {PluginCommands} from 'molstar/lib/mol-plugin/commands';
import {PluginContext} from 'molstar/lib/mol-plugin/context';
import {PluginLayoutControlsDisplay} from 'molstar/lib/mol-plugin/layout';
import {Color as msColor} from 'molstar/lib/mol-util/color';
import {getPdbHelper, IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {StructureComponentRef} from 'molstar/lib/mol-plugin-state/manager/structure/hierarchy-state';
import {BuiltInTrajectoryFormat, BuiltInTrajectoryFormats} from 'molstar/lib/mol-plugin-state/formats/trajectory';

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
   // -- Behaviour --
  showSelectedRowsLigands = 'showSelectedRowsLigands',
  showCurrentRowLigand = 'showCurrentRowLigand',
  showMouseOverRowLigand = 'showMouseOverRowLigand',
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

export type LigandMapItem = { rowIdx: number, structureRefs: Array<string> | null };
export type LigandMap = { selected: LigandMapItem[], current: LigandMapItem | null, hovered: LigandMapItem | null };

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

  // -- Behaviour --
  [PROPS.showSelectedRowsLigands]: boolean;
  [PROPS.showCurrentRowLigand]: boolean;
  [PROPS.showMouseOverRowLigand]: boolean;
  // propsEngine = new class {
  //
  // }(this);

  constructor() {
    super();
    this.helpUrl = '/help/visualize/viewers/biostructure';

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

    // -- Behaviour --
    this.showSelectedRowsLigands = this.bool(PROPS.showSelectedRowsLigands, false,
      {category: PROPS_CATS.BEHAVIOUR});
    this.showCurrentRowLigand = this.bool(PROPS.showCurrentRowLigand, true,
      {category: PROPS_CATS.BEHAVIOUR});
    this.showMouseOverRowLigand = this.bool(PROPS.showMouseOverRowLigand, true,
      {category: PROPS_CATS.BEHAVIOUR});

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
    const applyProperty = async (_propName: string, _value: any) => {
      if (!this.viewer) throw new Error('viewer does not exists');

      const plugin: PluginContext = this.viewer.plugin;
      switch (property.name) {
      case PROPS.layoutShowLog: {
        //plugin.layout.setProps({layoutShowLog: value;});

        //PluginCommands.Layout.Update(plugin, );
        const _k = 11;
      }
        break;
      case PROPS.layoutShowControls: {
        // eslint-disable-next-line new-cap
        await PluginCommands.Layout.Update(plugin, {state: {showControls: this.layoutShowControls}});
        const _k = 11;
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
    case PROPS.showSelectedRowsLigands:
    case PROPS.showCurrentRowLigand:
    case PROPS.showMouseOverRowLigand:
      this.rebuildViewLigands();
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

    this.viewPromise = this.viewPromise.then(async () => { // onTableAttached
      superOnTableAttached();
      await this.setData();
    });
  }

  override detach(): void {
    if (this.setDataInProgress) return;
    _package.logger.debug('MolstarViewer.detach(), ');

    const superDetach = super.detach.bind(this);
    this.detachPromise = this.detachPromise.then(async () => { // detach
      await this.viewPromise;
      if (this.viewed) {
        await this.destroyView('detach');
        this.viewed = false;
      }
      superDetach();
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
      if (!this.ligandColumnName) {
        const molCol: DG.Column | null = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
        if (molCol)
          this.ligandColumnName = molCol.name;
      }
    }).then(async () => {
      if (!this.viewed) {
        await this.buildView('setData');
        this.viewed = true;
      }
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
    const df: DG.DataFrame = this.dataFrame;

    await this.viewer.loadStructureFromData(this.pdbStr, 'pdb', false);
    const plugin: PluginContext = this.viewer.plugin;
    this.viewSubs.push(plugin.commands.subscribe(PluginCommands.Layout.Update, () => {

    }));

    this.viewSubs.push(df.onSelectionChanged.subscribe(this.dataFrameOnSelectionChanged.bind(this)));
    this.viewSubs.push(df.onCurrentRowChanged.subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
    this.viewSubs.push(df.onMouseOverRowChanged.subscribe(this.dataFrameOnMouseOverRowChanged.bind(this)));

    await this.buildViewLigands();
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
    fileEl.addEventListener('change', async (_event) => {
      const _k = 11;
      if (fileEl.files != null && fileEl.files.length == 1) {
        const [pdbStr, _pdbHelper]: [string, IPdbHelper] = await Promise.all([
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

  private rootOnSizeChanged(_value: any): void {
    _package.logger.debug('MolstarViewer.rootOnSizeChanged() ');
    this.calcSize();
  }

  private dataFrameOnSelectionChanged(_value: any): void {
    _package.logger.debug('BiostructureViewer.dataFrameOnSelectionChanged() ');
    this.rebuildViewLigands();
  }

  private dataFrameOnCurrentRowChanged(_value: any): void {
    _package.logger.debug('MolstarViewer.dataFrameOnCurrentRowChanged() ');
    this.rebuildViewLigands();
  }
  private dataFrameOnMouseOverRowChanged(_value: any): void {
    _package.logger.debug('BiostructureViewer.dataFrameOnMouseOverRowChanged() ');
    this.rebuildViewLigands();
  }

  // -- Ligands routines --

  private ligands: LigandMap = {selected: [], current: null, hovered: null};
  private ligandsPromise: Promise<void> = Promise.resolve();

  /** Unify get mol* component key/ref, not static for performance */
  getCompKey(comp: StructureComponentRef): string | null {
    return comp.cell.sourceRef ?? null; // comp.version
  }

  private getLigandStrOfRow(rowIdx: number): string {
    const ligandMol: string = this.dataFrame.get(this.ligandColumnName, rowIdx);
    const ligandStr: string = ligandMol + '$$$$';
    // const ligandBlob: Blob = new Blob([ligandStr], {type: 'text/plain'});
    // return ligandBlob;
    return ligandStr;
  }
  private rebuildViewLigands(): void {
    this.viewPromise = this.viewPromise.then(async () => {
      await this.destroyViewLigands();
      await this.buildViewLigands();
    });
  }

  private async destroyViewLigands(): Promise<void> {
    if (!this.viewer) throw new Error('The viewer is not created'); // return; // There is not PDB data
    if (!this.ligandColumnName) return;

    const allLigands: LigandMapItem[] = [
      ...this.ligands.selected,
      ...(this.ligands.current ? [this.ligands.current] : []),
      ...(this.ligands.hovered ? [this.ligands.hovered] : [])
    ];

    const refRemovingPromises: Promise<void>[] = [];
    for (const ligand of allLigands) {
      if (!ligand.structureRefs) continue;
      for (const lingandRef of ligand.structureRefs) {
        refRemovingPromises.push(this.viewer!.plugin.commands.dispatch(PluginCommands.State.RemoveObject,
          {state: this.viewer!.plugin.state.data,
            ref: lingandRef}));
      }
    }
    for (const ligand of allLigands) ligand.structureRefs = null; // unbind with this.stage.compList
  }

  /** Builds up ligands on the stage view */
  private async buildViewLigands(): Promise<void> {
    if (!this.viewer) throw new Error('The mol* viewer is not created'); // return; // There is not PDB data
    if (!this.ligandColumnName) return;

    this.ligands.selected = !this.showSelectedRowsLigands ? [] :
      wu(this.dataFrame.selection.getSelectedIndexes())
        .map((selRowIdx) => { return {rowIdx: selRowIdx, structureRefs: null}; })
        .toArray();
    this.ligands.current = !this.showCurrentRowLigand ? null :
      this.dataFrame.currentRowIdx >= 0 ? {rowIdx: this.dataFrame.currentRowIdx, structureRefs: null} : null;
    this.ligands.hovered = !this.showMouseOverRowLigand ? null :
      this.dataFrame.mouseOverRowIdx >= 0 ? {rowIdx: this.dataFrame.mouseOverRowIdx, structureRefs: null} : null;

    /** Adds ligand and returns component key */
    const addLigandOnStage = async (rowIdx: number, _color: DG.Color | null): Promise<Array<string> | null> => {
      const plugin = this.viewer!.plugin;
      const ligandLabel: string = `<Ligand at row ${rowIdx}>`;
      const ligandStr = this.getLigandStrOfRow(rowIdx);
      const _moldata = await plugin.builders.data.rawData({data: ligandStr, label: 'moldata'});
      const _moltraj = await plugin.builders.structure.parseTrajectory(
        _moldata, 'sdf');
      const _model = await plugin.builders.structure.createModel(_moltraj);
      const _structure = await plugin.builders.structure.createStructure(_model);
      const _component = await plugin.builders.structure.tryCreateComponentStatic(
        _structure, 'ligand', {label: ligandLabel}
      );
      await plugin.builders.structure.hierarchy.applyPreset(_moltraj, 'default',
        {representationPreset: 'polymer-and-ligand'});
      // this.ligandRefSet.add(_moldata.ref);
      // this.ligandRefSet.add(_moltraj.ref);
      // this.ligandRefSet.add(model.ref);
      // this.ligandRefSet.add(structure.ref);
      // this.ligandRefSet.add(a!.ref);
      return [_moldata.ref, _moltraj.ref, _model.ref, _structure.ref, _component!.ref];
    };

    const selCount = this.ligands.selected.length;
    for (const [selectedLigand, selI] of wu.enumerate(this.ligands.selected)) {
      const color =
        this.showCurrentRowLigand || this.showMouseOverRowLigand ?
          (selCount > 1 ? DG.Color.selectedRows : null) :
          (selCount > 1 ? DG.Color.scaleColor(selI, 0, selCount, 0.5) : null);

      selectedLigand.structureRefs = await addLigandOnStage(selectedLigand.rowIdx, color);
    }
    if (this.ligands.current) {
      const color = this.showSelectedRowsLigands ? DG.Color.currentRow : null;

      this.ligands.current.structureRefs = await addLigandOnStage(this.ligands.current.rowIdx, color);
    }
    if (this.ligands.hovered) {
      // TODO: color hovered ligand
      const color =
        this.showSelectedRowsLigands || this.showCurrentRowLigand ?
          DG.Color.mouseOverRows : null;
      this.ligands.hovered.structureRefs = await addLigandOnStage(this.ligands.hovered.rowIdx, color);
    }
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
  await initViewer()
    .then((v: RcsbViewer) => {
      v.loadPdbId(pdbId);
    });
  //v.handleResize();
}

/**
 * Creates an instance of Mol* viewer.
 *
 * @param {string} data Data in PDB
 */
export async function byData(data: string, name: string = 'Mol*', format: BuiltInTrajectoryFormat = 'pdb') {
  await initViewer(name)
    .then(async (viewer: RcsbViewer) => {
      // detecting format by data content
      let binary: boolean = false;
      if (data.includes('mmcif_pdbx.dic')) {
        format = 'mmcif';
        binary = false;
      }

      await viewer.loadStructureFromData(data, format, binary);
    });
  //v.handleResize();
}

export async function viewMolstarUI(content: string, name?: string, format?: BuiltInTrajectoryFormat): Promise<void> {
  await byData(content, name, format);
}

/** Creates view with Molstar viewer to preview Biostructure (PDB)
 * returns the view immidiately, but the viewer is created asynchronously and promise
 * for that is returned separately which resolves once the viewer is initialized.
*/
export function previewMolstarUI(file: DG.FileInfo): { view: DG.View, loadingPromise: Promise<void> } {
  const builtinFormats = BuiltInTrajectoryFormats.map((obj) => obj[0]) as string[];
  const extendedFormats = ['cif', 'mcif'];
  if (!isSupportedFormat()) {
    grok.shell.error(`Unsupported format: ${file.extension}`);
    throw new Error(`Unsupported format: ${file.extension}`);
  }

  const formatLoader = extendedFormats.includes(file.extension) ? 'mmcif' : file.extension;

  const view = DG.View.create();
  const viewer = new RcsbViewer(view.root, castProps(defaults));
  const subs: Unsubscribable[] = [];
  subs.push(ui.onSizeChanged(view.root).subscribe((_value: any) => {
    viewer.handleResize();
  }));
  subs.push(grok.events.onViewRemoved.subscribe((evtView) => {
    if (evtView.id === view.id)
      for (const sub of subs) sub.unsubscribe();
  }));

  function isSupportedFormat() {
    return [...builtinFormats, ...extendedFormats].includes(file.extension);
  }
  async function loadString(data: string) {
    const binary: boolean = false;
    await viewer.loadStructureFromData(data, formatLoader as BuiltInTrajectoryFormat, binary);
  }

  async function loadBytes(bytes: any) {
    const binary: boolean = false;
    await viewer.loadStructureFromData(bytes, formatLoader as BuiltInTrajectoryFormat, binary);
  }

  const loadingPromise = new Promise<void>(async (resolve, reject) => {
    try {
      if (isSupportedFormat())
        await loadString(await file.readAsString());
      else
        await loadBytes(await file.readAsBytes());
      resolve();
    } catch (e) {
      reject(e);
    }
  });

  return {view, loadingPromise};
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
