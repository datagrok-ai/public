import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import {FilterSources, WebLogoViewer, PROPS as wlPROPS} from '../viewers/web-logo-viewer';
import {IVdRegionsViewer, VdRegion, VdRegionType} from '@datagrok-libraries/bio/src/vd-regions';
import {PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {Unsubscribable} from 'rxjs';

const vrt = VdRegionType;

// Positions of regions for numbering schemes
// http://www.bioinf.org.uk/abs/info.html

// const imgtRegions: VdRegion[] = [
//   new VdRegion(vrt.FR, 'FR1', 'Light', 1, '1', '26'),
//   new VdRegion(vrt.FR, 'FR1', 'Heavy', 1, '1', '26'),
//
//   new VdRegion(vrt.CDR, 'CDR1', 'Light', 2, '27', '38'), // 27-32
//   new VdRegion(vrt.CDR, 'CDR1', 'Heavy', 2, '27', '38'), // 27-32
//
//   new VdRegion(vrt.FR, 'FR2', 'Light', 3, '39', '55'),
//   new VdRegion(vrt.FR, 'FR2', 'Heavy', 3, '39', '55'),
//
//   new VdRegion(vrt.CDR, 'CDR2', 'Light', 4, '56', '65'),
//   new VdRegion(vrt.CDR, 'CDR2', 'Heavy', 4, '56', '65'),
//
//   new VdRegion(vrt.FR, 'FR3', 'Light', 5, '66', '104'),
//   new VdRegion(vrt.FR, 'FR3', 'Heavy', 5, '66', '104'),
//
//   new VdRegion(vrt.CDR, 'CDR3', 'Light', 6, '105', '117'),
//   new VdRegion(vrt.CDR, 'CDR3', 'Heavy', 6, '105', '117'),
//
//   new VdRegion(vrt.FR, 'FR4', 'Light', 7, '118', null/*127*/),
//   new VdRegion(vrt.FR, 'FR4', 'Heavy', 7, '118', null/*128*/),
// ];

/** Viewer with tabs based on description of chain regions.
 *  Used to define regions of an immunoglobulin LC.
 */
export class VdRegionsViewer extends DG.JsViewer implements IVdRegionsViewer {
  private viewed: boolean = false;

  // private regionsDf: DG.DataFrame;
  private regionsFg: DG.FilterGroup | null = null;
  // private regionsTV: DG.TableView;
  private regionsRoot: HTMLElement | null = null;

  private isOpened: boolean = false;
  private panelNode: DG.DockNode | null = null;

  public regions: VdRegion[] = [];
  public regionTypes: string[];
  public chains: string[];
  public sequenceColumnNamePostfix: string;

  public skipEmptyPositions: boolean;
  public positionWidth: number;
  public positionHeight: string;

  constructor() {
    super();

    // To prevent ambiguous numbering scheme in MLB
    this.regionTypes = this.stringList('regionTypes', [vrt.CDR],
      {choices: Object.values(vrt).filter((t) => t != vrt.Unknown)});
    this.chains = this.stringList('chains', ['Heavy', 'Light'],
      {choices: ['Heavy', 'Light']});
    this.sequenceColumnNamePostfix = this.string('sequenceColumnNamePostfix', 'chain sequence');

    this.skipEmptyPositions = this.bool('skipEmptyPositions', false);
    this.positionWidth = this.float('positionWidth', 16);
    this.positionHeight = this.string('positionHeight', PositionHeight.Entropy,
      {choices: Object.keys(PositionHeight)});
  }

  public async init() {
    //#region regionsDF with filter
    // this.regionsDf = DG.DataFrame.fromObjects(this.regions);
    // this.regionsDf.rows.filter((row) => row.name == 'CDR1');
    // // To available options /
    // this.regionsFg = (await this.regionsDf.plot.fromType(DG.VIEWER.FILTERS, {
    //   // columnNames: ['name',],
    //   showFilterCountsIndication: false,
    //   showHeader: false,
    //   showSearchBox: false,
    //   filters: [
    //     {type: DG.FILTER_TYPE.CATEGORICAL, column: 'type', label: 'Region name', showHistogram: false},
    //     {type: DG.FILTER_TYPE.CATEGORICAL, column: 'name', label: 'Region type', showHistogram: false},
    //   ],
    //   title: 'Regions filter',
    //   showTitle: true,
    //   description: 'Filter for regions of multiple alignment by IMGT nomenclature',
    //   someProperty: 'Hello',
    // })) as DG.FilterGroup;

    //#endregion regionsDF with filter

    // this.mlbView.dockManager.dock(this.regionsFg.root, DG.DOCK_TYPE.LEFT, rootNode, 'Filter regions', 0.2);

    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
    this.subs.push(rxjs.fromEvent<MouseEvent>(this.root, 'mousemove').subscribe(this.rootOnMouseMove.bind(this)));

    // await this.buildView('init'); // init
  }

  public override async onTableAttached() {
    const superOnTableAttached = super.onTableAttached.bind(this);
    this.viewPromise = this.viewPromise.then(async () => {
      superOnTableAttached();
      if (!this.viewed) {
        await this.buildView('onTableAttached'); // onTableAttached
        this.viewed = true;
      }
    });
  }

  public override onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('Bio: VdRegionsViewer.onPropertyChanged() property is null');
      return;
    }

    if (property) {
      switch (property.name) {
      case 'regionTypes':
        break;
      case 'chains':
        break;
      case 'sequenceColumnNamePostfix':
        break;
        // for (let orderI = 0; orderI < this.logos.length; orderI++) {
        //   for (let chainI = 0; chainI < this.chains.length; chainI++) {
        //     const chain: string = this.chains[chainI];
        //     this.logos[orderI][chain].setOptions({skipEmptyPositions: this.skipEmptyPositions});
        //   }
        // }
        // this.calcSize();
      }
    }

    switch (property.name) {
    case 'skipEmptyPositions':
    case 'positionWidth':
    case 'positionHeight':
      window.setTimeout(
        async () => {
          await this.setData(this.dataFrame, this.regions);
        }, 0 /* next event cycle */);
      break;
    }
  }

  // -- Data --

  // TODO: .onTableAttached is not calling on dataFrame set, onPropertyChanged  also not calling
  public async setData(mlbDf: DG.DataFrame, regions: VdRegion[]) {
    console.debug('Bio: VdRegionsViewer.setData()');
    this.viewPromise = this.viewPromise.then(async () => {
      if (this.viewed) {
        await this.destroyView('setData'); // setData
        this.viewed = false;
      }
    });

    this.regions = regions;
    this.dataFrame = mlbDf;

    this.viewPromise = this.viewPromise.then(async () => {
      if (!this.viewed) {
        await this.buildView('setData'); // setData
        this.viewed = true;
      }
    });
  }

  override detach() {
    const superDetach = super.detach.bind(this);
    // window.setTimeout(async () => {
    //   if (this.viewed) {
    //     await this.destroyView('detach'); // detach
    //     this.viewed = false;
    //   }
    //   superDetach();
    // }, 0 /* next event cycle */);
    this.viewPromise = this.viewPromise.then(async () => {
      if (this.viewed) {
        await this.destroyView('detach'); // detach
        this.viewed = false;
      }
      superDetach();
    });
  }

  // -- View --

  private viewPromise: Promise<void> = Promise.resolve();

  private host: HTMLElement | null = null;
  private filterSourceInput: DG.InputBase<boolean | null> | null = null;
  private mainLayout: HTMLTableElement | null = null;
  private logos: { [chain: string]: WebLogoViewer }[] = [];

  private viewSubs: Unsubscribable[] = [];

  private async destroyView(purpose: string): Promise<void> {
    // TODO: Unsubscribe from and remove all view elements
    console.debug(`Bio: VdRegionsViewer.destroyView( mainLayout = ${!this.mainLayout ? 'none' : 'value'} ), ` +
      `purpose = '${purpose}'`);
    if (this.filterSourceInput) {
      //
      ui.empty(this.filterSourceInput.root);
    }

    if (this.mainLayout != null) {
      // this.root.removeChild(this.host);
      this.mainLayout.remove();
      this.host!.remove();
      this.host = null;
      this.mainLayout = null;
    }

    for (const sub of this.viewSubs) sub.unsubscribe();
  }

  private async buildView(purpose: string): Promise<void> {
    console.debug(`Bio: VdRegionsViewer.buildView() begin, ` + `purpose = '${purpose}'`);

    const colNames: { [chain: string]: string } = Object.assign({},
      ...this.chains.map((chain) => ({[chain]: `${chain} ${this.sequenceColumnNamePostfix}`})));

    const regionsFiltered: VdRegion[] = this.regions.filter((r: VdRegion) => this.regionTypes.includes(r.type));

    const orderList: number[] = Array.from(new Set(regionsFiltered.map((r) => r.order))).sort();

    this.logos = [];

    for (let orderI = 0; orderI < orderList.length; orderI++) {
      const regionChains: { [chain: string]: WebLogoViewer } = {};
      for (const chain of this.chains) {
        const region: VdRegion | undefined = regionsFiltered
          .find((r) => r.order == orderList[orderI] && r.chain == chain);
        regionChains[chain] = (await this.dataFrame.plot.fromType('WebLogo', {
          sequenceColumnName: colNames[chain],
          startPositionName: region!.positionStartName,
          endPositionName: region!.positionEndName,
          fixWidth: true,
          skipEmptyPositions: this.skipEmptyPositions,
          positionWidth: this.positionWidth,
          positionHeight: this.positionHeight,
        })) as unknown as WebLogoViewer;
      }
      // WebLogo creation fires onRootSizeChanged event even before control being added to this.logos
      this.logos[orderI] = regionChains;
    }

    // ui.tableFromMap()
    // DG.HtmlTable.create()
    this.mainLayout = ui.table(
      this.chains,
      (chain) => {
        const elements = [
          // This is chain label
          ...(orderList.length > 0 ? [ui.div(chain, {
            style: {
              transform: 'rotate(-90deg)',
              font: '12px Roboto, Roboto Local, sans-serif',
              textAlign: 'center',
              width: '16px',
              marginTop: '24px',
              marginLeft: '6px',
            }
          })] : []),
          // List with controls for regions
          ...[...Array(orderList.length).keys()].map((orderI) => {
            const wl: WebLogoViewer = this.logos[orderI][chain];
            wl.root.style.height = '100%';

            const resDiv = ui.div([wl.root]/*`${chain} ${regionsFiltered[rI]}`*/, {
              style: {
                // height: '100%',
                marginTop: '4px',
                marginBottom: '4px',
              }
            });

            return resDiv;
          })];
        return elements;
      },
      ['', ...[...Array(orderList.length).keys()].map(
        (orderI: number) => regionsFiltered.find(
          (r: VdRegion) => r.order == orderList[orderI] && r.chain == this.chains[0]
        )!.name || 'Name')]
    );
    this.mainLayout.className = 'mlb-vd-regions-viewer-table2';
    // this.mainLayout.style.background = '#EEEEFF';
    // this.mainLayout.style.height = '100%';
    // this.mainLayout.style.border = '1px solid black';

    this.filterSourceInput = ui.boolInput('', false, this.filterSourceInputOnValueChanged.bind(this));
    this.filterSourceInput.root.style.position = 'absolute';
    this.filterSourceInput.root.style.left = '10px';
    this.filterSourceInput.root.style.top = '-3px';
    ui.tooltip.bind(this.filterSourceInput.root, 'Check to filter sequences for selected VRs');

    const color: string = `#ffbb${Math.ceil(Math.random() * 255).toString(16)}`;
    this.host = ui.div([this.mainLayout, this.filterSourceInput!.root],
      {/*style: {backgroundColor: color}*/});
    this.root.appendChild(this.host);
    this.root.style.overflowX = 'auto';

    this.calcSize();

    console.debug('Bio: VdRegionsViewer.buildView() end');
  }

  private calcSize() {
    const logoHeight = (this.root.clientHeight - 54) / this.chains.length;

    const maxHeight: number = Math.min(logoHeight,
      Math.max(...this.logos.map((wlDict) =>
        Math.max(...Object.values(wlDict).map((wl) => wl.maxHeight))))
    );

    for (let orderI = 0; orderI < this.logos.length; orderI++) {
      for (let chainI = 0; chainI < this.chains.length; chainI++) {
        const chain: string = this.chains[chainI];
        this.logos[orderI][chain].root.style.height = `${maxHeight}px`;
      }
    }
  }

  // -- Handle events --

  private rootOnSizeChanged(args: any): void {
    this.calcSize();
  }

  private rootOnMouseMove(e: MouseEvent) {
    // ui.tooltip.show('text', e.x + 8, e.y + 8,);
    // console.log(`onMouseMoveRoot.( x: ${e.x}, y: ${e.y} )`);
  }

  private filterSourceInputOnValueChanged(): void {
    const filterSource: FilterSources = this.filterSourceInput!.value == true ?
      FilterSources.Selected : FilterSources.Filtered;

    for (let orderI = 0; orderI < this.logos.length; orderI++) {
      for (let chainI = 0; chainI < this.chains.length; chainI++) {
        const chain: string = this.chains[chainI];
        const wl: DG.JsViewer = this.logos[orderI][chain];
        wl.setOptions({[wlPROPS.filterSource]: filterSource});
      }
    }
  }
}
