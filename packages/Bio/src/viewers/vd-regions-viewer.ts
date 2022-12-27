import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IVdRegionsViewer, PositionHeight, VdRegion, VdRegionType, WebLogoViewer} from '@datagrok-libraries/bio';

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


  public get df(): DG.DataFrame {
    return this.dataFrame;
  }

  // TODO: .onTableAttached is not calling on dataFrame set, onPropertyChanged  also not calling
  public async setDf(value: DG.DataFrame, regions: VdRegion[]) {
    console.debug('VdRegionsViewer.setDf()');
    await this.destroyView();
    this.regions = regions;
    this.dataFrame = value;
    await this.buildView();
  }

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
    // rxjs.fromEvent(this.root, 'mousemove').subscribe(this.onMouseMoveRoot.bind(this));
    this.root.addEventListener('mousemove', this.onMouseMoveRoot.bind(this));

    await this.buildView();
  }

  public override async onTableAttached() {
    await this.init();
  }

  public override async onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);
    if (property) {
      switch (property.name) {
      case 'regionTypes':
        break;
      case 'chains':
        break;
      case 'sequenceColumnNamePostfix':
        break;
      case 'skipEmptyPositions':
        // for (let orderI = 0; orderI < this.logos.length; orderI++) {
        //   for (let chainI = 0; chainI < this.chains.length; chainI++) {
        //     const chain: string = this.chains[chainI];
        //     this.logos[orderI][chain].setOptions({skipEmptyPositions: this.skipEmptyPositions});
        //   }
        // }
        // this.calcSize();
        await this.destroyView();
        await this.buildView();
        break;
      case 'positionWidth':
        await this.destroyView();
        await this.buildView();
        break;

      case 'positionHeight':
        await this.destroyView();
        await this.buildView();
        break;
      }
    }
  }


  public async reset() {

  }

  public async open(mlbView: DG.TableView) {
    if (!this.isOpened) {
      this.isOpened = true;
      this.panelNode = mlbView.dockManager.dock(this.root, DG.DOCK_TYPE.TOP, null, 'Regions', 0.2);
    }
  }

  public async show(mlbView: DG.TableView) {

  }

  // #region -- Handle controls' events --

  private resizing: boolean = false;

  private rootOnSizeChanged(args: any): void {
    this.calcSize();
  }

  // #endregion

  //#region -- View --
  private host: HTMLElement | null = null;
  private mainLayout: HTMLTableElement | null = null;
  private logos: { [chain: string]: WebLogoViewer }[] = [];

  private async destroyView(): Promise<void> {
    // TODO: Unsubscribe from and remove all view elements
    console.debug(`VdRegionsViewer.destroyView( mainLayout = ${!this.mainLayout ? 'none' : 'value'} )`);
    if (this.mainLayout != null) {
      // this.root.removeChild(this.host);
      this.mainLayout.remove();
      this.host!.remove();
      this.host = null;
      this.mainLayout = null;
    }
  }

  private async buildView(): Promise<void> {
    console.debug('VdRegionsViewer.buildView() start');

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

    const color: string = `#ffbb${Math.ceil(Math.random() * 255).toString(16)}`;
    this.host = ui.box(this.mainLayout,
      {/*style: {backgroundColor: color}*/});
    this.root.appendChild(this.host);
    this.root.style.overflowX = 'auto';

    this.calcSize();

    console.debug('VdRegionsViewer.buildView() end');
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

  //#endregion -- View --

  private onMouseMoveRoot(e: MouseEvent) {
    // ui.tooltip.show('text', e.x + 8, e.y + 8,);
    // console.log(`onMouseMoveRoot.( x: ${e.x}, y: ${e.y} )`);
  }
}
