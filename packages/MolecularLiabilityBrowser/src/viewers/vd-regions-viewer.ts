import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

import {VdRegionType, VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';

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
//
// const chothiaRegions: VdRegion[] = [
//   new VdRegion(vrt.FR, 'FR1', 'Light', 1, '1', '23'),
//   new VdRegion(vrt.FR, 'FR1', 'Heavy', 1, '1', '25'),
//
//   new VdRegion(vrt.CDR, 'L1', 'Light', 2, '24', '34'),
//   new VdRegion(vrt.CDR, 'H1', 'Heavy', 2, '26', '32'),
//
//   new VdRegion(vrt.FR, 'FR2', 'Light', 3, '35', '49'),
//   new VdRegion(vrt.FR, 'FR2', 'Heavy', 3, '33', '51'),
//
//   new VdRegion(vrt.CDR, 'L2', 'Light', 4, '50', '56'),
//   new VdRegion(vrt.CDR, 'H2', 'Heavy', 4, '52', '56'),
//
//   new VdRegion(vrt.FR, 'FR3', 'Light', 5, '57', '88'),
//   new VdRegion(vrt.FR, 'FR3', 'Heavy', 5, '57', '94'),
//
//   new VdRegion(vrt.CDR, 'L3', 'Light', 6, '89', '97'),
//   new VdRegion(vrt.CDR, 'H3', 'Heavy', 6, '95', '102'),
//
//   new VdRegion(vrt.FR, 'FR4', 'Light', 7, '98', null/*107*/),
//   new VdRegion(vrt.FR, 'FR4', 'Heavy', 7, '103', null/*113*/),
// ];
//
// const kabatRegions: VdRegion[] = [
//   new VdRegion(vrt.FR, 'FR1', 'Light', 1, '1', '23'),
//   new VdRegion(vrt.FR, 'FR1', 'Heavy', 1, '1', '37'),
//
//   new VdRegion(vrt.CDR, 'L1', 'Light', 2, '24', '42'),
//   new VdRegion(vrt.CDR, 'H1', 'Heavy', 2, '38', '42'),
//
//   new VdRegion(vrt.FR, 'FR2', 'Light', 3, '43', '56'),
//   new VdRegion(vrt.FR, 'FR2', 'Heavy', 3, '43', '55'),
//
//   new VdRegion(vrt.CDR, 'L2', 'Light', 4, '57', '72'),
//   new VdRegion(vrt.CDR, 'H2', 'Heavy', 4, '56', '74'),
//
//   new VdRegion(vrt.FR, 'FR3', 'Light', 5, '73', '106'),
//   new VdRegion(vrt.FR, 'FR3', 'Heavy', 5, '75', '106'/*??*/),
//
//   new VdRegion(vrt.CDR, 'L3', 'Light', 6, '107', '138'),
//   new VdRegion(vrt.CDR, 'H3', 'Heavy', 6, '107' /*??*/, '138'),
//
//   new VdRegion(vrt.FR, 'FR4', 'Light', 7, '139', null /*107*/),
//   new VdRegion(vrt.FR, 'FR4', 'Heavy', 7, '139', null /*113*/),
// ];

// const numberingSchemeRegions: { [name: string]: VdRegion[] } = {
//   'imgt': imgtRegions,
//   'chothia': chothiaRegions,
//   'kabat': kabatRegions,
// };

/** Viewer with tabs based on description of chain regions.
 *  Used to define regions of an immunoglobulin LC.
 */
export class VdRegionsViewer extends DG.JsViewer {
  // private readonly regions: VdRegion[] = imgtRegions;

  // private regionsDf: DG.DataFrame;
  private regionsFg: DG.FilterGroup;
  // private regionsTV: DG.TableView;
  private regionsRoot: HTMLElement;

  private isOpened: boolean;
  private panelNode: DG.DockNode;

  public regions: VdRegion[] = [];
  // public numberingScheme: string;
  public regionTypes: string[];
  public chains: string[];
  public heavyChainSequenceColumnName: string;
  public lightChainSequenceColumnName: string;

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
    // this.numberingScheme = this.string('numberingScheme', 'imgt',
    //   {choices: Object.keys(numberingSchemeRegions)});
    this.regionTypes = this.stringList('regionTypes', [vrt.CDR],
      {choices: Object.values(vrt).filter((t) => t != vrt.Unknown)});
    this.chains = this.stringList('chains', ['Heavy', 'Light'],
      {choices: ['Heavy', 'Light']});
    this.heavyChainSequenceColumnName = this.string('heavyChainSequenceColumnName', null);
    this.lightChainSequenceColumnName = this.string('lightChainSequenceColumnName', null);
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

    // this.root = ui.box(/*[this.heavyDiv, this.lightDiv]*/);
    // this.root.style.background = '#FFEEEE';

    // const rootNode: DG.DockNode = this.mlbView.dockManager.dock(this.root, DG.DOCK_TYPE.DOWN);

    // this.root.appendChild(this.heavyDiv);
    // this.root.appendChild(this.lightDiv);

    // this.mlbView.dockManager.dock(this.regionsFg.root, DG.DOCK_TYPE.LEFT, rootNode, 'Filter regions', 0.2);

    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
    // rxjs.fromEvent(this.root, 'mousemove').subscribe(this.onMouseMoveRoot.bind(this));
    this.root.addEventListener('mousemove', this.onMouseMoveRoot.bind(this));

    await this.buildView();
    // this.regionsDf.onRowsFiltered.subscribe(this.onRowsFilteredRegionsDf.bind(this));
  }

  public override async onTableAttached() {
    await this.init();
  }

  public override async onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);
    if (property) {
      switch (property.name) {
      case 'numberingScheme':
        // Ignore change property numberingScheme to prevent buildView() double async call
        // await this.destroyView();
        // await this.buildView();
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
  private host: HTMLElement = null;
  private mainLayout: HTMLTableElement = null;
  private logos: { [chain: string]: WebLogo }[] = [];

  private async destroyView(): Promise<void> {
    // TODO: Unsubscribe from and remove all view elements
    console.debug(`VdRegionsViewer.destroyView( mainLayout = ${!this.mainLayout ? 'none' : 'value'} )`);
    if (this.mainLayout != null) {
      // this.root.removeChild(this.host);
      this.mainLayout.remove();
      this.host.remove();
      this.host = null;
      this.mainLayout = null;
    }
  }

  private async buildView(): Promise<void> {
    console.debug('VdRegionsViewer.buildView() start');

    const colNames: { [chain: string]: string } = Object.assign({},
      ...['Heavy', 'Light'].map((chain) => ({[chain]: `${chain} chain sequence`})));

    // const scheme = this.numberingScheme ||
    //   this.dataFrame.col(colNames['Heavy']).getTag('numberingScheme');

    // const regionsFiltered: VdRegion[] = scheme in numberingSchemeRegions ?
    //   numberingSchemeRegions[scheme].filter((r: VdRegion) => this.regionTypes.includes(r.type)) : [];
    const regionsFiltered: VdRegion[] = this.regions.filter((r: VdRegion) => this.regionTypes.includes(r.type));

    const orderList: number[] = Array.from(new Set(regionsFiltered.map((r) => r.order))).sort();

    this.logos = [];

    for (let orderI = 0; orderI < orderList.length; orderI++) {
      const regionChains: { [chain: string]: WebLogo } = {};
      for (const chain of this.chains) {
        const region: VdRegion = regionsFiltered.find((r) => r.order == orderList[orderI] && r.chain == chain);

        regionChains[chain] = (await this.dataFrame.plot.fromType('WebLogo', {
          sequenceColumnName: colNames[chain],
          startPositionName: region.positionStartName,
          endPositionName: region.positionEndName,
          fixWidth: true,
        })) as unknown as WebLogo;
      }
      // WebLogo creation fires onRootSizeChanged event even before control being added to this.logos
      this.logos[orderI] = regionChains;
    }

    // ui.tableFromMap()
    // DG.HtmlTable.create()
    this.mainLayout = ui.table(
      ['Heavy', 'Light'],
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
            const wl: WebLogo = this.logos[orderI][chain];
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
      [null, ...[...Array(orderList.length).keys()].map(
        (orderI: number) => regionsFiltered.find(
          (r: VdRegion) => r.order == orderList[orderI] && r.chain == 'Heavy'
        ).name || 'Name')]
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
    const logoHeight = (this.root.clientHeight - 54) / 2;

    const maxHeight: number = Math.min(logoHeight,
      Math.max(...this.logos.map((wlDict) =>
        Math.max(...Object.values(wlDict).map((wl) => wl.maxHeight))))
    );

    for (let orderI = 0; orderI < this.logos.length; orderI++) {
      this.logos[orderI]['Heavy'].root.style.height = `${maxHeight}px`;
      this.logos[orderI]['Light'].root.style.height = `${maxHeight}px`;
    }
  }

  //#endregion -- View --

  private onMouseMoveRoot(e: MouseEvent) {
    // ui.tooltip.show('text', e.x + 8, e.y + 8,);
    // console.log(`onMouseMoveRoot.( x: ${e.x}, y: ${e.y} )`);
  }
}
