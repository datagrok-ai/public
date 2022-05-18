import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

import {WebLogo} from '@datagrok-libraries/bio';


export enum VdRegionType {
  Unknown = 'unknown',
  FR = 'framework',
  CDR = 'cdr',
}

/** Describes V-DOMAIN (IG and TR) region (of multiple alignment)
 * https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
 */
export class VdRegion {
  type: VdRegionType;
  name: string;
  startPositionName: string;
  endPositionName: string;

  /**
   * start and position are strings because they correspond position names as column names in ANARCI output
   * @param {VdRegionType} type  Type of the region
   * @param {string} name  Name of the region
   * @param {string} startPositionName  Region start position (inclusive)
   * @param {string} endPositionName  Region end position (inclusive)
   */
  constructor(type: VdRegionType, name: string, startPositionName: string, endPositionName: string) {
    this.type = type;
    this.name = name;
    this.startPositionName = startPositionName;
    this.endPositionName = endPositionName;
  }

  public toString(): string {
    return `type: '${this.type}', name: '${this.name}', ` +
      `start: '${this.startPositionName}', end: '${this.endPositionName}'`;
  }
}

export class RegionViewer {
  private readonly region: VdRegion;

  //#region -- Design --
  private root: HTMLElement;

  //#endregion -- Design --

  constructor(region: VdRegion) {
    this.region = region;
  }

  public init() {
    const label: string = `start: ${this.region.startPositionName}, end: ${this.region.endPositionName}`;
    this.root = ui.div([ui.inlineText([label])]);
  }
}

/** Viewer with tabs based on description of chain regions.
 *  Used to define regions of an immunoglobulin LC.
 */
export class RegionSelectorViewer {
  private readonly regions: VdRegion[] = [
    new VdRegion(VdRegionType.FR, 'FR1', '1', '26'),
    new VdRegion(VdRegionType.CDR, 'CDR1', '27', '38'),
    new VdRegion(VdRegionType.FR, 'FR2', '39', '55'),
    new VdRegion(VdRegionType.CDR, 'CDR2', '56', '65'),
    new VdRegion(VdRegionType.FR, 'FR3', '66', '104'),
    new VdRegion(VdRegionType.CDR, 'CDR3', '105', '117'),
    new VdRegion(VdRegionType.FR, 'FR4', '118', '128'),
  ];

  private mlbView: DG.TableView;

  private regionsDf: DG.DataFrame;
  private regionsFg: DG.FilterGroup;
  // private regionsTV: DG.TableView;
  private regionsRoot: HTMLElement;

  private isOpened: boolean;
  private panelNode: DG.DockNode;

  //#region -- Design --
  private root: HTMLElement;

  private regionTabs: DG.TabControl;

  // private heavyDiv: HTMLDivElement;
  // private lightDiv: HTMLDivElement;

  private heavyWL: DG.Widget;
  //private heavyWL: DG.Widget;
  // private lightWL: WebLogo;

  //#endregion -- Design --

  constructor() { }

  public async init(mlbView: DG.TableView) {
    this.mlbView = mlbView;

    //#region regionsDF with filter
    this.regionsDf = DG.DataFrame.fromObjects(this.regions);
    this.regionsDf.rows.filter((row) => row.name == 'CDR1');
    // To available options /
    this.regionsFg = (await this.regionsDf.plot.fromType(DG.VIEWER.FILTERS, {
      // columnNames: ['name',],
      showFilterCountsIndication: false,
      showHeader: false,
      showSearchBox: false,
      filters: [
        {type: DG.FILTER_TYPE.CATEGORICAL, column: 'type', label: 'Region name', showHistogram: false},
        {type: DG.FILTER_TYPE.CATEGORICAL, column: 'name', label: 'Region type', showHistogram: false},
      ],
      title: 'Regions filter',
      showTitle: true,
      description: 'Filter for regions of multiple alignment by IMGT nomenclature',
      someProperty: 'Hello',
    })) as DG.FilterGroup;

    //#endregion regionsDF with filter

    this.heavyWL = await this.mlbView.dataFrame.plot.fromType('WebLogo', {sequenceColumnName: 'Heavy chain sequence'});
    this.heavyWL.root.style.background = 'window';
    // this.heavyWL.root.style.width = '100%';
    this.heavyWL.root.style.height = '100%';

    // this.heavyDiv = ui.div([ui.inlineText(['Heavy blue']), this.heavyWL.root], {style: {width: '20%;'}});
    // this.heavyDiv.style.background = 'blue';
    // this.heavyDiv.appendChild(this.heavyWL.root);

    // this.lightDiv = ui.div([ui.inlineText(['Light yellow'])], {style: {width: '20%;'}});
    // this.lightDiv.style.background = 'yellow';

    this.root = ui.box(/*[this.heavyDiv, this.lightDiv]*/);
    // this.root.style.background = '#FFEEEE';

    const rootNode: DG.DockNode = this.mlbView.dockManager.dock(this.root, DG.DOCK_TYPE.TOP);

    // this.root.appendChild(this.heavyDiv);
    // this.root.appendChild(this.lightDiv);

    this.mlbView.dockManager.dock(this.regionsFg.root, DG.DOCK_TYPE.LEFT, rootNode, 'Filter regions', 0.2);

    ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this));
    // rxjs.fromEvent(this.root, 'mousemove').subscribe(this.onMouseMoveRoot.bind(this));
    this.root.addEventListener('mousemove', this.onMouseMoveRoot.bind(this));

    await this.buildView();
    this.regionsDf.onRowsFiltered.subscribe(this.onRowsFilteredRegionsDf.bind(this));
  }

  public async reset() {

  }

  public async open(mlbView: DG.TableView) {
    if (!this.isOpened) {
      this.isOpened = true;
      this.panelNode = mlbView.dockManager.dock(this.root, DG.DOCK_TYPE.TOP, null, 'Regions', 0.2);
    }
  }

  public async close(mlbView: DG.TableView) {
    if (this.isOpened) {
      mlbView.dockManager.close(this.panelNode);
      this.isOpened = false;
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
  private view: HTMLElement[] = [];
  private mainLayout: HTMLTableElement = null;
  private logos: { [chain: string]: WebLogo }[] = [];

  private async destroyView(): Promise<void> {
    // TODO: Unsubscribe from and remove all view elements
    if (this.mainLayout != null) {
      this.root.removeChild(this.mainLayout);
      this.mainLayout = null;
    }
    console.log('RegionsSelectorViewer.destroyView()');
  }

  private async buildView(): Promise<void> {
    console.log('RegionSelectorViewer.buildView()');

    // const regionsIndices = this.regionsDf.filter.getSelectedIndexes();
    const chainList = ['Heavy', 'Light'];

    const regionsFiltered = Array.from(this.regionsDf.filter.getSelectedIndexes()).map(
      (rI) => { return this.regions[rI]; }, this);
    this.logos = [];

    for (let regionI = 0; regionI < regionsFiltered.length; regionI++) {
      this.logos[regionI] = {};
      for (const chain of chainList) {
        const region: VdRegion = regionsFiltered[regionI];
        this.logos[regionI][chain] = (await this.mlbView.dataFrame.plot.fromType('WebLogo', {
          sequenceColumnName: `${chain} chain sequence`,
          startPositionName: region.startPositionName,
          endPositionName: region.endPositionName,
        })) as WebLogo;
      }
    }

    // ui.tableFromMap()
    // DG.HtmlTable.create()
    this.mainLayout = ui.table(
      ['Heavy', 'Light'],
      (chain) => {
        const elements = [
          // This is chain label
          ...(regionsFiltered.length > 0 ? [ui.div(chain, {
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
          ...[...Array(regionsFiltered.length).keys()].map((rI) => {
            const wl: WebLogo = this.logos[rI][chain];
            wl.root.style.height = '100%';

            const resDiv = ui.div([wl.root]/*`${chain} ${regionsFiltered[rI]}`*/, {
              style: {
                height: '100%',
                marginTop: '4px',
                marginBottom: '4px',
              }
            });

            return resDiv;
          })];
        return elements;
      },
      [, ...regionsFiltered.map((r) => r.name)]);
    this.mainLayout.className = 'mlb-vd-regions-viewer-table2';
    // this.mainLayout.style.background = '#EEEEFF';
    // this.mainLayout.style.height = '100%';
    // this.mainLayout.style.border = '1px solid black';

    this.root.appendChild(this.mainLayout);
    this.root.style.overflowX = 'auto';

    for (let rI = 0; rI < this.logos.length; rI++) {
      await this.logos[rI]['Heavy'].init();
      await this.logos[rI]['Light'].init();
    }
    this.calcSize();
  }

  private calcSize() {
    const logoHeight = (this.root.clientHeight - 54) / 2;

    for (let rI = 0; rI < this.logos.length; rI++) {
      this.logos[rI]['Heavy'].root.style.height = `${logoHeight}px`;
      this.logos[rI]['Light'].root.style.height = `${logoHeight}px`;
    }
  }

  //#endregion -- View --

  private onMouseMoveRoot(e: MouseEvent) {
    // ui.tooltip.show('text', e.x + 8, e.y + 8,);
    // console.log(`onMouseMoveRoot.( x: ${e.x}, y: ${e.y} )`);
  }

  //#region -- Events handling of Regions filter --

  private async onRowsFilteredRegionsDf(): Promise<void> {
    console.debug('regionsDf.onRowsFiltered()');
    // Here we should destroy previous view and build new
    // if (this.view !== null) {
    //   for (const el of this.view) {
    //     console.debug('Destroy viewer');
    //   }
    // }

    await this.destroyView();
    await this.buildView();
  }

  //#endregion -- Events handling of Regions filter --
}
