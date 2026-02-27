import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {fromEvent, Observable, Subject, Unsubscribable} from 'rxjs';

import {testEvent} from '@datagrok-libraries/utils/src/test';
import {
  IVdRegionsViewer,
  VdRegion, VdRegionType,
  VdRegionsProps, VdRegionsPropsDefault,
} from '@datagrok-libraries/bio/src/viewers/vd-regions';
import {PromiseSyncer} from '@datagrok-libraries/bio/src/utils/syncer';
import {FilterSources, IWebLogoViewer, PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';


import {WebLogoViewer, PROPS as wlPROPS} from '../viewers/web-logo-viewer';

import {_package} from '../package';

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

export enum PROPS_CATS {
  STYLE = 'Style',
  BEHAVIOR = 'Behavior',
  LAYOUT = 'Layout',
  DATA = 'Data',
}

export enum PROPS {
  // -- Data --
  skipEmptyPositions = 'skipEmptyPositions',
  regionTypes = 'regionTypes',
  chains = 'chains',

  // -- Layout --
  fitWidth = 'fitWidth',
  positionWidth = 'positionWidth',
  positionHeight = 'positionHeight',

  // -- Behavior --
  filterSource = 'filterSource',
}

const defaults: VdRegionsProps = VdRegionsPropsDefault;

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
  public regionTypes: VdRegionType[];
  public chains: string[];
  // public sequenceColumnNamePostfix: string;

  public skipEmptyPositions: boolean;
  public fitWidth: boolean;
  public positionWidth: number;
  public positionHeight: PositionHeight;

  public filterSource: FilterSources;

  constructor() {
    super();

    // -- Data --
    this.skipEmptyPositions = this.bool(PROPS.skipEmptyPositions, defaults.skipEmptyPositions,
      {category: PROPS_CATS.DATA});

    // To prevent ambiguous numbering scheme in MLB
    this.regionTypes = this.stringList(PROPS.regionTypes, defaults.regionTypes, {
      category: PROPS_CATS.DATA, choices: Object.values(vrt).filter((t) => t != vrt.Unknown)
    }) as VdRegionType[];
    this.chains = this.stringList(PROPS.chains, defaults.chains,
      {category: PROPS_CATS.DATA, choices: ['Heavy', 'Light']});

    // -- Layout --
    this.fitWidth = this.bool(PROPS.fitWidth, defaults.fitWidth,
      {category: PROPS_CATS.LAYOUT});
    this.positionWidth = this.float(PROPS.positionWidth, defaults.positionWidth, {
      category: PROPS_CATS.LAYOUT, editor: 'slider', min: 0, max: 64,
      description: 'Internal WebLogo viewers property width of position.'
    });
    this.positionHeight = this.string(PROPS.positionHeight, defaults.positionHeight,
      {category: PROPS_CATS.LAYOUT, choices: Object.keys(PositionHeight)}) as PositionHeight;

    // -- Behavior --
    this.filterSource = this.string(PROPS.filterSource, defaults.filterSource,
      {category: PROPS_CATS.BEHAVIOR, choices: Object.values(FilterSources)}) as FilterSources;

    this.viewSyncer = new PromiseSyncer(_package.logger);
  }

  private static viewerCounter: number = -1;
  private readonly viewerId: number = ++VdRegionsViewer.viewerCounter;

  private viewerToLog(): string { return `VdRegionsViewer<${this.viewerId}>`; }

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

    this.subs.push(fromEvent<MouseEvent>(this.root, 'mousemove').subscribe(this.rootOnMouseMove.bind(this)));

    // await this.buildView('init'); // init
  }

  override detach() {
    const logPrefix = `${this.viewerToLog()}.detach()`;
    const superDetach = super.detach.bind(this);
    this.viewSyncer.sync(`${logPrefix}`, async () => { // detach
      if (this.setDataInProgress) return; // check setDataInProgress synced
      if (this.viewed) {
        await this.destroyView('detach');
        this.viewed = false;
      }
      superDetach();
    });
  }

  override onTableAttached() {
    super.onTableAttached();
    this.setData(this.regions);
  }

  public override onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('Bio: VdRegionsViewer.onPropertyChanged() property is null');
      return;
    }

    switch (property.name) {
      case PROPS.regionTypes:
      case PROPS.chains:
        this.setData(this.regions);
        break;
    }

    switch (property.name) {
      case PROPS.skipEmptyPositions:
        for (let orderI = 0; orderI < this.logos.length; ++orderI) {
          for (const chain of this.chains)
            this.logos[orderI][chain].setOptions({[wlPROPS.skipEmptyPositions]: this.skipEmptyPositions});
        }
        this.calcSize();
        break;

      case PROPS.fitWidth:
      case PROPS.positionWidth:
        this.calcSize();
        break;

      case PROPS.positionHeight:
        for (let orderI = 0; orderI < this.logos.length; ++orderI) {
          for (const chain of this.chains)
            this.logos[orderI][chain].setOptions({[wlPROPS.positionHeight]: this.positionHeight});
        }
        this.calcSize();
        break;

      case PROPS.filterSource:
        this.filterSourceInput.value = this.filterSource;
        break;

      default:
        this.setData(this.regions); // onPropertyChanged
        break;
    }
  }

  // -- Data --

  // private static viewerCount = 0;
  // private viewerId: number = ++VdRegionsViewer.viewerCount;
  // private setDataInCount: number = 0;

  // TODO: .onTableAttached is not calling on dataFrame set, onPropertyChanged  also not calling
  public setData(regions: VdRegion[]) {
    const logPrefix = `${this.viewerToLog()}.setData()`;
    // const setDataInId = ++this.setDataInCount;
    _package.logger.debug(`${logPrefix}, in, ` +
      // `viewerId = ${this.viewerId}, setDataInId = ${setDataInId}, ` +
      `regions.length = ${regions.length}`
    );

    this.viewSyncer.sync(`${logPrefix}`, async () => { // setData
      // _package.logger.debug('Bio: VdRegionsViewer.setData(), in sync, ' +
      //   `viewerId = ${this.viewerId}, setDataInId = ${setDataInId}, ` +
      //   `regions.length = ${regions.length}`);
      if (!this.setDataInProgress) this.setDataInProgress = true; else return; // check setDataInProgress synced
      // _package.logger.debug('Bio: VdRegionsViewer.setData(), start, ' +
      //   `viewerId = ${this.viewerId}, setDataInId = ${setDataInId}, ` +
      //   `regions.length = ${regions.length}`);
      try {
        if (this.viewed) {
          // _package.logger.debug('Bio: VdRegionsViewer.setData(), destroyView, ' +
          //   `viewerId = ${this.viewerId}, setDataInId = ${setDataInId}, ` +
          //   `regions.length = ${regions.length}`);
          await this.destroyView('setData');
          this.viewed = false;
        }

        // -- Data --
        this.regions = regions;

        if (!this.viewed) {
          // _package.logger.debug('Bio: VdRegionsViewer.setData(), buildView, ' +
          //   `viewerId = ${this.viewerId}, setDataInId = ${setDataInId}, ` +
          //   `regions.length = ${regions.length}`);
          await this.buildView('setData');
          this.viewed = true;
        }
      } finally {
        // _package.logger.debug('Bio: VdRegionsViewer.setData(), finally, ' +
        //   `viewerId = ${this.viewerId}, setDataInId = ${setDataInId}, ` +
        //   `regions.length = ${regions.length}`);
        this.setDataInProgress = false;
      }
    });
  }

  // -- View --

  private viewSyncer: PromiseSyncer;
  private setDataInProgress: boolean = false;

  private host: HTMLElement | null = null;
  private filterSourceInput!: DG.InputBase<FilterSources | null>;
  private mainLayout: HTMLTableElement | null = null;
  private logos: { [chain: string]: WebLogoViewer }[] = [];

  private viewSubs: Unsubscribable[] = [];

  private async destroyView(purpose: string): Promise<void> {
    // TODO: Unsubscribe from and remove all view elements
    _package.logger.debug(`Bio: VdRegionsViewer.destroyView( mainLayout = ${!this.mainLayout ? 'none' : 'value'} ), ` +
      `purpose = '${purpose}', this.regions.length = ${this.regions.length}`);
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
    _package.logger.debug(`Bio: VdRegionsViewer.buildView() begin, ` +
      `purpose = '${purpose}', this.regions.length = ${this.regions.length}`);

    const regionsFiltered: VdRegion[] = this.regions.filter((r: VdRegion) => this.regionTypes.includes(r.type));
    const orderList: number[] = Array.from(new Set(regionsFiltered.map((r) => r.order))).sort();

    const logoPromiseList: Promise<[number, string, WebLogoViewer]>[] = [];
    for (let orderI = 0; orderI < orderList.length; orderI++) {
      for (const chain of this.chains) {
        const region: VdRegion | undefined = regionsFiltered
          .find((r) => r.order == orderList[orderI] && r.chain == chain);
        logoPromiseList.push((async (): Promise<[number, string, WebLogoViewer]> => {
          const wl: WebLogoViewer = await this.dataFrame.plot.fromType('WebLogo', {
            sequenceColumnName: region!.sequenceColumnName,
            startPositionName: region!.positionStartName,
            endPositionName: region!.positionEndName,
            fixWidth: true,
            skipEmptyPositions: this.skipEmptyPositions,
            positionWidth: this.positionWidth,
            positionHeight: this.positionHeight,
            filterSource: this.filterSource,
          }) as unknown as WebLogoViewer;
          wl.onSizeChanged.subscribe(() => { this.calcSize(); });
          return [orderI, chain, wl];
        })());
      }
    }
    const logoList: [number, string, WebLogoViewer][] = await Promise.all(logoPromiseList);
    // Fill in this.logos with created viewers
    this.logos = new Array(orderList.length);
    for (let orderI = 0; orderI < orderList.length; ++orderI)
      this.logos[orderI] = {};
    for (const [orderI, chain, wl] of logoList) {
      this.logos[orderI][chain] = wl;
      this.viewSubs.push(wl.onFreqsCalculated.subscribe(() => { this.calcSize(); }));
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
            },
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
              },
            });

            return resDiv;
          })];
        return elements;
      },
      ['', ...[...Array(orderList.length).keys()].map(
        (orderI: number) => regionsFiltered.find(
          (r: VdRegion) => r.order == orderList[orderI] && r.chain == this.chains[0],
        )!.name || 'Name')],
    );
    this.mainLayout.className = 'mlb-vd-regions-viewer-table2';
    // this.mainLayout.style.background = '#EEEEFF';
    // this.mainLayout.style.height = '100%';
    // this.mainLayout.style.border = '1px solid black';

    this.filterSourceInput = ui.input.choice<FilterSources>('Data source', {value: this.filterSource,
      items: Object.values(FilterSources), onValueChanged: this.filterSourceInputOnValueChanged.bind(this)});
    this.filterSourceInput.root.style.position = 'absolute';
    this.filterSourceInput.root.style.right = '9px';
    this.filterSourceInput.root.style.top = '-4px';
    //this.filterSourceInput.setTooltip('Check to filter sequences for selected VRs'); // TODO: GROK-13614
    //ui.tooltip.bind(this.filterSourceInput.input, 'Check to filter sequences for selected VRs');

    const _color: string = `#ffbb${Math.ceil(Math.random() * 255).toString(16)}`;
    this.host = ui.div([this.mainLayout, this.filterSourceInput!.root],
      {/*style: {backgroundColor: color}*/});
    this.root.appendChild(this.host);
    this.root.style.overflowX = 'auto';

    this.calcSize();
    this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));

    _package.logger.debug('Bio: VdRegionsViewer.buildView() end');
  }

  private calcSizeRequested: boolean = false;

  private calcSize() {
    _package.logger.debug(`Bio: VdRegionsViewer.calcSize(), start`);
    const calcSizeInt = (): void => {
      // Postponed calcSizeInt can result call after the viewer has been closed (on tests)
      if (!this.host) return;

      const logoHeight = (this.root.clientHeight - 54) / this.chains.length;
      let totalPos: number = 0;
      for (let orderI = 0; orderI < this.logos.length; orderI++) {
        for (const chain of this.chains) {
          const wl = this.logos[orderI][chain];
          wl.root.style.height = `${logoHeight}px`;
        }

        totalPos += Math.max(...this.chains.map((chain) => this.logos[orderI][chain].Length));
      }

      if (this.fitWidth) {
        if (this.logos.length > 0 && totalPos > 0) {
          const leftPad = 22/* Chain label */;
          const rightPad = 6 + 6 + 1;
          const logoMargin = 8 + 1;
          const fitPositionWidth =
            (this.root.clientWidth - leftPad - (this.logos.length - 1) * logoMargin - rightPad) / totalPos;

          for (let orderI = 0; orderI < this.logos.length; orderI++) {
            for (const chain of this.chains) {
              const wl = this.logos[orderI][chain];
              wl.setOptions({[wlPROPS.positionWidth]: (fitPositionWidth - wl.positionMarginValue)});
              wl.root.style.width = `${fitPositionWidth * wl.Length}px`;
            }
          }
        }
        this.host.style.setProperty('overflow', 'hidden', 'important');
      } else {
        for (let orderI = 0; orderI < this.logos.length; orderI++) {
          for (const chain of this.chains)
            this.logos[orderI][chain].setOptions({[wlPROPS.positionWidth]: this.positionWidth});
        }
        this.host.style.removeProperty('overflow');
      }

      if (this.positionWidth === 0)
        this.host!.style.setProperty('overflow-x', 'hidden', 'important');
      else
        this.host!.style.removeProperty('overflow-x');
    };

    if (!this.calcSizeRequested) {
      this.calcSizeRequested = true;
      window.setTimeout(() => {
        calcSizeInt();
        this.calcSizeRequested = false;
      }, 0 /* next event cycle */);
    }
  }

  // -- Handle events --

  private rootOnSizeChanged(_args: any): void {
    this.calcSize();
  }

  private rootOnMouseMove(_e: MouseEvent) {
    // ui.tooltip.show('text', e.x + 8, e.y + 8,);
    // console.log(`onMouseMoveRoot.( x: ${e.x}, y: ${e.y} )`);
  }

  private filterSourceInputOnValueChanged(): void {
    const logPrefix = `${this.viewerToLog()}.filterSourceInputOnValueChanged()`;
    const filterSourceValue = this.filterSourceInput.value;
    // Using promise to prevent 'Bad state: Cannot fire new event. Controller is already firing an event'
    this.viewSyncer.sync(`${logPrefix}`, async () => {
      if (this.filterSource !== filterSourceValue) {
        this.props.getProperty(PROPS.filterSource).set(this, filterSourceValue); // to update value in property panel

        for (let orderI = 0; orderI < this.logos.length; orderI++) {
          for (let chainI = 0; chainI < this.chains.length; chainI++) {
            const chain: string = this.chains[chainI];
            const wl: DG.JsViewer = this.logos[orderI][chain];
            wl.setOptions({[wlPROPS.filterSource]: this.filterSource});
          }
        }
      }
    });
  }

  // -- IRenderer --

  private _onRendered: Subject<void> = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    const logPrefix = `${this.viewerToLog()}.invalidate(${caller ? ` <- ${caller} ` : ''})`;
    // Put the event trigger in the tail of the synced calls queue.
    this.viewSyncer.sync(`${logPrefix}`, async () => {
      // update view / render
      // TODO: Request re-render
      this._onRendered.next();
    });
  }

  async awaitRendered(timeout: number | undefined = 5000): Promise<void> {
    await testEvent(this.onRendered, () => {}, () => {
      this.invalidate();
    }, timeout);

    // Rethrow stored syncer error (for test purposes)
    const viewErrors = this.viewSyncer.resetErrors();
    if (viewErrors.length > 0) throw viewErrors[0];
  }
}
