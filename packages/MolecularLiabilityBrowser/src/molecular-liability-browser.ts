import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {DataLoader, FilterPropertiesType, JsonType, NumsType, ObsPtmType, PdbType} from './utils/data-loader';
import {Subscription} from 'rxjs';
import {Aminoacids} from '@datagrok-libraries/bio/src/aminoacids';
import {VdRegionsViewer} from './viewers/vd-regions-viewer';
import {TreeBrowser} from './mlb-tree';
import {TreeAnalyzer} from './utils/tree-stats';
import {MiscMethods} from './viewers/misc';
import {TwinPviewer} from './viewers/twin-p-viewer';


export class MolecularLiabilityBrowser {
  dataLoader: DataLoader;

  // uriParser: HTMLAnchorElement;
  baseUri: string = '/apps/MolecularLiabilityBrowser/MolecularLiabilityBrowser/';
  urlParams: URLSearchParams = null;

  antigenDf: DG.DataFrame = null;
  pf: FilterPropertiesType = null;
  antigenName: string = null;
  schemeName: string = 'chothia';
  schemeChoices: string[] = ['imgt', 'chothia', 'kabat'];

  /** Current antibody selected molecule */
  vid: string = null;

  refDf: DG.DataFrame = null; // from file ptm_in_cdr.d42
  hChainDf: DG.DataFrame = null;
  lChainDf: DG.DataFrame = null;

  vids: string[] = null;
  vidsObsPTMs: string[] = null;

  allVids: DG.Column = null;
  allIds: DG.Column = null;
  idMapping: { [key: string]: string []; };

  subs: Subscription[] = [];

  mlbView: DG.TableView = null;
  regionsViewer: VdRegionsViewer = null;
  treeBrowser: TreeBrowser = null;
  twinPviewer: TwinPviewer;

  antigenInput: DG.InputBase = null;
  antigenPopup: HTMLElement = null;
  schemeInput: DG.InputBase = null;
  hideShowIcon: HTMLElement = null;

  constructor(dataLoader: DataLoader) {
    this.dataLoader = dataLoader;
  }

  async init(urlParams: URLSearchParams) {
    try {
      this.urlParams = urlParams;

      this.pf = this.dataLoader.filterProperties;

      this.antigenDf = await this.dataLoader.listAntigens();
      if (!this.antigenDf) {
        const msg: string = `Get antigen list failed.`;
        grok.shell.error(msg);
        throw new Error(msg);
      }

      this.antigenName = this.urlParams.get('antigen') || this.antigenDf.col('antigen').get(0);
      this.schemeName = this.schemeChoices.includes(this.urlParams.get('scheme')) ? this.urlParams.get('scheme') :
        this.schemeChoices[0];

      [this.vids, this.vidsObsPTMs] = (await Promise.all([
        this.dataLoader.getVids(),
        this.dataLoader.getObservedPtmVids()
      ]));

      await this.loadData();
      await this.setView();
      await this.setViewFilters();
      await this.setViewTwinPViewer();
    } finally {
      grok.events.onViewRemoved.subscribe((v) => {
        if (v.type === DG.VIEW_TYPE.TABLE_VIEW && (v as DG.TableView).dataFrame.id === this.mlbDf.id)
          this.subs.forEach((s) => s.unsubscribe());
      });
    }
  }

  setAntigenInput(agDf: DG.DataFrame): DG.InputBase {
    this.antigenInput = ui.stringInput('AG', this.antigenName, null, {tabindex: '1'} as unknown);

    const agCol: DG.Column = agDf.col('antigen');
    const gsCol: DG.Column = agDf.col('antigen_gene_symbol');
    const agDfGrid: DG.Grid = agDf.plot.grid({
      allowEdit: false,
      allowRowSelection: false,
      allowRowResizing: false,
      allowRowReordering: false,
      allowColReordering: false,
      allowBlockSelection: false,
      showRowHeader: false,
    });
    const idGCol: DG.GridColumn = agDfGrid.col('id');
    const agGCol: DG.GridColumn = agDfGrid.col('antigen');
    const ncbiGCol: DG.GridColumn = agDfGrid.col('antigen_ncbi_id');
    const gsGCol: DG.GridColumn = agDfGrid.col('antigen_gene_symbol');
    agGCol.name = 'Antigen';
    agGCol.width = 90;
    gsGCol.name = 'Gene symbol';
    gsGCol.width = 120;
    idGCol.visible = false;
    ncbiGCol.visible = false;
    agDfGrid.root.style.setProperty('width', '220px');
    this.antigenPopup = ui.div([agDfGrid.root], {tabindex: '2'} as unknown);

    agDf.onCurrentRowChanged.subscribe(() => {
      this.antigenPopup.hidden = true;

      const antigenName: string = agCol.get(agDf.currentRow.idx);
      window.setTimeout(this.onAntigenChanged.bind(this), 10, antigenName);
    });

    this.antigenInput.root.addEventListener('input', (event: Event) => {
      /* Here we should filter dataframe with antigens */
      agDf.filter.init((iRow: number) => {
        return agCol.get(iRow).includes(this.antigenInput.value) || gsCol.get(iRow).includes(this.antigenInput.value);
      });
    });
    // this.agInput.root.addEventListener('focusout', (event: FocusEvent) => {
    //   let k = 11;
    //   agDfGrid.root.focus();
    // });
    // this.agInput.root.addEventListener('keydown', (event: KeyboardEvent) => {
    //   if (event.code === '\t') {
    //     let k = 11;
    //   }
    // });

    this.antigenInput.root.addEventListener('mousedown', (event: MouseEvent) => {
      this.antigenPopup.hidden = false;
      ui.showPopup(this.antigenPopup, this.antigenInput.root);
    });

    return this.antigenInput;
  }

  setSchemeInput(): DG.InputBase {
    this.schemeInput = ui.choiceInput('Scheme', this.schemeName, this.schemeChoices, () => {
      // this.regionsViewer.setScheme(this.schemeInput.value);
      this.onSchemeChanged(this.schemeInput.value);
    });

    return this.schemeInput;
  }

  setHideShowIcon(): void {
    this.hideShowIcon = ui.tooltip.bind(ui.iconFA('eye', () => {
      if (this.hideShowIcon.classList.value.includes('fa-eye-slash'))
        grok.events.fireCustomEvent('showAllDock', null);
      else
        grok.events.fireCustomEvent('closeAllDock', null);
    }), 'show structure');
    this.hideShowIcon.classList.value = 'grok-icon fal fa-eye';

    grok.events.onCustomEvent('closeAllDock').subscribe((v) => {
      this.twinPviewer.close(this.mlbView);
      this.hideShowIcon.classList.value = 'grok-icon fal fa-eye-slash';
    });

    grok.events.onCustomEvent('showAllDock').subscribe((v) => {
      this.twinPviewer.open(this.mlbView);
      this.hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    });
  }


  setRibbonPanels(): void {
    this.mlbView.ribbonMenu.clear();

    this.setAntigenInput(this.antigenDf);
    this.setSchemeInput();
    this.setHideShowIcon();

    this.mlbView.setRibbonPanels([
      [this.antigenInput.root],
      [this.schemeInput.root],
      [],
      [this.hideShowIcon],
    ]);
  }

  static prepareDataMlbDf(df: DG.DataFrame, numberingScheme: string,
    hChainDf: DG.DataFrame, lChainDf: DG.DataFrame, pf: FilterPropertiesType): DG.DataFrame {
    for (const column of df.columns)
      column.name = column.name.replaceAll('_', ' ');

    // TODO: Obsolete filtering by antigen
    const mlbColumnsToMultiValue = ['antigen list', 'antigen ncbi id', 'antigen gene symbol'];
    for (const column of df.columns) {
      if (mlbColumnsToMultiValue.includes(column.name))
        column.setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, ',');
    }

    // TODO: Load chains sequences
    console.debug(`hChainDf: ${hChainDf.rowCount} rows, lChainDf: ${lChainDf.rowCount} rows`);
    [{name: 'Heavy', df: hChainDf}, {name: 'Light', df: lChainDf}].forEach((chain) => {
      const seqCol = this.mergeSequenceColumns(chain.df, chain.name);
      grok.data.joinTables(
        df, chain.df,
        ['v id'], ['Id'],
        (df.columns as DG.ColumnList).names(), [seqCol.name],
        DG.JOIN_TYPE.LEFT, true);

      // crutch, because grok.data.joinTables() loses right table columns tags
      df.col(seqCol.name).setTag('positionNames', chain.df.col(seqCol.name).getTag('positionNames'));
      df.col(seqCol.name).setTag('numberingScheme', numberingScheme);
    });

    // Adding columns after tableView breaks display (no columns show)
    const clonesCol = df.columns.addNewString('clones');

    //bands on plots for properties
    for (let i = 0; i < pf.names.length; i++) {
      df.col(pf.names[i])!.setTag(
        '.default-filter', '{ "min": ' + pf.yellowLeft[i] + ', "max": ' + pf.yellowRight[i] + ' }');
      // this.mlbTable.col(pf.names[i])!.setTag('.default-filter', JSON.stringify(
      //   {
      //     min: pf.yellowLeft[i],
      //     max: pf.yellowRight[i]
      //   }));
      df.col(pf.names[i])!.setTag('.charts', JSON.stringify([
        {
          title: 'BandYellowLeft', type: 'band',
          rule: `${pf.redLeft[i]}-${pf.yellowLeft[i]}`,
          color: '#FFD700', opacity: 15
        },
        {
          title: 'BandYellowRight', type: 'band',
          rule: `${pf.yellowRight[i]}-${pf.redRight[i]}`,
          color: '#FFD700', opacity: 15
        },
        {title: 'BandRedLeft', type: 'band', rule: `< ${pf.redLeft[i]}`, color: '#DC143C', opacity: 15},
        {title: 'BandRedRight', type: 'band', rule: `> ${pf.redRight[i]}`, color: '#DC143C', opacity: 15},
        {
          'title': 'TAP metrics', 'type': 'spline',
          'color': '#7570B3', 'width': 1, 'normalize-y': true, 'visible': true,
          'x': pf.plotsX[i], 'y': pf.plotsY[i],
        }
      ]));
    }

    return df;
  }

  static prepareDataTreeDf(df: DG.DataFrame, treeColumnName: string = 'TREE') {
    function _modifyTreeNodeIds(nwk: string): string {
      if (TreeAnalyzer.newickRegEx.test(nwk.trim()))
        return nwk.replaceAll(/([^|,:()]+)\|([^|,:()]+)\|([^|,:()]+)\|([^|,:()]+)/g, '$3');

      return nwk;
    }

    const treeCol = df.col(treeColumnName);
    const trees = treeCol.toList().map((v) => _modifyTreeNodeIds(v));
    (df.columns as DG.ColumnList).replace(treeColumnName, DG.Column.fromStrings(treeColumnName, trees));

    return df;
  }

  /** Builds multiple alignment sequences from monomers in positions
   *  as virtual (calculated) column in source DataFrame.
   * @param {DG.DataFrame} df  DataFrame with ANARCI results
   * @param {string} chain  Name of chain (used for result column name)
   * @param {number} startingColumnIndex  The first column with positions
   * @return {DG.Column}
   */
  static mergeSequenceColumns(df: DG.DataFrame, chain: string, startingColumnIndex?: number) {
    const positionRegExp = /^\d+[A-Z]*$/g;
    const columns: DG.ColumnList = df.columns;

    // All positions schemes contains position '1'
    const positionNames = startingColumnIndex !== void 0 ? columns.names().slice(startingColumnIndex) :
      columns.names().slice(columns.names().indexOf('1'));

    const positionColumns = positionNames.filter((v: string) => v.match(positionRegExp) !== null);
    const seqCol = columns.addNewVirtual(
      `${chain} chain sequence`,
      (i: number) => positionColumns.map((v) => df.get(v, i)).join('')
    );
    seqCol.semType = Aminoacids.SemTypeMultipleAlignment;

    const positionNamesTxt = positionNames.join(', '); /* Spaces are for word wrap */
    seqCol.setTag('positionNames', positionNamesTxt);

    return seqCol;
  }

  private _mlbDf: DG.DataFrame = DG.DataFrame.fromObjects([]);
  private _treeDf: DG.DataFrame = DG.DataFrame.fromObjects([]);

  get mlbDf(): DG.DataFrame { return this._mlbDf; }

  get treeDf(): DG.DataFrame { return this._treeDf; }

  async setData(mlbDf: DG.DataFrame, treeDf: DG.DataFrame): Promise<void> {
    await this.destroyView();
    // .catch((ex) => {
    //   console.error(`MolecularLiabilityBrowser.setData() > destroyView() error:\n${ex.toString()}`);
    // });
    this._mlbDf = mlbDf;
    this._treeDf = treeDf;
    await this.buildView()
      .catch((ex) => {
        console.error(`MolecularLiabilityBrowser.setData() > buildView() error:\n${ex ? ex.toString() : 'none'}`);
      });
  }

  /** Sets controls' layout. Called once from init(). */
  async setView(): Promise<void> {
    console.debug('MolecularLiabilityBrowser.setView()');
    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;
    this.setRibbonPanels();

    this.mlbView.name = 'Molecular Liability Browser';
    for (const column of this.mlbDf.columns) {
      const gridColumn: DG.GridColumn = this.mlbView.grid.columns.byName(column.name);
      gridColumn.name = column.name.replaceAll('_', ' ');
    }

    this.mlbView.grid.columns.byName('v id')!.width = 120;
    this.mlbView.grid.columns.byName('v id')!.cellType = 'html';


    //table visual polishing
    this.mlbView.grid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => {
      if (args.cell.isColHeader) {
        if (args.cell.gridColumn.visible) {
          const textSize = args.g.measureText(args.cell.gridColumn.name);
          args.g.fillText(args.cell.gridColumn.name, args.bounds.x +
            (args.bounds.width - textSize.width) / 2, args.bounds.y +
            (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent));
          args.g.fillStyle = '#4b4b4a';
        }
        args.preventDefault(); // this is required to prevent drawing headers of hidden columns
      }
    });

    this.mlbView.grid.onCellPrepare((gc) => {
      if (gc.isTableCell && gc.gridColumn.name === 'v id') {
        if (this.vids.includes(gc.cell.value.toString())) {
          gc.style.element = ui.divV(
            [ui.link(gc.cell.value.toString(), () => {
              this.vid = gc.cell.value;
              this.changeVid();
            })],
            {style: {position: 'absolute', top: 'calc(50% - 8px)', left: '5px'}});
        } else {
          gc.style.element = ui.divV([ui.label(gc.cell.value.toString())],
            {style: {position: 'absolute', top: 'calc(50% - 8px)', left: '5px'}});
        }
      }
    });
  }

  setViewFilters(): void {
    const filterList: { type: string, column?: string, label?: string }[] = [];

    for (const pfName of this.pf.names) {
      // this.mlbTable.columns.byName(pfName).width = 150;
      this.mlbView.grid.col(pfName).width = 150;
      filterList.push({type: 'histogram', column: pfName});
    }
    filterList.push({type: 'MolecularLiabilityBrowser:ptmFilter'});

    const filterView = this.mlbView.filters({filters: filterList});
    this.mlbView.dockManager.dock(filterView, DG.DOCK_TYPE.LEFT, null, 'Filters', 0.25);

    grok.events.onTooltipShown.subscribe((args) => {
      if (args.args.context instanceof DG.Column) {
        const colName: string = args.args.context.name;
        switch (colName) {
        case 'cdr length':
          args.args.element.innerHTML = this.pf.tooltips[0];
          break;

        case 'surface cdr hydrophobicity':
          args.args.element.innerHTML = this.pf.tooltips[1];
          break;

        case 'positive cdr charge':
          args.args.element.innerHTML = this.pf.tooltips[2];
          break;

        case 'negative cdr charge':
          args.args.element.innerHTML = this.pf.tooltips[3];
          break;

        case 'sfvcsp':
          args.args.element.innerHTML = this.pf.tooltips[4];
          break;

        default:
          console.error('MolecularLiabilityBrowser.setViewFilters() onToolTipShown() ' +
            `unexpected context name - ${colName}`);
          break;
        }
      }
    });
  }

  async setViewTwinPViewer(): Promise<void> {
    const [jsonStr, pdbStr, jsonNums, jsonStrObsPtm]: [JsonType, string, NumsType, ObsPtmType] = (
      await Promise.all([
        this.dataLoader.loadExample(this.vid),
        this.dataLoader.loadPdb(this.vid),
        this.dataLoader.realNums,
        ((): Promise<ObsPtmType> => {
          let res: Promise<ObsPtmType> = null;
          if (this.vidsObsPTMs.includes(this.vid))
            res = this.dataLoader.loadObsPtm(this.vid);
          return res;
        })()
      ]));

    this.twinPviewer = new TwinPviewer(this.dataLoader);
    this.twinPviewer.init(jsonStr, pdbStr, jsonStrObsPtm);

    await this.twinPviewer.show(this.mlbView);
    await this.twinPviewer.open(this.mlbView);
  }

  async destroyView(): Promise<void> {
    console.debug('MolecularLiabilityBrowser.destroyView()');
    // DG.TableView.dataFrame cannot be null
    // if (this.mlbView !== null)
    //    this.mlbView.dataFrame = null;

    this.idMapping = {};
    this.allIds = null;
    this.allVids = null;
  }

  async buildView(): Promise<void> {
    console.debug('MolecularLiabilityBrowser.buildView()');

    this.allVids = this.mlbDf.col('v id')!;
    this.allIds = this.mlbDf.col('gdb id mappings');

    this.idMapping = {};
    for (let i = 0; i < this.allVids.length; i++)
      this.idMapping[this.allVids.get(i)] = this.allIds.get(i).replaceAll(' ', '').split(',');

    if (this.mlbView === null) {
      // this.mlbView = grok.shell.addView();
      this.mlbView = grok.shell.addTableView(this.mlbDf);
    } else {
      this.mlbView.dataFrame = this.mlbDf;
    }

    if (this.treeBrowser === null) {
      this.treeBrowser = (await this.treeDf.plot.fromType('MlbTree', {})) as unknown as TreeBrowser;
      await this.treeBrowser.setData(this.treeDf, this.mlbDf);
      this.mlbView.dockManager.dock(this.treeBrowser, DG.DOCK_TYPE.RIGHT, null, 'Clone', 0.5);
    } else {
      await this.treeBrowser.setData(this.treeDf, this.mlbDf);
    }

    if (this.regionsViewer === null) {
      this.regionsViewer = (await this.mlbDf.plot.fromType('VdRegions', {
        numberingScheme: this.schemeName,
      })) as unknown as VdRegionsViewer;
      this.mlbView.dockManager.dock(this.regionsViewer, DG.DOCK_TYPE.DOWN, null, 'Regions', 0.3);
    } else {
      this.regionsViewer.numberingScheme = this.schemeName;
      // this.regionsViewer.setOptions({numberingScheme: this.schemeName});
      await this.regionsViewer.setDf(this.mlbDf);
    }

    this.updateView();

    this.mlbDf.onCurrentRowChanged.subscribe(this.onMLBGridCurrentRowChanged.bind(this));
  }

  /** Restores column hiding, sets view path after dataFrame replacement. */
  updateView() {
    console.debug('MolecularLiabilityBrowser.updateView()');
    const urlParamsTxt = Array.from(this.urlParams.entries())
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');

    this.mlbView.path = `${this.baseUri}?${urlParamsTxt}`;

    // Leonid instructed to hide the columns
    const mlbColumnsToHide = [].concat(...[
      ['cdr length', 'surface cdr hydrophobicity', 'positive cdr charge', 'negative cdr charge', 'SFvCSP'],
      ['antigen list', 'antigen ncbi id', 'antigen gene symbol'],
      ['Heavy chain sequence', 'Light chain sequence']
    ]);
    for (let colI = 0; colI < this.mlbView.grid.columns.length; colI++) {
      const gridColumn: DG.GridColumn = this.mlbView.grid.columns.byIndex(colI);
      if (gridColumn.column !== null && mlbColumnsToHide.includes(gridColumn.column.name))
        gridColumn.visible = false;
    }
  }

  onMLBGridCurrentRowChanged(args: any) {
    this.vid = this.mlbDf.currentRow['v id'];
    window.setTimeout(async () => { await this.changeVid(true); }, 10);
  }

  onAntigenChanged(antigenName: string): void {
    this.antigenName = antigenName;

    // Preventing firing the event
    if (this.antigenInput.value != this.antigenName)
      this.antigenInput.value = this.antigenName;

    this.urlParams.set('antigen', this.antigenName);

    window.setTimeout(async () => { await this.loadData(); }, 10);

    // const treeDf: DG.DataFrame = await this.dataLoader.getTreeByAntigen(this.antigenName);
    // this.treeBrowser.treeDf = treeDf;
  }

  onSchemeChanged(schemeName: string): void {
    this.schemeName = schemeName;

    // Preventing firing the event
    if (this.schemeInput.value != this.schemeName)
      this.schemeInput.value = this.schemeName;

    this.urlParams.set('scheme', this.schemeName);

    window.setTimeout(async () => { await this.loadData(); }, 10);
  }

  onVidChanged(vid: string): void {
    this.vid = vid;

    this.urlParams.set('vid', this.vid);

    window.setTimeout(async () => { await this.changeVid(); }, 10);
  }

  /** Loads MLB data and sets prepared DataFrame */
  async loadData(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Loading data...');
    try {
      // let mlbDf: DG.DataFrame;
      // let hChainDf: DG.DataFrame;
      // let lChainDf: DG.DataFrame;

      const t1 = Date.now();
      const [mlbDf, hChainDf, lChainDf, treeDf]: [DG.DataFrame, DG.DataFrame, DG.DataFrame, DG.DataFrame] =
        (await Promise.all([
          this.dataLoader.getMlbByAntigen(this.antigenName),
          this.dataLoader.getAnarci(this.schemeName, 'heavy', this.antigenName),
          this.dataLoader.getAnarci(this.schemeName, 'light', this.antigenName),
          this.dataLoader.getTreeByAntigen(this.antigenName),
        ]));
      const t2 = Date.now();
      console.debug(`MolecularLiabilityBrowser.loadMlbDf() load duration ${((t2 - t1) / 1000).toString()} s`);

      this.setData(
        MolecularLiabilityBrowser.prepareDataMlbDf(mlbDf, this.schemeName, hChainDf, lChainDf, this.pf),
        MolecularLiabilityBrowser.prepareDataTreeDf(treeDf));
      ;
      const t3 = Date.now();
      console.debug(`MolecularLiabilityBrowser.loadMlbDf() prepare ${((t3 - t2) / 1000).toString()} s`);
    } finally {
      pi.close();
    }
  }

  // --

  async changeVid(silent: boolean = false): Promise<void> {
    if (!this.vids.includes(this.vid)) {
      Object.keys(this.idMapping).every((vid) => {
        if (this.idMapping[vid].includes(this.vid)) {
          this.vid = vid;
          return false;
        }
        return true;
      });

      if (!this.vids.includes(this.vid)) {
        grok.shell.warning('No PDB data data for associated v id');
        return;
      }
    }

    // this.mlbView.path = `/Table/${this.vIdInput.value}`;
    //hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

    const [jsonStr, pdbStr, jsonNums, jsonStrObsPtm]: [JsonType, string, NumsType, ObsPtmType] = (
      await Promise.all([
        this.dataLoader.loadExample(this.vid),
        this.dataLoader.loadPdb(this.vid),
        this.dataLoader.realNums,
        ((): Promise<ObsPtmType> => {
          let res: Promise<ObsPtmType> = null;
          if (this.vidsObsPTMs.includes(this.vid))
            res = this.dataLoader.loadObsPtm(this.vid);
          return res;
        })()
      ]));

    const hNumberingStr = jsonNums.heavy_numbering;
    const lNumberingStr = jsonNums.light_numbering;

    const hNumbering = [];
    const lNumbering = [];

    for (let i = 0; i < hNumberingStr.length; i++)
      hNumbering.push(parseInt(hNumberingStr[i].replaceAll(' ', '')));

    for (let i = 0; i < lNumberingStr.length; i++)
      lNumbering.push(parseInt(lNumberingStr[i].replaceAll(' ', '')));

    jsonStr['map_H'] = hNumbering;
    jsonStr['map_L'] = lNumbering;

    await this.twinPviewer.reset(jsonStr, pdbStr, jsonStrObsPtm);
    await this.twinPviewer.show(this.mlbView);
    await this.twinPviewer.open(this.mlbView);

    pi.close();
  };
}
