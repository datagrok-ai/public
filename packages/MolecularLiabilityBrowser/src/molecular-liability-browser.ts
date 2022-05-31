import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {DataLoader} from './utils/data-loader';
import {Subscription} from 'rxjs';
import {Aminoacids} from '@datagrok-libraries/bio/src/aminoacids';


export class MolecularLiabilityBrowser {
  dataLoader: DataLoader;

  // uriParser: HTMLAnchorElement;
  baseUri: string = '/apps/MolecularLiabilityBrowser/MolecularLiabilityBrowser/';
  urlParams: URLSearchParams = null;

  antigenDf: DG.DataFrame = null;
  antigenName: string = null;
  schemeName: string = 'chothia';
  schemeChoices: string[] = ['imgt', 'chothia', 'kabat'];

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
  antigenInput: DG.InputBase = null;
  antigenPopup: HTMLElement = null;
  schemeInput: DG.InputBase = null;


  constructor(dataLoader: DataLoader) {
    this.dataLoader = dataLoader;
  }

  async init(urlParams: URLSearchParams) {
    try {
      this.urlParams = urlParams;

      this.antigenDf = await this.dataLoader.listAntigens();
      this.antigenName = this.urlParams.get('antigen') || this.antigenDf.col('antigen').get(0);

      this.schemeName = this.urlParams.get('scheme') || this.schemeChoices[0];

      await this.loadMlbDf();
      await this.setView();

      this.vids = await this.dataLoader.getVids();
      this.vidsObsPTMs = await this.dataLoader.getObservedPtmVids();

      // this.mlbView = grok.shell.addTableView(this.mlbDf);
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

  setRibbonPanels(): void {
    this.mlbView.ribbonMenu.clear();

    this.setAntigenInput(this.antigenDf);
    this.setSchemeInput();

    this.mlbView.setRibbonPanels([
      [this.antigenInput.root],
      [this.schemeInput.root],
      [],
    ]);
  }

  static prepareDataMlbDf(df: DG.DataFrame, hChainDf: DG.DataFrame, lChainDf: DG.DataFrame): DG.DataFrame {
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
    });

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

  get mlbDf(): DG.DataFrame {
    return this._mlbDf;
  }

  async setMlbDf(value: DG.DataFrame): Promise<void> {
    await this.destroyView();
    this._mlbDf = value;
    await this.buildView();
  }

  /** Sets controls' layout. Called once from init(). */
  setView(): void {
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
              // this.vIdInput.value = gc.cell.value;
              // this.changeVid();
            })],
            {style: {position: 'absolute', top: 'calc(50% - 8px)', left: '5px'}});
        } else {
          gc.style.element = ui.divV([ui.label(gc.cell.value.toString())],
            {style: {position: 'absolute', top: 'calc(50% - 8px)', left: '5px'}});
        }
      }
    });
  }

  private async destroyView(): Promise<void> {
    // DG.TableView.dataFrame cannot be null
    // if (this.mlbView !== null)
    //    this.mlbView.dataFrame = null;

    this.idMapping = {};
    this.allIds = null;
    this.allVids = null;
  }

  private async buildView(): Promise<void> {
    this.allVids = this.mlbDf.col('v id')!;
    this.allIds = this.mlbDf.col('gdb id mappings');

    this.idMapping = {};
    for (let i = 0; i < this.allVids.length; i++)
      this.idMapping[this.allVids.get(i)] = this.allIds.get(i).replaceAll(' ', '').split(',');

    if (this.mlbView === null)
      this.mlbView = grok.shell.addTableView(this.mlbDf);
    else
      this.mlbView.dataFrame = this.mlbDf;

    this.updateView();
  }

  /** Restores column hiding, sets view path after dataFrame replacement. */
  updateView() {
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

  async onAntigenChanged(antigenName: string): Promise<void> {
    this.antigenName = antigenName;

    // Preventing firing the event
    if (this.antigenInput.value != this.antigenName)
      this.antigenInput.value = this.antigenName;

    this.urlParams.set('antigen', this.antigenName);

    window.setTimeout(async () => { await this.loadMlbDf(); }, 10);

    // const treeDf: DG.DataFrame = await this.dataLoader.getTreeByAntigen(this.antigenName);
    // this.treeBrowser.treeDf = treeDf;
  }

  onSchemeChanged(schemeName: string): void {
    this.schemeName = schemeName;

    // Preventing firing the event
    if (this.schemeInput.value != this.schemeName)
      this.schemeInput.value = this.schemeName;

    this.urlParams.set('scheme', this.schemeName);

    window.setTimeout(async () => { await this.loadMlbDf(); }, 10);
  }

  /** Loads MLB data and sets prepared DataFrame */
  async loadMlbDf(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Loading data...');
    try {
      // let mlbDf: DG.DataFrame;
      // let hChainDf: DG.DataFrame;
      // let lChainDf: DG.DataFrame;

      const t1 = Date.now();
      const [mlbDf, hChainDf, lChainDf]: [DG.DataFrame, DG.DataFrame, DG.DataFrame] = (await Promise.all([
        this.dataLoader.getMlbByAntigen(this.antigenName),
        this.dataLoader.getAnarci(this.schemeName, 'heavy', this.antigenName),
        this.dataLoader.getAnarci(this.schemeName, 'light', this.antigenName),
      ]));
      const t2 = Date.now();
      console.debug(`MolecularLiabilityBrowser.loadMlbDf() load duration ${((t2 - t1) / 1000).toString()} s`);

      await this.setMlbDf(MolecularLiabilityBrowser.prepareDataMlbDf(mlbDf, hChainDf, lChainDf));
      const t3 = Date.now();
      console.debug(`MolecularLiabilityBrowser.loadMlbDf() prepare ${((t3 - t2) / 1000).toString()} s`);
    } finally {
      pi.close();
    }
  }
}
