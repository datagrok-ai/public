import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {TwinPviewer} from './viewers/twin-p-viewer';
import {TreeBrowserOld} from './tree-old';
import {Subscription} from 'rxjs';
import {
  DataLoader,
  JsonType,
  PdbType,
  NumsType,
  ObsPtmType,
} from './utils/data-loader';
import {DataLoaderFiles} from './utils/data-loader-files';
import {Aminoacids} from '@datagrok-libraries/bio/src/aminoacids';
import {IVdRegionsViewer} from '@datagrok-libraries/bio/src/viewers/vd-regions-viewer';

// import {WebLogo} from '@datagrok-libraries/bio';


export class MolecularLiabilityBrowserOld {
  dataLoader: DataLoader;

  constructor(dataLoader: DataLoader) {
    this.dataLoader = dataLoader;
  }

  mlbTable: DG.DataFrame;
  vids: string[];
  vidsObsPTMs: string[];
  allVids: DG.Column;
  allIds: DG.Column;
  idMapping: { [key: string]: string[]; };

  mlbView: DG.TableView;
  vIdInput: DG.InputBase;
  filterIcon: HTMLElement;
  filterIcon2: HTMLElement;
  selectionClearIcon: HTMLElement;
  queryIdIcon: HTMLElement;
  queryAntigenIcon: HTMLElement;
  hideShowIcon: HTMLElement;
  hideShowIcon2: HTMLElement;
  treesIcon: HTMLElement;

  regionsViewer: IVdRegionsViewer;
  twinPviewer: TwinPviewer;
  //compostionPviewer: CompostionPviewer;

  subs: Subscription[] = [];

  private changeVid = async (silent: boolean = false): Promise<void> => {
    if (!this.vids.includes(this.vIdInput.value)) {
      Object.keys(this.idMapping).every((vid) => {
        if (this.idMapping[vid].includes(this.vIdInput.value)) {
          this.vIdInput.value = vid;
          return false;
        }
        return true;
      });

      if (!this.vids.includes(this.vIdInput.value)) {
        grok.shell.warning('No PDB data data for associated v id');
        return;
      }
    }

    // #region Commented to replace from RepertoireBrowser
    // this.mlbView.path = `/Table/${this.vIdInput.value}`;
    //hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    // const pi = DG.TaskBarProgressIndicator.create('Creating NGL view');
    // #endregion

    // #region Code from RepertoireBrowser
    this.mlbView.path = `/Table/${this.vIdInput.value}`;
    //hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

    ////////////////////////////////////////////////////
    // const jsonStr: JsonType = await this.dataLoader.loadJson(this.vIdInput.value);
    //
    // const pdbStr: PdbType = await this.dataLoader.loadPdb(this.vIdInput.value);
    //
    // const jsonNums: NumsType = await this.dataLoader.loadRealNums(this.vIdInput.value);
    //
    // let jsonStrObsPtm: ObsPtmType = null;
    // if (this.vidsObsPTMs.includes(this.vIdInput.value))
    //   jsonStrObsPtm = await this.dataLoader.loadObsPtm(this.vIdInput.value);
    const [jsonStr, pdbStr, jsonNums, jsonStrObsPtm]:
      [JsonType, string, NumsType, ObsPtmType] = await this.dataLoader.load3D(this.vIdInput.value);

    ////////////////////////////////////////////////////
    // #region

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

    if (!this.twinPviewer) {
      this.twinPviewer = new TwinPviewer(this.dataLoader);
      this.twinPviewer.init(jsonStr, pdbStr, jsonStrObsPtm);
    } else {
      this.twinPviewer.reset(jsonStr, pdbStr, jsonStrObsPtm);
    }

    await this.twinPviewer.show(this.mlbView);
    await this.twinPviewer.open(this.mlbView);

    pi.close();
  };

  setPropertiesFilters(): void {
    //external data load
    const pf = this.dataLoader.filterProperties;

    //this.mlbView.grid.columns.byName('CDR Clothia').visible = false;

    //bands on plots for properties
    for (let i = 0; i < pf.names.length; i++) {
      this.mlbTable.col(pf.names[i])!.setTag(
        '.default-filter', '{ "min": ' + pf.yellowLeft[i] + ', "max": ' + pf.yellowRight[i] + ' }');
      // this.mlbTable.col(pf.names[i])!.setTag('.default-filter', JSON.stringify(
      //   {
      //     min: pf.yellowLeft[i],
      //     max: pf.yellowRight[i]
      //   }));
      this.mlbTable.col(pf.names[i])!.setTag('.charts', JSON.stringify([
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

      // this.mlbTable.col(pf.names[i])!.setTag('.charts', '[' +
      //   '{"type": "band", "title":"BandYellowLeft", "rule" : "' +
      //   pf.redLeft[i] + '-' + pf.yellowLeft[i] + '", "color": "#FFD700", "opacity": 15},' +
      //   '{"type": "band", "title":"BandYellowRight", "rule" : "' +
      //   pf.yellowRight[i] + '-' + pf.redRight[i] + '", "color": "#FFD700", "opacity": 15}, ' +
      //   '{"type": "band", "title":"BandRedLeft", "rule" : "< ' + pf.redLeft[i] +
      //   '", "color": "#DC143C", "opacity": 15}, ' +
      //   '{"type": "band", "title":"BandRedRight", "rule" : "> ' + pf.redRight[i] +
      //   '", "color": "#DC143C", "opacity": 15},' +
      //   '{ "type": "spline", "title": "TAP metrics", "y" : [' +
      //   pf.plotsY[i].toString() + '], "color": "#7570B3", "width": 1, "x" : [' +
      //   pf.plotsX[i].toString() + '], "normalize-y": true, "visible": true}' +
      //   ']');
    }

    const filterList: { type: string, column?: string, label?: string }[] = [];

    for (const pfName of pf.names) {
      // this.mlbTable.columns.byName(pfName).width = 150;
      this.mlbView.grid.col(pfName)!.width = 150;
      filterList.push({type: 'histogram', column: pfName});
    }
    filterList.push({type: 'MolecularLiabilityBrowser:ptmFilter'});
    filterList.push({type: DG.FILTER_TYPE.MULTI_VALUE, column: 'antigen list', label: 'antigen id'});
    filterList.push({type: DG.FILTER_TYPE.MULTI_VALUE, column: 'antigen gene symbol', label: 'antigen gene symbol'});

    const filterView = this.mlbView.filters({filters: filterList});

    grok.events.onTooltipShown.subscribe((args) => {
      if (args.args.context instanceof DG.Column) {
        switch (args.args.context.name) {
        case 'cdr_length':
          args.args.element.innerHTML = pf.tooltips[0];
          break;

        case 'surface_cdr_hydrophobicity':
          args.args.element.innerHTML = pf.tooltips[1];
          break;

        case 'positive_cdr_charge':
          args.args.element.innerHTML = pf.tooltips[2];
          break;

        case 'negative_cdr_charge':
          args.args.element.innerHTML = pf.tooltips[3];
          break;

        case 'sfvcsp':
          args.args.element.innerHTML = pf.tooltips[4];
          break;

        default:
          break;
        }
      }
    });
  }

  async setView(): Promise<void> {
    this.mlbView = grok.shell.addTableView(this.mlbTable);

    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;
    this.mlbView.ribbonMenu.clear();

    this.mlbView.name = 'Molecular Liability Browser';
    // for (const column of this.mlbTable.columns)
    //   column.name = column.name.replaceAll('_', ' ');
    for (const column of this.mlbTable.columns) {
      const gridColumn: DG.GridColumn = this.mlbView.grid.columns.byName(column.name)!;
      gridColumn.name = column.name.replaceAll('_', ' ');
    }

    this.mlbView.grid.columns.byName('v id')!.width = 120;
    this.mlbView.grid.columns.byName('v id')!.cellType = 'html';

    // Leonid instructed to hide the columns
    const mlbColumnsToHide: string[] = [
      ...['cdr length', 'surface cdr hydrophobicity', 'positive cdr charge', 'negative cdr charge', 'sfvcsp'],
      ...['antigen list', 'antigen ncbi id', 'antigen gene symbol'],
      ...['Heavy chain sequence', 'Light chain sequence']
    ];
    for (let colI = 0; colI < this.mlbView.grid.columns.length; colI++) {
      const gridColumn: DG.GridColumn = this.mlbView.grid.columns.byIndex(colI)!;
      if (gridColumn.column !== null && mlbColumnsToHide.includes(gridColumn.column.name))
        gridColumn.visible = false;
    }

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
              this.vIdInput.value = gc.cell.value;
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

  setVidInput(): void {
    this.vIdInput = ui.stringInput('VID', '');
    this.vIdInput.value = this.vids[0];

    this.vIdInput.root.addEventListener('keyup', (event) => {
      if (event.key === 'Enter')
        this.changeVid();
    });
  }

  // setFilterIcon(): void {
  //   this.filterIcon = ui.tooltip.bind(ui.iconFA('filter', () => {
  //     const a = this.mlbTable;
  //     a.rows.select((row) => this.vids.includes(row['v id'].toString()));
  //     grok.functions.call('CmdSelectionToFilter');
  //   }), 'filter data for 3D view');
  // }

  setFilterIcon2(): void {
    this.filterIcon2 = ui.tooltip.bind(ui.iconFA('filter', () => {
      const a = this.mlbTable;
      a.rows.select((row) => this.vidsObsPTMs.includes(row['v id'].toString()));
      grok.functions.call('CmdSelectionToFilter');
    }), 'Filter data with observed PTMs');
  }

  setSelectionClearIcon(): void {
    this.selectionClearIcon = ui.tooltip.bind(ui.iconFA('broom', () => {
      this.mlbTable.rows.select((i) => false);
    }), 'Selection clear');
  }

  setQueryIdIcon(): void {
    this.queryIdIcon = ui.tooltip.bind(ui.iconFA('layer-group', () => {
      //get all possible IDs
      const allIds = this.mlbTable.col('v id')!.toList().concat(
        this.mlbTable.col('gdb id mappings')!.toList().map((x) => x.replaceAll(' ', '').split(',')).flat());

      const idInput = ui.textInput('', '');
      //@ts-ignore
      idInput.input.placeholder = 'Paste your IDs here...';
      idInput.input.style.resize = 'none';

      ui.dialog({title: 'Filter by IDs'})
        .add(ui.box(idInput.input))
        .onOK(() => {
          const query = idInput.stringValue.replaceAll(' ', '').split(',');
          const missing = query.filter((id) => !allIds.includes(id));

          if (missing.length > 0) {
            for (let i = 0; i < missing.length; i++)
              grok.shell.warning(missing[i] + ' not found in the base');
          }

          this.mlbTable.rows.filter((row) => {
            const additionalIds = row['gdb id mappings'].replaceAll(' ', '').split(',');
            let idsIntersect = false;
            for (let i = 0; i < additionalIds.length; i++) {
              if (query.includes(additionalIds[i])) {
                idsIntersect = true;
                break;
              }
            }
            return query.includes(row['v id']) || idsIntersect;
          });
          this.mlbTable.filter.fireChanged();
        })
        .show();
    }), 'Filter by multiple id query');
  }

  setQueryAntigenIcon(): void {
    const QueryAntigenIconHandler = () => {
      const agIds = this.mlbTable.getCol('antigen list').toList()
        .map((x) => x.replaceAll(' ', '').split(',')).flat();
      const agNcbiIds = this.mlbTable.getCol('antigen ncbi id').toList()
        .map((x) => x.replaceAll(' ', '').split(',')).flat();
      const geneSymbols = this.mlbTable.getCol('antigen gene symbol').toList()
        .map((x) => x.replaceAll(' ', '').split(',')).flat();
      const allIds = [...agIds, ...agNcbiIds, ...geneSymbols];

      const txtInput = ui.textInput('', '');
      //@ts-ignore
      txtInput.input.placeholder = 'Paste antigen id or gene symbol here ...';
      txtInput.input.style.resize = 'none';

      ui.dialog({title: 'Filter by antigen id / antigen gene symbol'})
        .add(ui.box(txtInput.input))
        .onOK(() => {
          const query = txtInput.stringValue.replaceAll(' ', '').split(',').filter((v) => v != '');
          const missing = query.filter((id) => !allIds.includes(id));

          if (missing.length > 0) {
            for (let i = 0; i < missing.length; i++)
              grok.shell.warning(`Value '${missing[i]}' not found in the base.`);
          }

          if (query.length > 0) {
            this.mlbTable.rows.filter((row) => {
              const rowIds = [].concat(
                row['antigen gene symbol'].replace(' ', '').split(','),
                row['antigen ncbi id'].replace(' ', '').split(','),
                row['antigen list'].replace(' ', '').split(',')
              ).filter((v) => v != '');
              return rowIds.some((rowId) => query.includes(rowId));
            });
          } else {
            this.mlbTable.rows.filter((row) => true);
          }
          this.mlbTable.filter.fireChanged();
        })
        .show();
    };

    this.queryAntigenIcon = ui.tooltip.bind(
      ui.iconFA('key', QueryAntigenIconHandler.bind(this)), 'Filter by antigen id / antigen gene symbol query');
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

  /** Calculate common positions sequence names for both Heavy and Light chains.
   * This is a draft.
   * @param {string[]} heavyPositions
   * @param {string[]} lightPositions
   */
  mergeSequencePositionsForImgt(heavyPositions: string[], lightPositions: string[]) {
    const res: string[] = [];
    let heavyI = 0;
    let lightI = 0;
    const positionRe = /^(\d+)([A-Z]*)$/g;

    while (heavyI < heavyPositions.length && lightI < lightPositions.length) {
      const heavyPositionName: string = heavyPositions[heavyI];
      const lightPositionName: string = lightPositions[lightI];
      const heavyPosM: RegExpMatchArray = heavyPositionName.match(positionRe)!;
      const lightPosM: RegExpMatchArray = lightPositionName.match(positionRe)!;

      const heavyPosNumber: number = parseInt(heavyPosM[0]);
      const lightPosNumber: number = parseInt(lightPosM[0]);

      const heavyPosPostfix: string = heavyPosM[1];
      const lightPosPostfix: string = lightPosM[1];

      if (heavyPosNumber == lightPosNumber) {
        if (heavyPosPostfix == lightPosPostfix) {
          res.push(heavyPositionName);
          heavyI++;
          lightI++;
        }
      }
    }
  }

  /** Builds multiple alignment sequences from monomers in positions
   *  as virtual (calculated) column in source DataFrame.
   * @param {DG.DataFrame} df  DataFrame with ANARCI results
   * @param {string} chain  Name of chain (used for result column name)
   * @param {number} startingColumnIndex  The first column with positions
   * @return {DG.Column}
   */
  mergeSequenceColumns(df: DG.DataFrame, chain: string, startingColumnIndex?: number) {
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

  addLogoViewer(table: DG.DataFrame, seqCol: DG.Column, dockType: DG.DOCK_TYPE, node?: DG.DockNode) {
    // Hide heavy & light chains amino acid sequences visually
    this.mlbView.grid.columns.byName(seqCol.name)!.visible = false;

    // const webLogo: WebLogo = new WebLogo();
    const logo = this.mlbView.addViewer('WebLogo', {sequenceColumnName: seqCol.name});
    return this.mlbView.dockManager.dock(logo, dockType, node, `${seqCol.name} chain`, 0.22);
  }

  onMLBGridCurrentRowChanged(args: any) {
    this.vIdInput.value = this.mlbTable.currentRow['v id'];
    this.changeVid(true);
  }

  async init(urlVid: string) {
    const pi = DG.TaskBarProgressIndicator.create('Loading data...');

    // #region Commented due to using code from RepertoireBrowser
    // const vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
    // #endregion

    ////////////////////////////////////////////////////
    this.mlbTable = await this.dataLoader.loadMlbDf();

    for (const column of this.mlbTable.columns)
      column.name = column.name.replaceAll('_', ' ');

    const mlbColumnsToMultiValue = ['antigen list', 'antigen ncbi id', 'antigen gene symbol'];
    for (const column of this.mlbTable.columns) {
      if (mlbColumnsToMultiValue.includes(column.name))
        column.setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, ',');
    }

    // let vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
    this.vids = this.dataLoader.vids;
    this.vidsObsPTMs = this.dataLoader.vidsObsPtm;
    ////////////////////////////////////////////////////

    // // #region Code from RepertoireBrowser till pi.close()
    // const referenceDf = this.dataLoader.refDf;
    //
    // const tempDf = referenceDf.clone(null, ['v_id']);
    // (tempDf.columns as DG.ColumnList).addNewInt('index').init((i) => i);
    // const indexes = (this.mlbTable.clone(null, ['v id'])
    //   .join(tempDf, ['v id'], ['v_id'], [], ['index'], 'left', false)
    //   .getCol('index').getRawData() as Int32Array);
    // // #endregion

    //#region -- Build columns with multiple alignments --
    /* Firstly get positions list for every chain.
       Secondly merge them to common positions sequence.
       Thirdly build sequences by common po
     */

    // commented out to remove old methods from DataLoader
    // const hChainDf = await this.dataLoader.loadHChainDf();
    // const lChainDf = await this.dataLoader.loadLChainDf();
    // [{name: 'Heavy', df: hChainDf}, {name: 'Light', df: lChainDf}].forEach((chain) => {
    //   const seqCol = this.mergeSequenceColumns(chain.df, chain.name);
    //   grok.data.joinTables(
    //     this.mlbTable, chain.df,
    //     ['v id'], ['Id'],
    //     (this.mlbTable.columns as DG.ColumnList).names(), [seqCol.name],
    //     DG.JOIN_TYPE.LEFT, true);
    //
    //   // crutch, because grok.data.joinTables() loses right table columns tags
    //   this.mlbTable.col(seqCol.name).setTag('positionNames', chain.df.col(seqCol.name).getTag('positionNames'));
    // });

    //#endregion -- Build columns with multiple alignments --

    pi.close();

    grok.events.onViewRemoved.subscribe((v) => {
      if (v.type === DG.VIEW_TYPE.TABLE_VIEW && (v as DG.TableView).dataFrame.id === this.mlbTable.id)
        this.subs.forEach((s) => s.unsubscribe());
    });

    await this.setView();


    // const hNode = this.addLogoViewer(hChainDf, 'Heavy', DG.DOCK_TYPE.TOP);
    // this.addLogoViewer(lChainDf, 'Light', DG.DOCK_TYPE.DOWN, hNode);

    this.allVids = this.mlbTable.col('v id')!;
    this.allIds = this.mlbTable.col('gdb id mappings')!;
    this.idMapping = {};
    for (let i = 0; i < this.allVids.length; i++)
      this.idMapping[this.allVids.get(i)] = this.allIds.get(i).replaceAll(' ', '').split(',');

    if (this.mlbView) {
      const tempDf: DG.DataFrame = DG.DataFrame.fromObjects([{}])!;
      this.regionsViewer = (await tempDf.plot.fromType(
        'VdRegions', {skipEmptyPositions: true
        })) as unknown as IVdRegionsViewer;
      await this.regionsViewer.init();
    }

    const dfTree: DG.DataFrame = await this.dataLoader.loadTreeDf();
    if (dfTree) {
      const treeBrowser = new TreeBrowserOld();
      await treeBrowser.init(dfTree, this.mlbView);
    }

    this.setVidInput();
    //this.setFilterIcon();
    this.setFilterIcon2();
    this.setQueryIdIcon();
    this.setQueryAntigenIcon();
    this.setSelectionClearIcon();
    this.setHideShowIcon();

    this.mlbView.setRibbonPanels([
      [this.vIdInput.root],
      [this.filterIcon, this.filterIcon2, this.queryIdIcon, this.queryAntigenIcon],
      [this.selectionClearIcon],
      [this.hideShowIcon, this.hideShowIcon2],
    ]);

    if (urlVid != null)
      this.vIdInput.value = urlVid;

    this.changeVid();
    this.setPropertiesFilters();

    this.mlbTable.onFilterChanged.subscribe((_) => {
      // this.compostionPviewer.do(this.mlbView);
    });
    this.mlbTable.onCurrentRowChanged.subscribe(this.onMLBGridCurrentRowChanged.bind(this));
  }
}
