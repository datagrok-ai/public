import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {TwinPviewer} from './viewers/twin-p-viewer';
import {TreeBrowser} from './tree';
import {_package} from './package';
import {CustomFilter, StringMap} from './utils/filters';
import {AminoacidsWebLogo} from './viewers/web-logo';
import {genData} from './utils/example-generator';

import $ from 'cash-dom';

export class MolecularLiabilityBrowser {
  mlbTable: DG.DataFrame;
  vids: string[];
  vidsObsPTMs: string[];
  allVids: DG.Column;
  allIds: DG.Column;
  idMapping: { [key: string]: string[]; }

  mlbView: DG.TableView;
  vIdInput: DG.InputBase;
  filterIcon: HTMLElement;
  filterIcon2: HTMLElement;
  queryIcon: HTMLElement;
  hideShowIcon: HTMLElement;
  hideShowIcon2: HTMLElement;
  treesIcon: HTMLElement;

  twinPviewer: TwinPviewer;

  pf: PropertiesData;

  constructor() {
    //external data load
    this.pf = require('./externalData/properties.json') as PropertiesData;
  }

  changeVid(silent = false) {
    if (!this.vids.includes(this.vIdInput.value)) {
      Object.keys(this.idMapping).every((vid) => {
        if (this.idMapping[vid].includes(this.vIdInput.value)) {
          this.vIdInput.value = vid;
          return false;
        }
        return true;
      });

      if (!this.vids.includes(this.vIdInput.value)) {
        if (!silent)
          grok.shell.warning('No PDB data data for associated v id');
        return;
      }
    }

    this.mlbView.path = `/Table/${this.vIdInput.value}`;
    //hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

    ////////////////////////////////////////////////////
    const jsonStr = require('./examples/example.json');
    //let path = _package.webRoot + 'examples/example.pdb';
    const pdbStr: string = require('./examples/examplePDB.json').pdb;
    const jsonN = require('./examples/exampleNums.json');

    let jsonStrObsPtm = null;
    if (this.vidsObsPTMs.includes(this.vIdInput.value))
      jsonStrObsPtm = require('./examples/exampleOptm.json');
    ////////////////////////////////////////////////////
    const hNumberingStr = jsonN.heavy_numbering;
    const lNumberingStr = jsonN.light_numbering;

    const hNumbering = [];
    const lNumbering = [];

    for (let i = 0; i < hNumberingStr.length; i++)
      hNumbering.push(parseInt(hNumberingStr[i].replaceAll(' ', '')));

    for (let i = 0; i < lNumberingStr.length; i++)
      lNumbering.push(parseInt(lNumberingStr[i].replaceAll(' ', '')));

    jsonStr['map_H'] = hNumbering;
    jsonStr['map_L'] = lNumbering;

    if (!this.twinPviewer) {
      this.twinPviewer = new TwinPviewer();
      this.twinPviewer.init(jsonStr, pdbStr, jsonStrObsPtm);
    } else
      this.twinPviewer.reset(jsonStr, pdbStr, jsonStrObsPtm);

    this.twinPviewer.show(this.mlbView);
    this.twinPviewer.open(this.mlbView);
    this.hideShowIcon.classList.value = 'grok-icon fal fa-eye';

    pi.close();
  }

  setVidInput(): void {
    this.vIdInput = ui.stringInput('VID', '');
    this.vIdInput.value = this.vids[0];

    this.vIdInput.root.addEventListener('keyup', (event) => {
      if (event.key === 'Enter')
        this.changeVid();
    });
  }

  setFilterIcon(): void {
    this.filterIcon = ui.tooltip.bind(ui.iconFA('filter', () => {
      const a = this.mlbTable;
      a.rows.select((row) => this.vids.includes(row['v id'].toString()));
      grok.functions.call('CmdSelectionToFilter');
    }), 'filter data for 3D view');
  }

  setFilterIcon2(): void {
    this.filterIcon2 = ui.tooltip.bind(ui.iconFA('filter', () => {
      const a = this.mlbTable;
      a.rows.select((row) => this.vidsObsPTMs.includes(row['v id'].toString()));
      grok.functions.call('CmdSelectionToFilter');
    }), 'filter data with observed PTMs');
  }

  setQueryIcon(): void {
    this.queryIcon = ui.tooltip.bind(ui.iconFA('layer-group', () => {
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
    }), 'multiple id query');
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

  setHideShowIcon2(): void {
    this.hideShowIcon2 = ui.iconFA('eye', () => {
      if (this.hideShowIcon2.classList.value.includes('fa-eye-slash'))
        grok.events.fireCustomEvent('showAllDock2', null);
      else
        grok.events.fireCustomEvent('closeAllDock2', null);
    });
    this.hideShowIcon2.classList.value = 'grok-icon fal fa-eye-slash';
  }

  onFiltersTooltipShown(args: any) {
    if (args.args.context instanceof DG.Column) {
      switch (args.args.context.name) {
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
        break;
      }
    }
  }

  setPropertiesFilters(
    ptmMap: StringMap,
    cdrMap: StringMap,
    referenceDf: DG.DataFrame,
    indexes: Int32Array,
  ): void {
    //this.mlbView.grid.columns.byName('CDR Clothia').visible = false;

    //bands on plots for properties
    for (let i = 0; i < this.pf.names.length; i++) {
      this.mlbTable.col(this.pf.names[i])!.setTag(
        '.default-filter', '{ "min": ' + this.pf.yellowLeft[i] + ', "max": ' + this.pf.yellowRight[i] + ' }');
      this.mlbTable.col(this.pf.names[i])!.setTag('.charts', '[' +
        '{"type": "band", "title":"BandYellowLeft", "rule" : "' +
        this.pf.redLeft[i] + '-' + this.pf.yellowLeft[i] + '", "color": "#FFD700", "opacity": 15},' +
        '{"type": "band", "title":"BandYellowRight", "rule" : "' +
        this.pf.yellowRight[i] + '-' + this.pf.redRight[i] + '", "color": "#FFD700", "opacity": 15}, ' +
        '{"type": "band", "title":"BandRedLeft", "rule" : "< ' + this.pf.redLeft[i] +
          '", "color": "#DC143C", "opacity": 15}, ' +
        '{"type": "band", "title":"BandRedRight", "rule" : "> ' + this.pf.redRight[i] +
          '", "color": "#DC143C", "opacity": 15},' +
        '{ "type": "spline", "title": "TAP metrics", "y" : [' +
        this.pf.plotsY[i].toString() + '], "color": "#7570B3", "width": 1, "x" : [' +
        this.pf.plotsX[i].toString() + '], "normalize-y": true, "visible": true}' +
        ']');
    }

    for (let i = 0; i < this.pf.names.length; i++)
      this.mlbTable.columns.byName(this.pf.names[i]).width = 150;

    const filters = this.mlbView.addViewer(DG.VIEWER.FILTERS);
    const cFilter = new CustomFilter(this.mlbTable, ptmMap, cdrMap, referenceDf, indexes).create();

    $(this.mlbView.root).ready(() => $(filters.root).ready(() => filters.root.append(cFilter.root)));

    grok.events.onTooltipShown.subscribe(this.onFiltersTooltipShown.bind(this));
  }

  setView(vIdColName = 'v id'): void {
    this.mlbView = grok.shell.addTableView(this.mlbTable);

    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;
    this.mlbView.ribbonMenu.clear();

    this.mlbView.name = 'Molecular Liability Browser';
    for (const column of this.mlbTable.columns)
      column.name = column.name.replaceAll('_', ' ');
    this.mlbView.grid.columns.byName(vIdColName)!.width = 120;
    this.mlbView.grid.columns.byName(vIdColName)!.cellType = 'html';

    //table visual polishing

    this.mlbView.grid.onCellRender.subscribe(function(args) {
      if (args.cell.isColHeader) {
        const textSize = args.g.measureText(args.cell.gridColumn.name);
        args.g.fillText(args.cell.gridColumn.name, args.bounds.x +
          (args.bounds.width - textSize.width) / 2, args.bounds.y +
          (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent));
        args.g.fillStyle = '#4b4b4a';
        args.preventDefault();
      }
    });

    this.mlbView.grid.onCellPrepare((gc) => {
      if (gc.isTableCell && gc.gridColumn.name === 'v id') {
        if (this.vids.includes(gc.cell.value.toString())) {
          gc.style.element = ui.divV([
            ui.link(gc.cell.value.toString(), () => {
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

  async readSequences() {
    //this.seqTable = (await grok.data.loadTable(_package.webRoot + 'src/examples/h_out1.csv'));
    const [hChainDf, lChainDf] = genData(this.mlbTable.col('v id').toList());
    return [hChainDf, lChainDf];
  }

  mergeSequenceColumns(df: DG.DataFrame, chain: string, startingColumnIndex = 13) {
    const positionRegExp = /^\d+[A-Z]*$/g;
    const columns: DG.ColumnList = df.columns;
    const names = columns.names().slice(startingColumnIndex);
    const positionColumns = names.filter((v: string) => v.match(positionRegExp) !== null);
    const seqCol = columns.addNewVirtual(
      `${chain} chain sequence`,
      (i: number) => positionColumns.map((v) => df.get(v, i)).join(''),
    );
    seqCol.semType = AminoacidsWebLogo.residuesSet;
    return seqCol;
  }

  addLogoViewer(table: DG.DataFrame, chain: string, dockType: DG.DOCK_TYPE, node?: DG.DockNode) {
    const seqCol = this.mergeSequenceColumns(table, chain);

    grok.data.joinTables(
      this.mlbTable,
      table,
      ['v id'],
      ['Id'],
      (this.mlbTable.columns as DG.ColumnList).names(),
      [seqCol.name],
      DG.JOIN_TYPE.LEFT,
      true,
    );

    const logo = this.mlbView.addViewer('AminoacidsWebLogo');
    return this.mlbView.dockManager.dock(logo, dockType, node, `${chain} chain`, 0.2);
  }

  onMLBGridCurrentRowChanged(args: any) {
    this.vIdInput.value = this.mlbTable.currentRow['v id'];
    this.changeVid(true);
  }

  async init(urlVid: string | null) {
    const pi = DG.TaskBarProgressIndicator.create('Loading data...');

    ////////////////////////////////////////////////////
    this.mlbTable = (await grok.data.loadTable(_package.webRoot + 'src/examples/mlb.csv'));
    this.mlbTable.columns.remove('ngl');

    const [hChainDf, lChainDf] = await this.readSequences();

    for (const column of this.mlbTable.columns)
      column.name = column.name.replaceAll('_', ' ');

    // let vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
    this.vids = ['VR000000008', 'VR000000043', 'VR000000044'];
    this.vidsObsPTMs = ['VR000000044'];
    ////////////////////////////////////////////////////

    const ptmMap = JSON.parse(await _package.files.readAsText('ptm_map.json'));
    const cdrMap = JSON.parse(await _package.files.readAsText('cdr_map.json'));
    const referenceDf = (await _package.files.readBinaryDataFrames('ptm_in_cdr.d42'))[0];

    const tempDf = referenceDf.clone(null, ['v_id']);
    (tempDf.columns as DG.ColumnList).addNewInt('index').init((i) => i);
    const indexes = (this.mlbTable.clone(null, ['v id']).
      join(tempDf, ['v id'], ['v_id'], [], ['index'], 'left', false).getCol('index').getRawData() as Int32Array);
    pi.close();

    this.setView();

    const hNode = this.addLogoViewer(hChainDf, 'Heavy', DG.DOCK_TYPE.TOP);
    this.addLogoViewer(lChainDf, 'Light', DG.DOCK_TYPE.DOWN, hNode);

    this.allVids = this.mlbTable.col('v id')!;
    this.allIds = this.mlbTable.col('gdb id mappings')!;
    this.idMapping = {};
    for (let i = 0; i < this.allVids.length; i++)
      this.idMapping[this.allVids.get(i)] = this.allIds.get(i).replaceAll(' ', '').split(',');

    const path = _package.webRoot + 'src/examples/tree.csv';
    const dfTree = await grok.data.loadTable(path);

    if (dfTree) {
      const treeBrowser = new TreeBrowser();
      treeBrowser.init(dfTree, this.mlbView);
    }

    this.setVidInput();
    this.setFilterIcon();
    this.setFilterIcon2();
    this.setQueryIcon();
    this.setHideShowIcon();
    this.setHideShowIcon2();

    this.mlbView.setRibbonPanels([[this.vIdInput.root],
      [this.filterIcon, this.filterIcon2, this.queryIcon, this.hideShowIcon, this.hideShowIcon2]]);

    if (urlVid != null)
      this.vIdInput.value = urlVid;

    this.changeVid();
    this.setPropertiesFilters(ptmMap, cdrMap, referenceDf, indexes);
    this.mlbTable.onCurrentRowChanged.subscribe(this.onMLBGridCurrentRowChanged.bind(this));
  }
}

type PropertiesData = {
  source: string;
  names: string[];
  yellowLeft: number[];
  yellowRight: number[];
  redLeft: number[];
  redRight: number[];
  plotsX: number[][];
  plotsY: number[][];
  tooltips: string[];
}
