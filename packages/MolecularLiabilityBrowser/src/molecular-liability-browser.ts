import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {TwinPviewer} from './viewers/twin-p-viewer';
import {_package} from './package';
// import {TreeBrowser} from './tree';
import {Subscription} from 'rxjs';

export class MolecularLiabilityBrowser {
  mlbTable: DG.DataFrame;
  vids: string[];
  vidsObsPTMs: string[];
  allVids: DG.Column;
  allIds: DG.Column;
  idMapping: {[key: string]: string[];};

  mlbView: DG.TableView;
  vIdInput: DG.InputBase;
  filterIcon: HTMLElement;
  filterIcon2: HTMLElement;
  queryIcon: HTMLElement;
  hideShowIcon: HTMLElement;
  hideShowIcon2: HTMLElement;
  treesIcon: HTMLElement;

  twinPviewer: TwinPviewer;
  //compostionPviewer: CompostionPviewer;

  subs: Subscription[] = [];

  private changeVid = async (): Promise<void> => {
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

    this.mlbView.path = `/Table/${this.vIdInput.value}`;
    //hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    const pi = DG.TaskBarProgressIndicator.create('Creating NGL view');

    const pName = _package.nqName;
    const options = {'vid': this.vIdInput.value};
    const jsonStr = JSON.parse((await grok.functions.call(`${pName}:getJsonByVid`, options)).columns[0].get(0));
    const pdbStr = (await grok.functions.call(`${pName}:getPdbByVid`, options)).columns[0].get(0);
    const jsonN = JSON.parse((await grok.functions.call(`${pName}:getJsonComplementByVid`, options)).columns[0].get(0));
    let jsonStrObsPtm = null;
    if (this.vidsObsPTMs.includes(this.vIdInput.value))
      jsonStrObsPtm = JSON.parse((await grok.functions.call(`${pName}:getJsonObsByVid`, options)).columns[0].get(0));

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

    pi.close();
  };

  setPropertiesFilters(): void {
    interface PropertiesData {
      source: string;
      names: string[];
      yellowLeft: number[];
      yellowRight: number[];
      redLeft: number[];
      redRight: number[];
      plotsX: number[][];
      plotsY: number[][];
      tooltips: string[];
    };

    //external data load
    const pf = require('./externalData/properties.json') as PropertiesData;

    //this.mlbView.grid.columns.byName('CDR Clothia').visible = false;

    //bands on plots for properties
    for (let i = 0; i < pf.names.length; i++) {
      this.mlbTable.col(pf.names[i])!.setTag(
        '.default-filter', '{ "min": ' + pf.yellowLeft[i] + ', "max": ' + pf.yellowRight[i] + ' }');
      this.mlbTable.col(pf.names[i])!.setTag('.charts', '[' +
        '{"type": "band", "title":"BandYellowLeft", "rule" : "' +
        pf.redLeft[i] + '-' + pf.yellowLeft[i] + '", "color": "#FFD700", "opacity": 15},' +
        '{"type": "band", "title":"BandYellowRight", "rule" : "' +
        pf.yellowRight[i] + '-' + pf.redRight[i] + '", "color": "#FFD700", "opacity": 15}, ' +
        '{"type": "band", "title":"BandRedLeft", "rule" : "< ' + pf.redLeft[i] +
        '", "color": "#DC143C", "opacity": 15}, ' +
        '{"type": "band", "title":"BandRedRight", "rule" : "> ' + pf.redRight[i] +
        '", "color": "#DC143C", "opacity": 15},' +
        '{ "type": "spline", "title": "TAP metrics", "y" : [' +
        pf.plotsY[i].toString() + '], "color": "#7570B3", "width": 1, "x" : [' +
        pf.plotsX[i].toString() + '], "normalize-y": true, "visible": true}' +
        ']');
    }

    for (let i = 0; i < pf.names.length; i++)
      this.mlbTable.columns.byName(pf.names[i]).width = 150;


    //FIXME: filters appear separately
    this.mlbView.addViewer(DG.VIEWER.FILTERS);
    this.mlbView.filters({filters: [{type: 'MolecularLiabilityBrowser:mlbFilter'}]});

    grok.events.onTooltipShown.subscribe((args) => {
      if (args.args.context instanceof DG.Column) {
        switch (args.args.context.name) {
        case 'cdr length':
          args.args.element.innerHTML = pf.tooltips[0];
          break;

        case 'surface cdr hydrophobicity':
          args.args.element.innerHTML = pf.tooltips[1];
          break;

        case 'positive cdr charge':
          args.args.element.innerHTML = pf.tooltips[2];
          break;

        case 'negative cdr charge':
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

  setView(): void {
    this.mlbView = grok.shell.addTableView(this.mlbTable);

    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;
    this.mlbView.ribbonMenu.clear();

    this.mlbView.name = 'Molecular Liability Browser';
    for (const column of this.mlbTable.columns)
      column.name = column.name.replaceAll('_', ' ');

    this.mlbView.grid.columns.byName('v id')!.width = 120;
    this.mlbView.grid.columns.byName('v id')!.cellType = 'html';

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

  async init(urlVid) {
    const pi = DG.TaskBarProgressIndicator.create('Loading data...');

    this.mlbTable = (await grok.functions.call('MolecularLiabilityBrowser:GetMolecularLiabilityBrowser'));
    for (const column of this.mlbTable.columns)
      column.name = column.name.replaceAll('_', ' ');


    // 'ngl' column have been removed from query 2022-04
    // this.mlbTable.columns.remove('ngl');

    const vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
    this.vids = vidsRaw.columns[0].toList();
    this.vidsObsPTMs = (await grok.functions.call('MolecularLiabilityBrowser:getObservedPtmVids'))
      .columns[0].toList();

    pi.close();

    window.alert('Here we are!');

    grok.events.onViewRemoved.subscribe((v) => {
      if (v.type === DG.VIEW_TYPE.TABLE_VIEW && (v as DG.TableView).dataFrame.id === this.mlbTable.id)
        this.subs.forEach((s) => s.unsubscribe());
    });

    this.setView();
    this.allVids = this.mlbTable.col('v id')!;
    this.allIds = this.mlbTable.col('gdb id mappings')!;
    this.idMapping = {};
    for (let i = 0; i < this.allVids.length; i++)
      this.idMapping[this.allVids.get(i)] = this.allIds.get(i).replaceAll(' ', '').split(',');

    // const path = _package.webRoot + 'src/examples/tree.csv';
    // const dfTree = await grok.data.loadTable(path);
    // if (dfTree) {
    //   const treeBrowser = new TreeBrowser();
    //   await treeBrowser.init(dfTree, this.mlbView);
    // }

    // this.setVidInput();
    // //this.setFilterIcon();
    // this.setFilterIcon2();
    // this.setQueryIcon();
    // this.setHideShowIcon();

    // this.mlbView.setRibbonPanels([[this.vIdInput.root],
    //   [this.filterIcon, this.filterIcon2, this.queryIcon, this.hideShowIcon, this.hideShowIcon2]]);

    // if (urlVid != null)
    //   this.vIdInput.value = urlVid;


    // this.changeVid();
    // this.setPropertiesFilters();

    // this.mlbTable.onFilterChanged.subscribe((_) => {
    //   this.compostionPviewer.do(this.mlbView);
    // });
  }
}
