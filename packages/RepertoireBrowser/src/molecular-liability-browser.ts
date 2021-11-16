import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { TwinPviewer } from "./viewers/twin-p-viewer"

import { _package } from "./package";

export class MolecularLiabilityBrowser {

  mlbTable: DG.DataFrame;
  vids: string[];

  mlbView: DG.TableView;
  vIdInput: DG.InputBase;
  filterIcon: HTMLElement;
  queryIcon: HTMLElement;
  hideShowIcon: HTMLElement;

  twinPviewer: TwinPviewer;

  #changeVid = (): void => {
    if (!this.vids.includes(this.vIdInput.value)) {
      grok.shell.warning("No PDB data data for associated v id");
      return;
    }

    if(!this.twinPviewer){
      this.twinPviewer = new TwinPviewer();
    } else {
      this.twinPviewer.reset(this.mlbView);
    }


    this.mlbView.path = `/Table/${this.vIdInput.value}`;
    //hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    let pi = DG.TaskBarProgressIndicator.create('Creating NGL view');

    ////////////////////////////////////////////////////
    let jsonStr = require("./examples/example.json");
    //let path = _package.webRoot + 'examples/example.pdb';
    let pdbStr: string = require("./examples/examplePDB.json").pdb;
    let jsonN = require("./examples/exampleNums.json");

    //let jsonStr = JSON.parse((await grok.functions.call('MolecularLiabilityBrowser:getJsonByVid', {"vid" : vid})).columns[0].get(0));
    //let pdbStr = (await grok.functions.call('MolecularLiabilityBrowser:getPdbByVid', {"vid" : vid})).columns[0].get(0);
    //let jsonN = JSON.parse((await grok.functions.call('MolecularLiabilityBrowser:getJsonComplementByVid', {"vid" : vid})).columns[0].get(0));
    ////////////////////////////////////////////////////
    let hNumberingStr = jsonN.heavy_numbering;
    let lNumberingStr = jsonN.light_numbering;

    let hNumbering = [];
    let lNumbering = [];

    for (let i = 0; i < hNumberingStr.length; i++)
      hNumbering.push(parseInt(hNumberingStr[i].replaceAll(" ", "")));

    for (let i = 0; i < lNumberingStr.length; i++)
      lNumbering.push(parseInt(lNumberingStr[i].replaceAll(" ", "")));

    jsonStr["map_H"] = hNumbering;
    jsonStr["map_L"] = lNumbering;

    this.twinPviewer.twin(this.mlbView, jsonStr, pdbStr);

    pi.close();
  };

  #setPropertiesFilters = (): void => {
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
    const pf = require("./externalData/properties.json") as PropertiesData;

    //bands on plots for properties  
    for (let i = 0; i < pf.names.length; i++) {

      this.mlbTable.col(pf.names[i])!.setTag('.default-filter', '{ "min": ' + pf.yellowLeft[i] + ', "max": ' + pf.yellowRight[i] + ' }');
      this.mlbTable.col(pf.names[i])!.setTag('.charts', '[' +
        '{"type": "band", "title":"BandYellowLeft", "rule" : "' + pf.redLeft[i] + '-' + pf.yellowLeft[i] + '", "color": "#FFD700", "opacity": 0.15},' +
        '{"type": "band", "title":"BandYellowRight", "rule" : "' + pf.yellowRight[i] + '-' + pf.redRight[i] + '", "color": "#FFD700", "opacity": 0.15}, ' +
        '{"type": "band", "title":"BandRedLeft", "rule" : "< ' + pf.redLeft[i] + '", "color": "#DC143C", "opacity": 0.15}, ' +
        '{"type": "band", "title":"BandRedRight", "rule" : "> ' + pf.redRight[i] + '", "color": "#DC143C", "opacity": 0.15},' +
        '{ "type": "spline", "title": "TAP metrics", "y" : [' + pf.plotsY[i].toString() + '], "color": "#7570B3", "width": 1, "x" : [' + pf.plotsX[i].toString() + '], "normalize-y": true, "visible": true}' +
        ']');
    }

    for (let i = 0; i < pf.names.length; i++)
      this.mlbTable.columns.byName(pf.names[i]).width = 150;

    //table visual polishing
    for (let column of this.mlbTable.columns)
      column.name = column.name.replaceAll("_", " ");

    this.mlbView.addViewer(DG.VIEWER.FILTERS);

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

  #setView = (): void => {
    this.mlbView = grok.shell.addTableView(this.mlbTable);

    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;
    this.mlbView.ribbonMenu.clear();

    this.mlbView.name = "Molecular Liability Browser";
    this.mlbView.grid.columns.byName('ngl')!.width = 100;
    this.mlbView.grid.columns.byName("ngl")!.cellType = 'html';

    this.mlbView.grid.onCellRender.subscribe(function (args) {
      if (args.cell.isColHeader) {
        let textSize = args.g.measureText(args.cell.gridColumn.name);
        args.g.fillText(args.cell.gridColumn.name, args.bounds.x + (args.bounds.width - textSize.width) / 2, args.bounds.y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent));
        args.g.fillStyle = '#4b4b4a';
        args.preventDefault();
      }
    });

    this.mlbView.grid.onCellPrepare((gc) => {
      if (gc.isTableCell && gc.gridColumn.name === 'ngl' && this.vids.includes(gc.cell.value.toString())) {
        //debugger;
        gc.style.element = ui.divV([
          ui.button("View", () => {
            this.vIdInput.value = gc.cell.value;
            this.#changeVid();
          })
        ]);
      }
    });
  }

  #setVidInput = (): void => {
    this.vIdInput = ui.stringInput('VID', '');
    this.vIdInput.value = this.vids[0];

    this.vIdInput.root.addEventListener("keyup", (event) => {
      if (event.key === 'Enter')
        this.#changeVid();
    });
  }

  #setFilterIcon = (): void => {
    this.filterIcon = ui.tooltip.bind(ui.iconFA('filter', () => {
      this.mlbTable.rows.select((row) => this.vids.includes(row.get('v id').toString()));
      grok.functions.call('CmdSelectionToFilter');
    }), 'filter data for NGL viewer');
  }

  #setQueryIcon = (): void => {
    this.queryIcon = ui.tooltip.bind(ui.iconFA('layer-group', () => {
      //get all possible IDs
      let allIds = this.mlbTable.col("v id")!.toList().concat(
        this.mlbTable.col("gdb id mappings")!.toList().map(x => x.replaceAll(" ", "").split(",")).flat()
      );

      let idInput = ui.textInput("", "");
      //@ts-ignore
      idInput.input.placeholder = 'Paste your IDs here...';
      idInput.input.style.resize = 'none';

      ui.dialog({ title: 'Filter by IDs' })
        .add(ui.box(idInput.input))
        .onOK(() => {
          let query = idInput.stringValue.replaceAll(" ", "").split(",");
          let missing = query.filter(id => !allIds.includes(id));

          if (missing.length > 0) {
            for (let i = 0; i < missing.length; i++) {
              grok.shell.warning(missing[i] + ' not found in the base');
            }
          }

          this.mlbTable.rows.filter((row) => {
            let additionalIds = row.get('gdb id mappings').replaceAll(" ", "").split(",");
            let idsIntersect = false;
            for (let i = 0; i < additionalIds.length; i++) {
              if (query.includes(additionalIds[i])) {
                idsIntersect = true;
                break;
              }
            }
            return query.includes(row.get('v id')) || idsIntersect;
          });
          this.mlbTable.filter.fireChanged();
        })
        .show();

    }), 'multiple id query');
  }

  #setHideShowIcon = (): void => {
    this.hideShowIcon = ui.iconFA('eye', () => {
      if (this.hideShowIcon.classList.value.includes('fa-eye-slash'))
        grok.events.fireCustomEvent("showAllDock", null);
      else
        grok.events.fireCustomEvent("closeAllDock", null);
    });
    this.hideShowIcon.classList.value = 'grok-icon fal fa-eye';

    grok.events.onCustomEvent("closeAllDock").subscribe((v) => {
      this.twinPviewer.reset(this.mlbView);
      this.hideShowIcon.classList.value = 'grok-icon fal fa-eye-slash';
    });

    grok.events.onCustomEvent("showAllDock").subscribe((v) => {
      this.#changeVid();
      this.hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    });
  }

  async init(urlVid: string | null) {

    ////////////////////////////////////////////////////
    this.mlbTable = (await grok.data.loadTable(_package.webRoot + 'src/examples/mlb.csv'));
    for (let column of this.mlbTable.columns)
      column.name = column.name.replaceAll("_", " ");
    // let vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
    this.vids = ["VR000000008", "VR000000043", "VR000000044"];
    ////////////////////////////////////////////////////

    this.#setView();
    this.#setPropertiesFilters();
    this.#setVidInput();
    this.#setFilterIcon();
    this.#setQueryIcon();
    this.#setHideShowIcon();

    this.mlbView.setRibbonPanels([[this.vIdInput.root], [this.filterIcon, this.queryIcon, this.hideShowIcon]]);
    
    if (urlVid != null) 
      this.vIdInput.value = urlVid;

    this.#changeVid();
  }
}
