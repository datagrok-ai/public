import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { MolecularLiabilityBrowserPanel } from "./molecular-liability-browser-panel.js"
import { NglMethods } from "./ngl.js"
import { PvizMethods } from "./pViz.js"
import { MiscMethods } from "./misc.js"


import { _package } from "./package";

export class MolecularLiabilityBrowser {

  inputs;
  openPanels;
  ngl;
  pViz;

  #setPropertiesFilters = function (mlbTable) {
    //external data load
    const pf = require("./externalData/properties.json");

    //bands on plots for properties
    for (let i = 0; i < pf.names.length; i++) {
      mlbTable.col(pf.names[i]).setTag('.default-filter', '{ "min": ' + pf.yellowLeft[i] + ', "max": ' + pf.yellowRight[i] + ' }');
      mlbTable.col(pf.names[i]).setTag('.charts', '[' +
        '{"type": "band", "title":"BandYellowLeft", "rule" : "' + pf.redLeft[i] + '-' + pf.yellowLeft[i] + '", "color": "#FFD700", "opacity": 0.15},' +
        '{"type": "band", "title":"BandYellowRight", "rule" : "' + pf.yellowRight[i] + '-' + pf.redRight[i] + '", "color": "#FFD700", "opacity": 0.15}, ' +
        '{"type": "band", "title":"BandRedLeft", "rule" : "< ' + pf.redLeft[i] + '", "color": "#DC143C", "opacity": 0.15}, ' +
        '{"type": "band", "title":"BandRedRight", "rule" : "> ' + pf.redRight[i] + '", "color": "#DC143C", "opacity": 0.15},' +
        '{ "type": "spline", "title": "TAP metrics", "y" : [' + pf.plotsY[i].toString() + '], "color": "#7570B3", "width": 1, "x" : [' + pf.plotsX[i].toString() + '], "normalize-y": true, "visible": true}' +
        ']');
    }

    for (let i = 0; i < pf.names.length; i++)
      mlbTable.columns.byName(pf.names[i]).width = 150;

    //table visual polishing
    for (let column of mlbTable.columns)
      column.name = column.name.replaceAll("_", " ");

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

  async init(urlVid) {

    ////////////////////////////////////////////////////
    let mlbTable = (await grok.data.loadTable(_package.webRoot + 'src/examples/mlb.csv'));
    for (let column of mlbTable.columns)
      column.name = column.name.replaceAll("_", " ");
    // let vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
    let vids = ["VR000000008", "VR000000043", "VR000000044"];
    ////////////////////////////////////////////////////

    this.#setPropertiesFilters(mlbTable);

    let a : number = 5;
    let view = grok.shell.addTableView(mlbTable);
    view.name = "Molecular Liability Browser";

    view.grid.columns.byName('ngl').width = 100;
    view.grid.columns.byName("ngl").cellType = 'html';

    view.grid.onCellRender.subscribe(function (args) {
      if (args.cell.isColHeader) {
        let textSize = args.g.measureText(args.cell.gridColumn.name);
        args.g.fillText(args.cell.gridColumn.name, args.bounds.x + (args.bounds.width - textSize.width) / 2, args.bounds.y + (textSize.fontBoundingBoxAscent + textSize.fontBoundingBoxDescent));
        args.g.fillStyle = '#4b4b4a';
        args.preventDefault();
      }
    });

    view.grid.onCellPrepare(function (gc) {
      if (gc.isTableCell && gc.gridColumn.name === 'ngl' && vids.includes(gc.cell.value.toString())) {
        //debugger;
        gc.style.element = ui.divV([
          ui.button("View", () => {
            ci.value = gc.cell.value;
            changeVid();
          })
        ]);
      }
    });

    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;

    view.addViewer(DG.VIEWER.FILTERS);

    let ci = ui.stringInput('VID', '');
    ci.value = vids[0];;

    let hideShowIcon = ui.iconFA('eye', () => {
      if (hideShowIcon.classList.value.includes('fa-eye-slash'))
        grok.events.fireCustomEvent("showAllDock", null);
      else
        grok.events.fireCustomEvent("closeAllDock", null);
    });
    hideShowIcon.classList.value = 'grok-icon fal fa-eye';

    let changeVid = () => {

      if (!!this.ngl)
        this.ngl.stage.removeAllComponents();
      if (!!this.inputs && !!this.inputs.sequence_tabs)
        view.dockManager.close(this.inputs.sequence_tabs);
      if (!!this.openPanels)
        this.openPanels.forEach((p) => view.dockManager.close(p));
      if (!vids.includes(ci.value)) {
        grok.shell.warning("No PDB data data for associated v id");
        return;
      }
      view.path = `/Table/${ci.value}`;
      hideShowIcon.classList.value = 'grok-icon fal fa-eye';
      this.init2(view, ci.value);
    };

    changeVid();

    ci.root.addEventListener("keyup", (event) => {
      if (event.key === 'Enter')
        changeVid();
    });

    let nglFilterIcon = ui.tooltip.bind(ui.iconFA('filter', () => {
      mlbTable.rows.select((row) => vids.includes(row['v id'].toString()));
      grok.functions.call('CmdSelectionToFilter');
    }), 'filter data for NGL viewer');

    //get all possible IDs
    let allIds = mlbTable.col("v id").toList().concat(
      mlbTable.col("gdb id mappings").toList().map(x => x.replaceAll(" ", "").split(",")).flat()
    );

    let queryIcon = ui.tooltip.bind(ui.iconFA('layer-group', () => {
      let idInput = ui.textInput("", "");
      //idInput.input.placeholder = 'Paste your IDs here...';
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

          mlbTable.rows.filter((row) => {
            let additionalIds = row['gdb id mappings'].replaceAll(" ", "").split(",");
            let idsIntersect = false;
            for (let i = 0; i < additionalIds.length; i++) {
              if (query.includes(additionalIds[i])) {
                idsIntersect = true;
                break;
              }
            }
            return query.includes(row['v id']) || idsIntersect;
          });
          mlbTable.filter.fireChanged();
        })
        .show();

    }), 'multiple id query');

    view.ribbonMenu.clear();
    view.setRibbonPanels([[ci.root], [nglFilterIcon, queryIcon, hideShowIcon]]);
    if (urlVid != null) {
      ci.value = urlVid;
      changeVid();
    }

    grok.events.onCustomEvent("closeAllDock").subscribe((v) => {
      if (!!this.ngl)
        this.ngl.stage.removeAllComponents();
      if (!!this.inputs && !!this.inputs.sequence_tabs)
        view.dockManager.close(this.inputs.sequence_tabs);
      if (!!this.openPanels)
        this.openPanels.forEach((p) => view.dockManager.close(p));

      hideShowIcon.classList.value = 'grok-icon fal fa-eye-slash';
    });

    grok.events.onCustomEvent("showAllDock").subscribe((v) => {
      changeVid();
      hideShowIcon.classList.value = 'grok-icon fal fa-eye';
    });

  }


  async init2(view, vid) {

    let pi = DG.TaskBarProgressIndicator.create('Creating NGL view');

    ////////////////////////////////////////////////////
    let jsonStr = require("./examples/example.json");
    //let path = _package.webRoot + 'examples/example.pdb';
    let pdbStr = require("./examples/examplePDB.json").pdb;
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


    this.inputs = new MolecularLiabilityBrowserPanel();
    this.openPanels = await this.inputs.init(view, jsonStr);
    this.ngl = new NglMethods();
    await this.ngl.init(view, this.inputs, pdbStr, jsonStr);
    let pViz = new PvizMethods();
    this.openPanels.push(await pViz.init(view, this.inputs, this.ngl, jsonStr));


    this.inputs.repChoice.onChanged(async () => {
      await pViz.loadSequence(this.inputs, 'H', jsonStr, true)
      await pViz.loadSequence(this.inputs, 'L', jsonStr, true)
    });
    this.inputs.cdr_scheme.onChanged(async () => {
      pViz.pVizParams.cdrMap = pViz.cdrMapping(this.inputs.cdr_scheme.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.paratopes.onChanged(async () => {
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.ptm_choices.onChanged(async () => {
      pViz.pVizParams.ptmMap = pViz.ptmMapping(this.inputs.ptm_choices.value, this.inputs.ptm_prob.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.ptm_motif_choices.onChanged(async () => {
      pViz.pVizParams.ptmMotifsMap = pViz.ptmMotifsMapping(this.inputs.ptm_motif_choices.value, this.inputs.ptm_prob.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.ptm_prob.onChanged(async () => {
      pViz.pVizParams.ptmMap = pViz.ptmMapping(this.inputs.ptm_choices.value, this.inputs.ptm_prob.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });

    pi.close();
  }
}