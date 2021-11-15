import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import {MolecularLiabilityBrowserPanel} from "./molecular-liability-browser-panel.js"
import {NglMethods} from "./ngl.js"
import {PvizMethods} from "./pViz.js"
import {MiscMethods} from "./misc.js"

import properties from "./externalData/properties.json";
import json from "./examples/example.json";
import jsonNumbering from "./examples/exampleNums.json";
import jsonPDB from "./examples/examplePDB.json";

import {_package} from "./package.ts";

export class LaunchBrowser {

    async init(urlVid) {

        let pf = properties;
        //all continuous data operations
        let conts = ['cdr length', 'surface cdr hydrophobicity', 'positive cdr charge', 'negative cdr charge', 'sfvcsp'];
        let yellowLeft = [43, 116.35, -100, -100, -5.7];
        let yellowRight = [55, 174.4, 1.26, 1.86, 100];
        let redLeft = [37, 89.9, -100, -100, -19.5];
        let redRight = [63, 237.8, 3.73, 4.25, 100];
        let plots_x = [
            [-2.41379310344828,0,2.41379310344828,4.82758620689655,7.24137931034483,9.6551724137931,12.0689655172414,14.4827586206897,16.8965517241379,19.3103448275862,21.7241379310345,24.1379310344828,26.551724137931,28.9655172413793,31.3793103448276,33.7931034482759,36.2068965517241,38.6206896551724,41.0344827586207,43.448275862069,45.8620689655172,48.2758620689655,50.6896551724138,53.1034482758621,55.5172413793103,57.9310344827586,60.3448275862069,62.7586206896552,65.1724137931034,67.5862068965517,70,72.4137931034483],
            [-8.19835862068966,0,8.19835862068966,16.3967172413793,24.595075862069,32.7934344827586,40.9917931034483,49.1901517241379,57.3885103448276,65.5868689655172,73.7852275862069,81.9835862068965,90.1819448275862,98.3803034482759,106.578662068966,114.777020689655,122.975379310345,131.173737931034,139.372096551724,147.570455172414,155.768813793103,163.967172413793,172.165531034483,180.363889655172,188.562248275862,196.760606896552,204.958965517241,213.157324137931,221.355682758621,229.55404137931,237.7524,245.95075862069],
            [-2.41379310344828,0,2.41379310344828,4.82758620689655,7.24137931034483,9.6551724137931,12.0689655172414,14.4827586206897,16.8965517241379,19.3103448275862,21.7241379310345,24.1379310344828,26.551724137931,28.9655172413793,31.3793103448276,33.7931034482759,36.2068965517241,38.6206896551724,41.0344827586207,43.448275862069,45.8620689655172,48.2758620689655,50.6896551724138,53.1034482758621,55.5172413793103,57.9310344827586,60.3448275862069,62.7586206896552,65.1724137931034,67.5862068965517,70,72.4137931034483],
            [-2.41379310344828,0,2.41379310344828,4.82758620689655,7.24137931034483,9.6551724137931,12.0689655172414,14.4827586206897,16.8965517241379,19.3103448275862,21.7241379310345,24.1379310344828,26.551724137931,28.9655172413793,31.3793103448276,33.7931034482759,36.2068965517241,38.6206896551724,41.0344827586207,43.448275862069,45.8620689655172,48.2758620689655,50.6896551724138,53.1034482758621,55.5172413793103,57.9310344827586,60.3448275862069,62.7586206896552,65.1724137931034,67.5862068965517,70,72.4137931034483],
            [-21.6034482758621,-18.5172413793103,-15.4310344827586,-12.3448275862069,-9.25862068965517,-6.17241379310345,-3.08620689655173,-1.77635683940025e-15,3.08620689655172,6.17241379310344,9.25862068965517,12.3448275862069,15.4310344827586,18.5172413793103,21.6034482758621,24.6896551724138,27.7758620689655,30.8620689655172,33.948275862069,37.0344827586207,40.1206896551724,43.2068965517241,46.2931034482759,49.3793103448276,52.4655172413793,55.551724137931,58.6379310344828,61.7241379310345,64.8103448275862,67.8965517241379,70.9827586206897,74.0689655172414]   
        ];
        let plots_y = [
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000786120899972892,0.00157224179994578,0.0172946597994036,0.0416644076985634,0.122634860395771,0.0770398481973434,0.0652480346977502,0.0605313092979127,0.0149362970994849,0.0094334507996747,0.00235836269991868,0.000786120899972894,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0.000231452818133918,0.000694358454401754,0.00092581127253567,0.00786939581655321,0.0143500747243029,0.0222194705408561,0.0240710930859274,0.0185162254507134,0.0164331500875082,0.00717503736215143,0.00439760354454444,0.00138871690880351,0.00115726409066959,0.00138871690880351,0.000925811272535672,0,0,0,0.000231452818133918,0],
            [0,0.385199240986717,0.0283003523990241,0.000786120899972892,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0.341176470588235,0.0715370018975332,0.00157224179994578,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0.00184452948596991,0.00245937264795989,0.00122968632397994,0.0043039021339298,0.0104523337538295,0.0344312170714384,0.0842335131926261,0.0682475909808869,0.044883550825268,0.0301273149375086,0.021519510669649,0.00983749059183955,0.00368905897193983,0.00122968632397994,0.0043039021339298,0,0.00122968632397995,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        ];

        ////////////////////////////////////////////////////
        let ids = (await grok.data.loadTable(_package.webRoot + 'src/examples/mlb.csv'));
        for (let column of ids.columns)
            column.name = column.name.replaceAll("_", " ");
        // let vidsRaw = (await grok.functions.call('MolecularLiabilityBrowser:getVids'));
        let vids = ["VR000000008", "VR000000043","VR000000044"];
        ////////////////////////////////////////////////////

        //let path = _package.webRoot + 'files/' + 'Tree_demo.csv';
        //let treeData = (await grok.data.loadTable(path));
        //let treeBrowser = new TreeBrowser();

        //bands on plots for properties
        for(let i = 0;i < conts.length; i++){
            ids.col(conts[i]).setTag('.default-filter', '{ "min": ' + yellowLeft[i] + ', "max": ' + yellowRight[i] + ' }');
            ids.col(conts[i]).setTag('.charts', '[' +
            '{"type": "band", "title":"BandYellowLeft", "rule" : "' + redLeft[i] + '-' + yellowLeft[i] + '", "color": "#FFD700", "opacity": 0.15},' +
            '{"type": "band", "title":"BandYellowRight", "rule" : "' + yellowRight[i] + '-' + redRight[i] + '", "color": "#FFD700", "opacity": 0.15}, ' +
            '{"type": "band", "title":"BandRedLeft", "rule" : "< ' + redLeft[i] + '", "color": "#DC143C", "opacity": 0.15}, ' +
            '{"type": "band", "title":"BandRedRight", "rule" : "> ' + redRight[i] + '", "color": "#DC143C", "opacity": 0.15},' +
            '{ "type": "spline", "title": "TAP metrics", "y" : [' + plots_y[i].toString() + '], "color": "#7570B3", "width": 1, "x" : [' + plots_x[i].toString() + '], "normalize-y": true, "visible": true}' +
            ']');
        }

        //table visual polishing
        for (let column of ids.columns)
            column.name = column.name.replaceAll("_", " ");

        let view = grok.shell.addTableView(ids);
        view.name = "Molecular Liability Browser";

        for(let i = 0;i < conts.length; i++)
            view.grid.columns.byName(conts[i]).width = 150;

        view.grid.columns.byName('ngl').width = 100;
        view.grid.columns.byName("ngl").cellType = 'html';

        view.grid.onCellRender.subscribe(function (args) {
            if(args.cell.isColHeader){
              let textSize = args.g.measureText(args.cell.gridColumn.name);
              args.g.fillText(args.cell.gridColumn.name, args.bounds.x + (args.bounds.width - textSize.width)/2, args.bounds.y + (textSize.fontBoundingBoxAscent+textSize.fontBoundingBoxDescent));
              args.g.fillStyle = '#4b4b4a';
              args.preventDefault();
            }
        });

        view.grid.onCellPrepare(function (gc) {
            if (gc.isTableCell && gc.gridColumn.name === 'ngl' && vids.includes(gc.cell.value.toString())) {
              //debugger;
                gc.style.element = ui.divV([
                    ui.button("View", ()=> {
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

        let hideShowIcon = ui.iconFA('eye', ()=>{
            if(hideShowIcon.classList.value.includes('fa-eye-slash'))
              grok.events.fireCustomEvent("showAllDock");
            else
              grok.events.fireCustomEvent("closeAllDock");
        });
        hideShowIcon.classList = 'grok-icon fal fa-eye';

        let changeVid = () => {

            if (!!this.ngl) 
                this.ngl.stage.removeAllComponents();
            if (!!this.inputs && !!this.inputs.sequence_tabs) 
                view.dockManager.close(this.inputs.sequence_tabs);
            if (!!this.openPanels) 
                this.openPanels.forEach((p) =>  view.dockManager.close(p));
            if(!vids.includes(ci.value)){
                grok.shell.warning("No PDB data data for associated v id");
                return;
            }
            view.path = `/Table/${ci.value}`;
            hideShowIcon.classList = 'grok-icon fal fa-eye';
            this.init2(view, ci.value);
        };

        changeVid();

        ci.root.addEventListener("keyup", (event) => {
            // Number 13 is the "Enter" key on the keyboard
            if (event.keyCode === 13) 
              changeVid();
        });

        let nglFilterIcon = ui.tooltip.bind(ui.iconFA('filter', ()=>{ 
            ids.rows.select((row) =>  vids.includes(row['v id'].toString()));
            grok.functions.call('CmdSelectionToFilter');
        }), 'filter data for NGL viewer');
        
        // let treeViewIcon = ui.tooltip.bind(ui.iconFA('trees', async ()=>{ 
        //     if (treeData) {
        //         await treeBrowser.init(treeData);
        //      }
        //      else{
        //          grok.shell.error("Tree data not loaded");
        //      }
        // }), 'get tree view');


        //get all possible IDs
        let allIds = ids.col("v id").toList().concat(
            ids.col("gdb id mappings").toList().map(x => x.replaceAll(" ","").split(",")).flat()
        );

        let queryIcon = ui.tooltip.bind(ui.iconFA('layer-group', ()=>{ 
            let idInput = ui.textInput("","");
            idInput.input.placeholder = 'Paste your IDs here...';
            idInput.input.style.resize='none';
        
            ui.dialog('Filter by IDs','help url ...')
            .add(ui.box(idInput.input))
            .onOK(()=>{
                let query = idInput.stringValue.replaceAll(" ","").split(",");
                let missing = query.filter(id => !allIds.includes(id));

                if (missing.length > 0){
                  for(let i = 0; i < missing.length; i++){
                    grok.shell.warning(missing[i] + " not found in the base");
                  }
                }
        
                ids.rows.filter((row) => 
                {
                  let additionalIds = row['gdb id mappings'].replaceAll(" ","").split(",");
                  let idsIntersect = false;
                  for(let i = 0; i < additionalIds.length; i++){
                    if(query.includes(additionalIds[i])){
                      idsIntersect = true;
                      break;
                    }
                  }
                  return query.includes(row['v id']) || idsIntersect;
                });
                ids.filter.fireChanged();
            })
            .show();
        
          }), 'multiple id query');

        $(document.getElementsByClassName("d4-ribbon-panel")[6]).empty();
        $(document.getElementsByClassName("d4-ribbon-panel")[5]).empty();
        view.ribbonMenu.clear();
        view.setRibbonPanels([[ci.root], [nglFilterIcon, queryIcon, hideShowIcon]]);
        if(urlVid != null){
            ci.value = urlVid;
            changeVid();
        }
        
        grok.events.onCustomEvent("closeAllDock").subscribe((v) => {
            if (!!this.ngl) 
                this.ngl.stage.removeAllComponents();
            if (!!this.inputs && !!this.inputs.sequence_tabs) 
                view.dockManager.close(this.inputs.sequence_tabs);
            if (!!this.openPanels) 
                this.openPanels.forEach((p) =>  view.dockManager.close(p));

            hideShowIcon.classList = 'grok-icon fal fa-eye-slash';
        });
      
        grok.events.onCustomEvent("showAllDock").subscribe((v) => {
            changeVid();
            hideShowIcon.classList = 'grok-icon fal fa-eye';
        });

        grok.events.onTooltipShown.subscribe((args) => {
            if (args.args.context instanceof DG.Column) {
                switch(args.args.context.name){
                    case 'cdr length': 
                    args.args.element.innerHTML = "Complementarity-determining region length<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding, length is expressed in number of residuals<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;
                  
                    case 'surface cdr hydrophobicity': 
                    args.args.element.innerHTML = "Complementarity-determining region surface hydrophobicity<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding, hydrophobicity is expressed as polar surface in A<sup>2</sup><br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;

                    case 'positive cdr charge': 
                    args.args.element.innerHTML = "Complementarity-determining region positive charge<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;

                    case 'negative cdr charge': 
                    args.args.element.innerHTML = "Complementarity-determining region negative charge<br>" +
                    "cdr is a variable part of antibody that is reponsible for antigen binding<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;

                    case 'sfvcsp': 
                    args.args.element.innerHTML = "Structural Fv Charge Symmetry Parameter<br>" +
                    "is a measure of aggregate-inducing electrostatic attraction of the entire variable region<br>" +
                    "Gray bins - the whole distribution of value in dataset<br>" +
                    "Blue bins - value distribution for selected rows<br>" +
                    "Amber region - 0-5 and 95-100 quantiles for phase-I therapeutic Fvs<br>" + 
                    "Raspberry region - <0 and >100 quantiles for phase-I therapeutic Fvs<br>" +
                    "Spline - empirical density for phase-I therapeutic Fvs";
                    break;
                  
                    default:
                    break;
                }
            }
        });

    }

    
    async init2(view, vid) {

        let pi = DG.TaskBarProgressIndicator.create('Creating NGL view');
        
        ////////////////////////////////////////////////////
        let jsonStr = json;
        //let path = _package.webRoot + 'examples/example.pdb';
        let pdbStr = jsonPDB.pdb;
        let jsonN = jsonNumbering;

        //let jsonStr = JSON.parse((await grok.functions.call('MolecularLiabilityBrowser:getJsonByVid', {"vid" : vid})).columns[0].get(0));
        //let pdbStr = (await grok.functions.call('MolecularLiabilityBrowser:getPdbByVid', {"vid" : vid})).columns[0].get(0);
        //let jsonN = JSON.parse((await grok.functions.call('MolecularLiabilityBrowser:getJsonComplementByVid', {"vid" : vid})).columns[0].get(0));
        ////////////////////////////////////////////////////
        let hNumberingStr = jsonN.heavy_numbering;
        let lNumberingStr = jsonN.light_numbering;

        let hNumbering =[];
        let lNumbering =[];

        for(let i =0; i < hNumberingStr.length;i++)
            hNumbering.push(parseInt(hNumberingStr[i].replaceAll(" ", "")));
        
        for(let i =0; i < lNumberingStr.length;i++)
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
            MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs, this.inputs.paratopes);
        });
        this.inputs.paratopes.onChanged(async () => {
            await pViz.loadSequence(this.inputs, 'H', jsonStr)
            await pViz.loadSequence(this.inputs, 'L', jsonStr)
            MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs, this.inputs.paratopes);
        });
        this.inputs.ptm_choices.onChanged(async () => {
            pViz.pVizParams.ptmMap = pViz.ptmMapping(this.inputs.ptm_choices.value, this.inputs.ptm_prob.value, jsonStr)
            await pViz.loadSequence(this.inputs, 'H', jsonStr)
            await pViz.loadSequence(this.inputs, 'L', jsonStr)
            MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs, this.inputs.paratopes);
        });
        this.inputs.ptm_motif_choices.onChanged(async () => {
            pViz.pVizParams.ptmMotifsMap = pViz.ptmMotifsMapping(this.inputs.ptm_motif_choices.value, this.inputs.ptm_prob.value, jsonStr)
            await pViz.loadSequence(this.inputs, 'H', jsonStr)
            await pViz.loadSequence(this.inputs, 'L', jsonStr)
            MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs, this.inputs.paratopes);
        });
        this.inputs.ptm_prob.onChanged(async () => {
            pViz.pVizParams.ptmMap = pViz.ptmMapping(this.inputs.ptm_choices.value, this.inputs.ptm_prob.value, jsonStr)
            await pViz.loadSequence(this.inputs, 'H', jsonStr)
            await pViz.loadSequence(this.inputs, 'L',jsonStr )
            MiscMethods.setDockSize(view, this.inputs.ngl_node, this.inputs.sequence_tabs, this.inputs.paratopes);
        });
    
        pi.close();
    }
}