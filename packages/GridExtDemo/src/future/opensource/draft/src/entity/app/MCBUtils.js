import {GridUtils} from "../../utils/GridUtils";
import {SemType} from "../SemType";
import {DGApp} from "./DGApp";
//import * as jt from "jsTree";
var $ = require( "jquery" );
import "../../../styles/jstree.css";
//var jtree = require("jstree");

import "../../ui/jstree/jstree.js";
//import {initJSTreeModule} from "../../ui/jstree/jstree";

export class MCBUtils {
       constructor()
       {
           throw new Error("Cannot create instances of this class");
       }
}

      /*
MCBUtils.shuffleArray = function(array) {

    for (var i = array.length - 1; i > 0; i--)
    {
        var j = Math.floor(Math.random() * (i + 1));
        var temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    let hhh = 0;
}
MCBUtils.generateRandomInteger = function(nIntFr, nIntTo) {

    let nInt = Math.floor(Math.random() * (nIntTo - nIntFr + 1)) + nIntFr;
    return nInt;
}

MCBUtils.generateRandomDate = function(nYearFr, nYearTo) {

    let nD =  Math.floor(Math.random() * 28);
    let nM = Math.floor(Math.random() * 12);
    let nY = Math.floor(Math.random() * (nYearTo - nYearFr + 1)) + nYearFr;

    let nTime = new Date(nY, nM, nD).getTime();
    return nTime;
}

MCBUtils.generateRandomRecentDate = function() {
    let nMonth = 10;//Math.floor(Math.random() * 7);
    let nDay =  1 + Math.floor(Math.random() * 8);
    let nTimeLastTested = new Date(2020, nMonth, nDay).getTime();
    return nTimeLastTested;
}
      */


























/*
MCBUtils.empty = function(){
    let sss = 0;
}
*/








MCBUtils.Test = async function(obPackage) {

    let nColCount = 100;
    let grid = grok.shell.addTableView(grok.data.demo.randomWalk(nRowCount, nColCount)).grid;


    return;

    //let dtt = new DateTime(new Date());
    grok.data.loadTable("http://localhost:8080/api/packages/published/files/GNF/admin/chembl_100k.csv")
        .then(async function(dframe){

            let colStructure = dframe.getCol("smiles");
            const nRecordCount = colStructure.length;
            let arMolFiles = new Array(nRecordCount);
            let strMolFile = null;
            for(var nR=0; nR<nRecordCount; ++nR) {
                strMolFile = colStructure.get(nR);
                arMolFiles[nR] = strMolFile;
            }
           /*
            let service = null;
            try{service = await MolServiceNew.createMolService(MolServiceChemEngine.RDKIt, arMolFiles, false);}
            catch(e)
            {
                let fghfhh= 0;
            }

            let arFlags = new Array(nRecordCount);
            await service.structureSearch(arFlags, arMolFiles[0]);
             */
            let ddd = 0;
        });//t => grok.shell.addTableView(t)
    //let dfs = await grok.functions.eval('OpenServerFile("https://dev.datagrok.ai/files/demo.testjobs.files.demofiles/chem/chembl/chembl_100k.csv")');


    return;
    let SAMPLE_SMILES =[
        "CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3",
        "COc1ccc2cc(ccc2c1)C(C)C(=O)Oc3ccc(C)cc3OC",
        "COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3",
        "CC(C(=O)NCCS)c1cccc(c1)C(=O)c2ccccc2",
        "FC(F)(F)c1ccc(OC2CCNCC2)cc1",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCCc3ccccc3",
        "COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO[N+](=O)[O-]",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO",
        "CN1CCC(CC1)Oc2ccc(cc2)C(F)(F)F",
        "COc1cc(C)ccc1OC(=O)C(C)c2ccc(CC(C)C)cc2",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)OCCCc2cccnc2",
        "COc1ccc(\C=N\NC(=N)N)c(Cl)c1OC",
        "Nc1ncnc2c1c(Br)cn2[C@@H]3OC[C@@H](O)[C@H]3O",
        "CNc1ncnc2c1c(I)cn2[C@@H]3O[C@H](C)[C@@H](O)[C@H]3O",
        "CN1CCC(O)(CC1)c2ccccc2",
        "OC(COc1ccccc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
        "OC(COc1ccc(Cl)cc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
        "OC(COc1ccc(Br)cc1)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
        "COC(=O)c1ccccc1OCC(O)CN2CCC(CC2)Oc3ccc(cc3)C(F)(F)F",
        "CCC1(CC)CC(CCNC(=O)c2ccc(OC)cc2)OC1=O"];






    let nRowCount = 10000;
    let col = DG.Column.fromType(DG.TYPE.STRING, "Mol", nRowCount);
    //col.semType = "Molecule";
    let arMols = new Array(nRowCount);
    let arSmiles = new Array(nRowCount);
    let strSmiles = null;
    var molFrag = rdKitModule.get_mol(SAMPLE_SMILES[0]);
    //molFrag.delete();

    for(var n=0; n<nRowCount; ++n) {
        try {
            strSmiles = SAMPLE_SMILES[n%SAMPLE_SMILES.length];
            col.set(n, strSmiles);
            arSmiles[n] = strSmiles;

            var mol = rdKitModule.get_mol(arSmiles[n]);
            arMols.push(mol);
            //mol.delete();


            //let match = mol.get_substruct_match(molFrag);
            //let ffggfgfg = 0;

        } catch (e) {
            console.log(e.stack);
            let rdfg = 0;
        }
    }

    let nTimeStart = new Date().getTime();

    try {
        grok.chem.similarityScoring(col,  "COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3", {sorted: true}).then(function (dframeResult) {
            let nTimeEnd = new Date().getTime();
            console.log("Spent: " + (nTimeEnd - nTimeStart) / 1000);

            let rtyrtyr = 0;

        }).catch(function(error)
        {
            let err = error.message;
        });

    }catch(ee)
    {
        let gghj = 0;
    }



    return;

    let service = null;
    try{service = await MolService.createMolService(arSmiles, false);}
    catch(e)
    {
        let fghfhh= 0;
    }

    let arFlags = new Array(nRowCount);

    await service.structureSearch(arFlags, strSmiles);
    await service.similaritySearch(arFlags, strSmiles, 0.3);
    await service.tanimotoScores(arFlags, strSmiles);

    let cvbcvb = 0;
    return;


    /*
   let dtframe =  null;
   try{dtframe = DG.DataFrame.fromColumns([col]);}
   catch(e)
   {
       let a = 0;
   }  */

    nTimeStart = new Date().getTime();

    grok.chem.substructureSearch(col,  "COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3", false).then(function (dframeResult) {
        let nTimeEnd = new Date().getTime();
        console.log("Spent: " + (nTimeEnd - nTimeStart)/1000);

        let rtyrtyr = 0;
        // for(var n=0; n<nRecordCount; ++n)
        //{
        //  console.log(dframeResult.get("index", n))
        //}
    }).catch(e)
    {
        let err = e.message;
    };
    return;

    grok.chem.similarityScoring(col,  "COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3", {sorted: true}).then(function (dframeResult) {
        let nTimeEnd = new Date().getTime();
        console.log("Spent: " + (nTimeEnd - nTimeStart)/1000);

        let rtyrtyr = 0;
        // for(var n=0; n<nRecordCount; ++n)
        //{
        //  console.log(dframeResult.get("index", n))
        //}
    });



    return;


     /*
    let nColCount = 100;
    let ggrid = grok.shell.addTableView(grok.data.demo.randomWalk(nRowCount, nColCount)).grid;

    let button = null;
    let btn = ui.button('Click', function(){

        let eContent = ui.splitV([ui.divText("AAAAAAAAAAAAAAAAAAAA"), ui.divText("BBBBBBBBBBBBBBBBBBB")]);
        let popup = ui.showPopup(eContent, button, true);
    });

    button = btn;
    ui.dialog("Popup Bug Demo").add(btn).show();
       */
    return;












    for(var nVCol=0; nVCol<20; ++nVCol)
    {
        let colVirt;
        colVirt = dframee.columns.addNewVirtual("Parent " + nVCol, function (nRow) {
            return nRow;
        });


        GridUtils.setSemType(colVirt, new ActivityConcentrationSemType());
        //colVirt.semType.init(package);

        let arColNames = [];
        for(var nPri=0; nPri<5; ++nPri)
        {
            arColNames[nPri] = "#" + (nVCol*5 + nPri);
        }

        GridUtils.getSemType(colVirt).setPrimitiveColumnNames(arColNames);
    }



    let viewT = grok.shell.addTableView(dframee);
     grid = viewT.grid;


    grok.events.onViewerAdded.subscribe((args) => {

        let viewer = args.args.viewer;
        //let type = viewer.type;
        let root = viewer.root;
        //console.log("Event " + e);
        // alert("View Added: " + view.type);
        //if (type === DG.VIEWER.SCATTER_PLOT) {} else {}
        //this.layoutMain(view);

        viewer.onEvent('d4-column-combo-box-popup-show').subscribe((args) => {

            let viewer = args.args.viewer;
            let opts = viewer.getOptions();

            let eDiv = ui.div();
            eDiv.id = "datajtree";
            viewer.root.appendChild(eDiv);

            //alert("Column Chooser: " + args.args.selectorName);

            let cbb = args.args.comboBox;
            let dial = ui.dialog();

            try {
                MCBUtils.createColumnSelector(eDiv, dframee, function(strColName, bRollover){

                    if(args.args.selectorName === "xColumnName") {
                        viewer.setOptions({x: strColName});
                    }
                    else if(args.args.selectorName === "yColumnName") {
                        viewer.setOptions({y: strColName});
                    }

                    else if(args.args.selectorName === "colorColumnName") {
                        viewer.setOptions({colorColumnName: strColName});
                    }

                    else if(args.args.selectorName === "sizeColumnName") {
                        viewer.setOptions({sizeColumnName: strColName});
                    }

                        //}
                        //else if(viewer.type === DG.VIEWER.HISTOGRAM)
                    //{
                    else if(args.args.selectorName === "valueColumnName") {
                        viewer.setOptions({valueColumnName: strColName});
                    }

                    else if(args.args.selectorName === "splitColumnName") {
                        viewer.setOptions({splitColumnName: strColName});
                    }

                    else if(args.args.selectorName === "categoryColumnName") {
                        viewer.setOptions({splitColumnName: strColName});
                    }


                    if(!bRollover)
                        dial.close();
                });
            }
            catch(e)
            {
                let nnn = 0;
            }


            dial.add(eDiv).show();

            //let ccb = DG.ColumnComboBox.create(grok.data.demo.demog(), (c) => c.type === 'int');
            //let v = grok.shell.newView('demo: column combo box', [ccb]);

            args.preventDefault();
        });

    });


    let typee = null;
    try{typee = grid.type;}
    catch(e)
    {
        let msg = e.message;
    }
    /*
 ui.dialog()
     .add(ui.button('Scroll', function(){
         grid.scrollToPixels(250, 250);
         let nXScroll = grid.horzScroll.min;//has correct value 150
         let nYScroll = grid.vertScroll.min;//has wrong value 5
         })
     ).show();
      */
    return;
    /*
 let bClicked = false;


 let nRowCount = 1000;
 let nColCount = 100;
 let grid = grok.shell.addTableView(grok.data.demo.randomWalk(nRowCount, nColCount)).grid;

 let rowHeader = grid.columns.rowHeader;

 rxjs.fromEvent( grid.overlay, 'click').subscribe((e) => {

     let nButton = e.button;
     if (nButton === 0)
     {

         grid.scrollToPixels(150, 150);
         let nXScroll = grid.horzScroll.min;
         let nYScroll = grid.vertScroll.min;
         let ddd = 0;
    }
 });
      */


    return;


    function processErr(e)
    {
        let dfg = 0;
    }

    function processResults(dframe)
    {
        let b = false;
    }


    function simisearch(dframe)
    {
        let col = dframe.col("Cpd");
        let mol = null;
        let arIndex = null;
        let search = null;
        try{search = new OCLM.SSSearcherWithIndex();}
        catch(e)
        {
            let dd = 0;
        }

        let strSmiles = null;
        let nRowCount = dframe.rowCount;
        for(var n=0; n<nRowCount; ++n)
        {
            strSmiles = col.get(n);//"CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3";
            strSmiles = strSmiles.getSmiles();
            mol = OCL.Molecule.fromSmiles(strSmiles);
            arIndex = search.createIndexMy(mol);

            let v = 0;
        }


        let fffu = 0;



        //grok.shell.addTableView(dframe);
        //let p = grok.chem.similaritySearch(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', DG.SIMILARITY_METRIC.TANIMOTO, 1000, 0.3);
        //p = p.then(processResults);
    }


    //let dframe = MCBUtils.generateRandomCpdDataFrame(1000, package);
    simisearch(dframe);
    //grok.shell.addTableView(dframe);




    return;

    // let p = grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv');
    //p = p.then(simisearch);
    /*
p.then(molecules => grok.chem.similaritySearch(molecules.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1')
 .then(similar => grok.shell.addTableView(similar)));
      */
    let bClicked6 = false;
    return;




    /*
        rxjs.fromEvent( grid.overlay, 'click').subscribe((e) => {

            let nButton = e.button;
            if (nButton === 0)
            {
                let cell = grid.hitTest(e.offsetX, e.offsetY);
                if (cell.isColHeader) {
                    bClicked = true;

                    let colGrid = null;
                    for (var nCol = 0; nCol < nColCount; ++nCol) {
                        colGrid = grid.columns.byIndex(nCol);
                        colGrid.width = 75 + nCol;
                    }

                }
            }
        });

        grid.onCellRender.subscribe(function (args) {
              if(!bClicked)
                return;

            let nCol = args.cell.gridColumn.idx;

            if(args.cell.isColHeader) {
                console.log("Column Header " + nCol);
            }
        });


          return;


                   */
    /*
grid.onCellRender.subscribe(function (args) {

if(!bClicked)
  return;

let nRow = args.cell.gridRow;
let nCol = args.cell.gridColumn.idx;

if(nRow < 0 && nCol === 0)//my changes !args.cell.isTableCell && !args.cell.isRowHeader && !args.cell.isRowHeader)
{
 console.log("Top Left " + nRow + " " + nCol + " " + args.cell.isColHeader + " " + args.cell.isRowHeader + " " + args.cell.isTableCell)

}
else if(args.cell.isColHeader)//nRow < 0 && nCol > 0)
{
  console.log("Col Header " + nRow + " " + nCol + " " + args.cell.isColHeader + " " + args.cell.isRowHeader + " " + args.cell.isTableCell);
}
else if(args.cell.isRowHeader)//args.cell.cellType === "row header")
{
  console.log("Row Header " + nRow + " " + nCol + " " + args.cell.isColHeader + " " + args.cell.isRowHeader + " " + args.cell.isTableCell);
}
else if(args.cell.isTableCell)
{

}
//args.preventDefault();
});         */


    /*
    rxjs.fromEvent( grid.overlay, 'click').subscribe((e) => {

        let nButton = e.button;
        if (nButton === 0)
        {
            let cell = grid.hitTest(e.offsetX, e.offsetY);

                let bRH = cell.isRowHeader; //logically should return true, but returns false
                let bCH = cell.isColHeader; //logically should return true, but returns false

                let colGrid = cell.gridColumn;

                let nIdx = 0;
                try {nIdx = colGrid.idx;}
                catch (e) {
                    let strErr = e.message; //shouldn't happen in any case

                }
        }
    });*/



    return;

    //let colGrid = null;
    try {
        for (var n = 30; n < 50; ++n) {
            //colGrid = grid.columns.byIndex(n);
            //colGrid.visible = false;
            //colGrid.width = 165;
        }
    }
    catch(e)
    {
        let ddd = 0;
    }

    rxjs.fromEvent( grid.overlay, 'click').subscribe((e) => {

        let nButton = e.button;
        if (nButton === 0)
        {
            let cell = grid.hitTest(e.offsetX, e.offsetY);
            if (cell.isColHeader) {

                //grid.scrollToPixels(150, 170);
                grid.scrollToCell(1, 20);
            }
        }
    });


}

MCBUtils.Test1 = function() {

    let dframe = grok.data.demo.demog(1000);
    let order = new Int32Array(dframe.rowCount);
    for (let i = 0; i < dframe.rowCount; i++) {

        order[i] = dframe.rowCount - i - 1;
    }

    let bClicked = false;

    let view = grok.shell.addTableView(dframe);
    view.grid.onCellRender.subscribe(function (args) {

        if(args.cell.isTableCell) {
            if(!bClicked)
                return;

            let nCol = args.cell.gridColumn.idx;
            let nRow = args.cell.gridRow;
            console.log(nCol + " " + nRow);
        }
    });


    rxjs.fromEvent( view.grid.overlay, 'click').subscribe((e) => {
        let cell = view.grid.hitTest(e.offsetX, e.offsetY);

        let nButton = e.button;
        if (nButton === 0)
        {
            let cell = view.grid.hitTest(e.offsetX, e.offsetY);
            if (cell.isColHeader) {
                bClicked = true;
                view.grid.setRowOrder(order);
            }
        }
    });



    //alert("before");

    //alert("after");

    let aaa = 0;
    /*
    let t = DG.DataFrame.fromColumns([
        DG.Column.fromList('int', 'int', [1, 2, 3,8,9]),
        DG.Column.fromList('double', 'double', [1.1, 2.1, 3.1, ,8.9]),
        DG.Column.fromList('string', 'string', ["a", "b", "c", "d","t"])
    ]);

    let view = grok.shell.addTableView(t);

    view.grid.setOptions({
        colHeaderHeight: 75,
        rowHeight: 150
    });

    //view.grid.columns.byIndex(0).width = 200;
    view.grid.columns.rowHeader.width = 200;
    */
}



MCBUtils.createColumnSelectorData = function(dframe, strPattern) {
    const lstCols = dframe.columns.toList();
    const nColCount = lstCols.length;
    let col = null;
    let colPrim = null;
    let typeSem = null;
    let node = null;
    let nChildCount = -1;
    let strParamName = null;
    let arChildren = null;
    let arPrimColNames = null;
    let strIconURL = null;
    const arNodes = [];

    for (var nCol = 0; nCol < nColCount; ++nCol) {
        col = lstCols[nCol];
        typeSem = GridUtils.getSemType(col);
        if (!(typeSem instanceof SemType) || !DGApp.isVirtual(col))
            continue;

        if(strPattern !== undefined && !col.name.toLowerCase().includes(strPattern.toLowerCase()))
            continue;

        node = {};
        node.text = col.name;
        node.id = "symtype" + nCol;

        if(typeSem.getIconURL() !== null)
            node.icon = typeSem.getIconURL();
        //node.children = true;


        arPrimColNames = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);

        nChildCount = arPrimColNames === null ? 0 : arPrimColNames.length;
        arChildren = new Array(nChildCount);
        for(var nChild=0; nChild<nChildCount; ++nChild)
          {
                colPrim = dframe.columns.byName(arPrimColNames[nChild]);
                typeSem = GridUtils.getSemType(colPrim);
                strIconURL = typeSem.getIconURL();
                if(strIconURL === null)
                    strIconURL = "jstree-file";

                strParamName = SemType.getChildColHeaderValue(colPrim);//)arPrimColNames[nChild];
                arChildren[nChild] = {text: strParamName, icon: strIconURL, data: arPrimColNames[nChild]};
          }
            node.children = arChildren;

        arNodes.push(node);
    }

    return arNodes;
}


MCBUtils.updateColumnSelector = function(eDiv, data) {
    let idDiv = eDiv.id;

    let str = "#" + idDiv;
    let ejtree = $(str);
    let jtree = ejtree.jstree(true);
    jtree.settings.core.data = data;
    jtree.refresh();
}

MCBUtils.createColumnSelector = function(eDiv, dframe, fnOnSelect) {

    const arNodes = MCBUtils.createColumnSelectorData(dframe);
    const idDiv = eDiv.id;
    const str = "#" + idDiv;
    const jtree = $(str).jstree({
        'core': {
            'data': arNodes,
            'multiple': false
        }
    });

    jtree.on('select_node.jstree', function(e, data) {
        //const strColName = data.node.text;

        if(fnOnSelect !== undefined && data.node.parent !== "#") {
            const strColName = data.node.data;
            fnOnSelect(strColName, false);
        }

    }).on('hover_node.jstree', function(e, data) {

        if(fnOnSelect !== undefined && data.node.parent !== "#") {
            const strColName = data.node.data;
            fnOnSelect(strColName, true);
        }

    })

    return arNodes;

    //$('#datajtree')
    /*
    $("'#" + datajtree + "'").jstree({
        'core': {
            'data': function (node, cb) {
                if (node.id === "#") {

                    let nodeOne = {};
                    nodeOne.text = "Root 1";
                    nodeOne.id = 1;
                    nodeOne.children = true;

                    let nodeTwo = {};
                    nodeTwo.text = "Root 2";
                    nodeTwo.id = 2;
                    nodeTwo.children = true;


                    cb([nodeOne, nodeTwo]);
                } else {
                    cb(["Child 1 for Node " + node.id, "Child 2 for Node " + node.id]);
                }
            }
        }
    });*/

}

