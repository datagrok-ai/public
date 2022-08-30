import {ClassMap} from "../../lang/ClassMap";
import {RDKitCpdEntityRenderer} from "../gnf/cpd/ui/rdkit/RDKitCpdEntityRenderer";
import {SemType} from "../SemType";
import {DefaultSemEntityRenderer} from "../ui/DefaultSemEntityRenderer";
import {CpdSemType} from "../gnf/cpd/CpdSemType";
import {InVivoSemType} from "../gnf/invivo/InVivoSemType";
import {InVivoEntityRenderer} from "../gnf/invivo/ui/InVivoEntityRenderer";
import {InSilicoSemType} from "../gnf/insilico/InSilicoSemType";
import {InSilicoEntityRenderer} from "../gnf/insilico/ui/InSilicoEntityRenderer";
import {IC50SemType} from "../gnf/ic50/IC50SemType";
import {IC50EntityRenderer} from "../gnf/ic50/ui/IC50EntityRenderer";
import {FlowTextRenderer} from "../ui/FlowTextRenderer";
import {EmptyButtonGridColumnRenderer} from "../ui/EmptyButtonGridColumnRenderer";
import {GridUtils} from "../../utils/GridUtils";

export class SemEntityAppFramework
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}

class SemEntityAppFrameworkImpl
{

}

SemEntityAppFrameworkImpl.onGridAdded = async function (viewerGrid, obPackage)
{
    let mapGridCellRenderers = new ClassMap();

    let rendererCpd = await RDKitCpdEntityRenderer.create();
    mapGridCellRenderers.set(SemType, new DefaultSemEntityRenderer());
    mapGridCellRenderers.set(CpdSemType, rendererCpd);
    mapGridCellRenderers.set(InVivoSemType, new InVivoEntityRenderer());
    mapGridCellRenderers.set(InSilicoSemType, new InSilicoEntityRenderer());
    mapGridCellRenderers.set(IC50SemType, new IC50EntityRenderer());

    let mapGridHeaderRenderers = new ClassMap();
    mapGridHeaderRenderers.set(SemType, new FlowTextRenderer());
    mapGridHeaderRenderers.set(CpdSemType, new EmptyButtonGridColumnRenderer());

    let dframe = viewerGrid.dataFrame;

    viewerGrid.dart.m_ctxSort = null;

    let strURLWR = obPackage.webRoot;
    console.log("url: " + strURLWR);

    let parent = viewerGrid.root.parentElement;
     /*
    var img = null;
    try {img = GridUtils.loadImage(strURLWR + "/images/adjust_rows_height.png");}
    catch (e) {
        let err0 = 0;
    }

    let btnH = ui.button("H", function() {
        GridUtils.adjustRowSizeToColumnWidth(viewerGrid);
    }, "Adjust Rows Height");

    if(img !== null)
        btnH.innerHTML = img.outerHTML;

    try {img = GridUtils.loadImage(strURLWR + "dist/images/fit_zoom.png");}
    catch (e) {}


    let btnF = ui.button("F", function() {
        GridUtils.autoFitColumns(viewerGrid, mapGridCellRenderers);
    }, "Zoom to Fit");

    if(img !== null)
        btnF.innerHTML = img.outerHTML;


    try {img = GridUtils.loadImage(strURLWR + "dist/images/default_zoom.png");}
    catch (e) {}

    let btnD = ui.button("D", function() {
        GridUtils.setDefaultColumnSize(viewerGrid, mapGridCellRenderers);
    }, "Default Zoom");

    if(img !== null)
        btnD.innerHTML = img.outerHTML;
*/


    let split = null;

    try{split = ui.splitH([btnH, btnF, btnD]);}
    catch(e)
    {
        let rrrrr = 0;
    }

    parent.appendChild(split);


    const nInsets = 2;
    const nDefaultColWidth = 75;
    const nMaxHeaderHeight = 100;


    //let eCanvas =  viewerGrid.canvas;
    var offscreen = new OffscreenCanvas(256, 256);
    var ctx = offscreen.getContext('2d');


   // var ctx = eCanvas.getContext("2d");
    var arLines = new Array();
    ctx.font = "13px Arial";


    let nHPrefHeader = GridUtils.calcPrefColHeaderHeight(ctx, dframe, nDefaultColWidth - (nInsets << 1), nInsets);
    if(nHPrefHeader > nMaxHeaderHeight)
        nHPrefHeader = nMaxHeaderHeight;


    viewerGrid.setOptions({
        colHeaderHeight: nHPrefHeader,
        //rowHeight: 150
    });


    //my changes viewerGrid.columns.rowHeader.width = 150;
    //viewerGrid.columns.byIndex(0).width = 250;


    let nColCpd = GridUtils.findGridColumnBySemType(viewerGrid, CpdSemType);
    let colGcpd = viewerGrid.columns.byIndex(nColCpd);  //my changes 1
    colGcpd.visible = false;

    GridUtils.initZoom(viewerGrid);
    GridUtils.initDefaultCellsSizes(viewerGrid, false, mapGridCellRenderers);


    let nTopLeftCount = 0;
    let nColHeaderCount = 0;
    let nRowHeaderCount = 0;

    viewerGrid.onCellRender.subscribe(function (args) {

        let nWCell = args.bounds.width;
        let nHCell = args.bounds.height;

        //let nCol = args.cell.gridColumn.idx;
        let nRow = args.cell.gridRow;
        let nCol = args.cell.gridColumn.idx;
    });

    }

SemEntityAppFramework.Init = async function(obPackage)
{
    grok.events.onViewerAdded.subscribe((args) => {
        let viewer = args.args.viewer;
        let type = viewer.type;
        if(type === DG.VIEWER.GRID)
        {
            viewer.o

            SemEntityAppFrameworkImpl.onGridAdded(viewer, obPackage);
        }

        let rrrrrr = 0;
    });

}
