/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

import {PinnedColumn} from "@datagrok-libraries/gridext/src/PinnedColumn";

import {GNFUtils} from "./entity/gnf/GNFUtils";
import {DGApp} from "./entity/app/DGApp";
import {ChemUtils} from "./chem/util/ChemUtils";
import {IntuenceDiscovery} from "./intuence/IntuenceDiscovery";
import {CpdSemType} from "./entity/gnf/cpd/CpdSemType";
import {CDFComposites} from "./service/nx/cdf/composite/CDFComposites";
import {AnalysisLoader} from "./service/nx/analysis/AnalysisLoader";
import {ActivityConcentrationSemType} from "./entity/nx/activity/ActivityConcentrationSemType";
import {ActivityPercentSemType} from "./entity/nx/activity/ActivityPercentSemType";
import {SemTypeConfig} from "./intuence/SemTypeConfig";
import {QNumDateSemType} from "./entity/qnum/QNumDateSemType";
import {TextDateSemType} from "./entity/text/TextDateSemType";
import {IdeaSemType} from "./entity/nx/idea/IdeaSemType";
import {ObjectSemType} from "./entity/object/ObjectSemType";
import {ProgressDialog} from "./ui/progress/ProgressDialog";
import {MolService} from "./chem/MolService";
import {MolServiceNew} from "./chem/MolServiceNew";
import {RDKitEngine} from "./chem/rdkit/RDKitEngine";
import {RDKitEngineNew} from "./chem/rdkit/RDKitEngineNew";
import {RDKitCpdGridRowHeaderRenderer} from "./entity/gnf/cpd/ui/rdkit/RDKitCpdGridRowHeaderRenderer";
import {RDKitCpdGridCellRenderer} from "./entity/gnf/cpd/ui/rdkit/RDKitCpdGridCellRenderer";
import {SemType} from "./entity/SemType";
import {ObjectGridCellRenderer} from "./entity/object/ui/ObjectGridCellRenderer";
import {InVivoSemType} from "./entity/gnf/invivo/InVivoSemType";
import {InVivoGridCellRenderer} from "./entity/gnf/invivo/ui/InVivoGridCellRenderer";
import {InSilicoSemType} from "./entity/gnf/insilico/InSilicoSemType";
import {InSilicoGridCellRenderer} from "./entity/gnf/insilico/ui/InSilicoGridCellRenderer";
import {IC50SemType} from "./entity/gnf/ic50/IC50SemType";
import {IC50GridCellRenderer} from "./entity/gnf/ic50/ui/IC50GridCellRenderer";
import {IdeaGridCellRenderer} from "./entity/nx/idea/ui/IdeaGridCellRenderer";
import {ActivityConcentrationGridCellRenderer} from "./entity/nx/activity/ui/ActivityConcentrationGridCellRenderer";
import {QNumDateGridCellRenderer} from "./entity/qnum/ui/QNumDateGridCellRenderer";
import {TextDateGridCellRenderer} from "./entity/text/ui/TextDateGridCellRenderer";
import {ButtonGridColumnHeaderRenderer} from "./entity/ui/ButtonGridColumnHeaderRenderer";
import {NullSemType} from "./entity/NullSemType";
import {EmptyButtonGridColumnRenderer} from "./entity/ui/EmptyButtonGridColumnRenderer";
import {GridRenderersConfig} from "./intuence/GridRenderersConfig";
import {CpdColumnHeaderRenderer} from "./entity/gnf/cpd/ui/CpdColumnHeaderRenderer";
import {GridUtils} from "./utils/GridUtils";
import {PrimitiveSemType} from "./entity/primitive/PrimitiveSemType";
import {PrimitiveGridCellRenderer} from "./entity/primitive/ui/PrimitiveGridCellRenderer";
import {CpdPrimitiveSemType} from "./entity/gnf/cpd/CpdPrimitiveSemType";
import {CpdPrimitiveGridCellRenderer} from "./entity/gnf/cpd/ui/CpdPrimitiveGridCellRenderer";
import {ButtonGridRowHeaderRenderer} from "./entity/ui/ButtonGridRowHeaderRenderer";
import {RangeSliderFilter} from "./ui/filters/RangeSliderFilter";
import {IC50Filter} from "./entity/gnf/ic50/ui/IC50Filter";
import {CheckBoxesFilter} from "./ui/filters/CheckBoxesFilter";
import {InSilicoFilter} from "./entity/gnf/insilico/ui/InSilicoFilter";
import {InVivoFilter} from "./entity/gnf/invivo/ui/InVivoFilter";
import {ComboBoxFilter} from "./ui/filters/ComboBoxFilter";
import {RecentTimeFilter} from "./ui/filters/RecentTimeFilter";
import {RandUtils} from "./utils/RandUtils";
import {MultilineTextFilter} from "./ui/filters/MultilineTextFilter";
import {CpdFilter} from "./entity/gnf/cpd/ui/CpdFilter";
import {ActivityConcentrationFilter} from "./entity/nx/activity/ui/ActivityConcentrationFilter";
import {QNumDateFilter} from "./entity/qnum/ui/QNumDateFilter";
import {TextDateFilter} from "./entity/text/ui/TextDateFilter";
import {ObjectFilter} from "./entity/object/ui/ObjectFilter";
import {IdeaFilter} from "./entity/nx/idea/ui/IdeaFilter";
import {FiltersUtils} from "./ui/filters/FiltersUtils";
import {RDKitDGGridRowHeaderRenderer} from "./entity/gnf/cpd/ui/rdkit/RDKitDGGridRowHeaderRenderer";
import {RDKitDGCpdGridCellRenderer} from "./entity/gnf/cpd/ui/rdkit/RDKitDGCpdGridCellRenderer";
import {MetaDataGridColHeaderRenderer} from "./intuence/MetaDataGridColHeaderRenderer";
import {ProgressIndicatorWrapper} from "./ui/progress/ProgressIndicatorWrapper";
import {FlowTileLayoutManager, TableViewLayoutManager} from "./ui/layout/TableViewLayoutManager";
import {DGUtils} from "./utils/DGUtils";
import {TextUtils} from "./utils/TextUtils";
import {TaskStatus} from "./concurrent/TaskStatus";
import { Observable } from 'rxjs';
import {CompoundFilter} from "./entity/nx/cpd/ui/CompoundFilter";
import {GridRowHeaderNew} from "./ui/grid/GridRowHeaderNew";
import * as PC from "./ui/grid/GridRowHeaderNew";
//import {tsconfig.json} from "./ui/grid/PinnedColumn";

export let _package = new DG.Package();

const _qnumBuf = new DataView(new ArrayBuffer(8));
DG.Qnum.getValue = function(x)
{
    _qnumBuf.setFloat64(0, x);
    let last =  _qnumBuf.getInt8(7) & 0xFC;
    _qnumBuf.setInt8(7, last);
    return  _qnumBuf.getFloat64(0);
};

//name: addFrozenColumn
//description: Adds a frozen column
//input: object colGrid
//output: object result
export function addFrozenColumn(colGrid) {
    let headerRows = null;
    try{headerRows = new GridRowHeaderNew(colGrid);}
    catch(e) {
        alert(e.message);
        throw e;
    }
    return headerRows;
}

//name: compoundFilter
//description: A filter to search by multiple IDs
//tagssss: filter
//output: filter result
export function compoundFilter() {
    return new CompoundFilter();
}

//name: getCompoundFilter
//description: Reurns the compound filter from a Filters viewer.
//input: object viewerFilters
//output: object result
export function getCompoundFilter(viewerFilters) {
    const dart = viewerFilters.toDart();
    let filter = dart.m_filterCpd;
    return filter === undefined ? null : filter;
}
//name: addCompoundFilter
//description: Adds a compound filter to a Filters viewer.
//input: object viewerFilters
//output: object result
export function addCompoundFilter(viewerFilters) {

    const dart = viewerFilters.toDart();
    let filter = dart.m_filterCpd;
    if(filter !== null && filter !== undefined)
        return filter;

    filter = compoundFilter();
    const drame = viewerFilters.dataFrame;
    filter.attach(drame);
    FiltersUtils.addFilter(filter, viewerFilters);

    dart.m_filterCpd = filter;

    return filter;
}


//name: removeCompoundFilter
//description: Removes a compound filter from a Filters viewer.
//input: object viewerFilters
//output: object result
export function removeCompoundFilter(viewerFilters) {

    const dart = viewerFilters.toDart();
    const filter = dart.m_filterCpd;
    if(filter === undefined || filter === null)
        return null;

    filter.detach();
    const b = FiltersUtils.removeFilter(filter, viewerFilters);
    if(b) {
        dart.m_filterCpd = null;
        return filter;
    }

    return null;
}


//name: Metadata OnDemand
//tags: viewerr
//output: viewer result
/*
export function metadataViewer()
{
    return new GridRowHeaderEmpty();
}*/

//name: Entity Filter
//description: Entity filters panel.
//tags: filterr
//output: filter result
/*
export function entityFilter()
{
    return new FilterPanel([RangeSliderFilter, IC50Filter]);
} */


//name: handshake
//description: Handshake Responder
//output: bool result
export function handshake() {
    return true;
}


//tags: app
//name: GNF MedChem Browser 5 Evaluation
export async function startApp() {

    const nRecordCount = 1000;
    const dframe = await GNFUtils.generateRandomDataFrame(nRecordCount, _package);

    const configRenderers = new GridRenderersConfig();
    //const rendererRowHeader= await RDKitCpdGridRowHeaderRenderer.create();
    //configRenderers.setRowHeaderRenderer(rendererRowHeader);
    //const rendererRowHeader= new OCLCpdGridRowHeaderRenderer();

    const rendererCpdCell = await RDKitCpdGridCellRenderer.create();
    //const rendererCpdCell = new OCLCpdGridCellRenderer();

    configRenderers.setCellRenderer(SemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(ObjectSemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(CpdSemType, rendererCpdCell);
    configRenderers.setCellRenderer(InVivoSemType, new InVivoGridCellRenderer());
    configRenderers.setCellRenderer(InSilicoSemType, new InSilicoGridCellRenderer());
    configRenderers.setCellRenderer(IC50SemType, new IC50GridCellRenderer());


    //const mapGridHeaderRenderers = new ClassMap();
    configRenderers.setColHeaderRenderer(SemType, new MetaDataGridColHeaderRenderer);//new ButtonGridColumnHeaderRenderer());
    //configRenderers.setColHeaderRenderer(CpdSemType, new CpdColumnHeaderRenderer());
    configRenderers.setColHeaderRenderer(NullSemType, new EmptyButtonGridColumnRenderer());

    configRenderers.setDeffaultColHeaderRenderer(new ButtonGridColumnHeaderRenderer());

    configRenderers.setColRowHeaderRenderer(new CpdColumnHeaderRenderer());


    const arFilters = [RangeSliderFilter, CpdFilter, InVivoFilter, InSilicoFilter, IC50Filter];

    const viewSS = grok.shell.addTableView(dframe);
    const grid = viewSS.grid;
    await DGApp.open(_package, viewSS, configRenderers, arFilters);

    const nColCpd = GridUtils.findGridColumnBySemType(grid, CpdSemType);
    const colGcpd = grid.columns.byIndex(nColCpd);  //my changes 1
    //my changes colGcpd.visible = false;


    //my changes FiltersUtils.openFilter(colGcpd.column);
}


//tags: app
//name: MedChem Table
export async function startIntuence() {

    window.devicePixelRatio = 1;

    const configSemTypes = new SemTypeConfig();
    //configSemTypes.add(CpdSemType);
    configSemTypes.add(QNumDateSemType);
    configSemTypes.add(TextDateSemType);
    configSemTypes.add(ActivityConcentrationSemType);
    configSemTypes.add(ActivityPercentSemType);
    configSemTypes.add(IdeaSemType);
    configSemTypes.add(ObjectSemType);

    const arFilterClasses = [CpdFilter, ActivityConcentrationFilter, IdeaFilter, QNumDateFilter, TextDateFilter, RangeSliderFilter, RecentTimeFilter, ComboBoxFilter, CheckBoxesFilter,
        MultilineTextFilter, ObjectFilter];

    const configRenderers = new GridRenderersConfig();

    //const mapGridCellRenderers = new ClassMap();
    //const rendererRowHeader= await RDKitDGGridRowHeaderRenderer.create();//RDKitCpdGridRowHeaderRenderer.create();
    // configRenderers.setRowHeaderRenderer(rendererRowHeader);
    //const rendererRowHeader= new OCLCpdGridRowHeaderRenderer();

    const rendererCpdCell = await RDKitDGCpdGridCellRenderer.create();
    //const rendererCpdCell = new OCLCpdGridCellRenderer();

    configRenderers.setCellRenderer(SemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(ObjectSemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(CpdSemType, rendererCpdCell);
    configRenderers.setCellRenderer(InVivoSemType, new InVivoGridCellRenderer());
    configRenderers.setCellRenderer(InSilicoSemType, new InSilicoGridCellRenderer());
    configRenderers.setCellRenderer(IC50SemType, new IC50GridCellRenderer());
    configRenderers.setCellRenderer(IdeaSemType, new IdeaGridCellRenderer());
    configRenderers.setCellRenderer(ActivityConcentrationSemType, new ActivityConcentrationGridCellRenderer());
    //configRenderers.setCellRenderer(ActivityConcentrationSemType, new ActivityConcentrationGridCellRenderer());
    configRenderers.setCellRenderer(QNumDateSemType, new QNumDateGridCellRenderer());
    configRenderers.setCellRenderer(TextDateSemType, new TextDateGridCellRenderer());
    configRenderers.setCellRenderer(CpdPrimitiveSemType, new CpdPrimitiveGridCellRenderer());

    let rendererDGRDKit = null;
    try{rendererDGRDKit = await grok.functions.call("chem:rdkitCellRenderer");}
    catch(e)
    {
        rendererDGRDKit = null;
    }
    configRenderers.setCellRenderer(PrimitiveSemType, new PrimitiveGridCellRenderer(rendererDGRDKit));

    //const mapGridHeaderRenderers = new ClassMap();
    configRenderers.setColHeaderRenderer(SemType, new MetaDataGridColHeaderRenderer());//new ButtonGridColumnHeaderRenderer());
    //configRenderers.setColHeaderRenderer(CpdSemType, new CpdColumnHeaderRenderer());
    configRenderers.setColHeaderRenderer(NullSemType, new EmptyButtonGridColumnRenderer());

    configRenderers.setDeffaultColHeaderRenderer(new ButtonGridColumnHeaderRenderer());
    //configRenderers.setColRowHeaderRenderer(new CpdColumnHeaderRenderer());

    const dialProgress = new ProgressDialog("Loading Analysis...", true);

    let arFrames = null;
    arFrames = await AnalysisLoader.loadAnalysisGUI(["f9cb0333-9dc8-4baa-bb78-ebf659c716c","fb59941a-cfa1-4cd6-b155-4350c37e29c8",
        "8e4a7844-68a6-47dc-9fdb-0384d100c0b3", "5fca74f5-55a3-4094-9ee6-6cc7bd236627", "18d9031e-5c09-4ca2-a537-4d665a7514e4",
        "335e1176-f093-4a67-bc7b-b0e091ea6996", "0b7740bc-c03e-4335-a0a2-e8329a74de70", "af93950a-d37e-445c-ad9a-41dccfb202de",
        "0e492156-6349-4b6d-9971-838a8fb3fb77", "31611fee-7ea1-4f05-8964-7837adb56a4a"], dialProgress);
     //arFrames = await AnalysisLoader.loadAnalysis("fb59941a-cfa1-4cd6-b155-4350c37e29c8");
    if(arFrames === null)//Cancel was clicked
    {
        await dialProgress.close();
        return;
    }

    const dframe = arFrames[0];

    const viewSS = grok.shell.addTableView(dframe);
    /*
    if(dframe.grid.d === undefined)
    {
        dframe.grid.d = {};
    }*/


    const bHasMetaData = DGApp.hasVirtualColumns(dframe);
    if(!bHasMetaData) {

        try {
            await IntuenceDiscovery.processAnalysisData(dframe, configSemTypes, _package, dialProgress);
        } catch (e) {
            await dialProgress.setUpperText("ERROR")
            let c = 0;
        }
    }
    await IntuenceDiscovery.open(viewSS, configRenderers, arFilterClasses, _package);
    await dialProgress.close();
     //FiltersUtils.openFilter(colCpd);
}



 /*
async function countToTen(s, progress) {
    for (let i = 0; i < 10; i++) {
        if (progress != null) {
            progress.update(i * 10, `${i * 10} % loaded Details does someting...`);
            progress.log(`${i * 10} % loaded Details does someting 1...`);
        }
        await new Promise(resolve => setTimeout(resolve, 1000));
    }
    return `result: ${s}`;
}  */


let m_handlerViewAddedHideVirtual = null;
let m_handlerViewLayoutHideVirtual = null;


//name: getMolService
//description: Returns molecular service for a specified DataFrame
//input: object dframe
//output: object result
export async function getMolService(dframe) {

    const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
    if(nColCpd < 0)
     return null;

    const colCpd = dframe.columns.byIndex(nColCpd);
    const typeSem = GridUtils.getSemType(colCpd);

    const service = typeSem.getMolService();
    return service;
}


//name: initEntityCentricData
//description: initializes virtual columns for the specified data frame
//input: object arDataFrames
//input: object progress
//output: bool result
export async function initEntityCentricData(arDataFrames, progress) {

    //alert("Commoncer Init: " + progress.constructor);
    const dframe = arDataFrames[0];

    //progress.update(0, "0% loaded Starting Meta Data.");
    //return true;


    const bHasMetaData = DGApp.hasVirtualColumns(dframe);

    if(!bHasMetaData) {

        const configSemTypes = new SemTypeConfig();
        //configSemTypes.add(CpdSemType);
        configSemTypes.add(QNumDateSemType);
        configSemTypes.add(TextDateSemType);
        configSemTypes.add(ActivityConcentrationSemType, "../images/dose_response.png");
        configSemTypes.add(IdeaSemType);
        configSemTypes.add(ObjectSemType);

        progress.update(0, "0% loaded Starting Meta Data.");
        const arDescriptors = await AnalysisLoader.populateMetaData(dframe);
        progress.update(99, "100& loaded. Meta Data Loaded.");
        progress.update(0, "0& loaded. Creating Virtual Columns");
        await IntuenceDiscovery.processAnalysisData(dframe, configSemTypes, _package, new ProgressIndicatorWrapper(progress));
/* moved to Id
        if(m_handlerViewAddedHideVirtual === null) {
            m_handlerViewAddedHideVirtual = grok.events.onViewAdded.subscribe((view) => {
                if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
                    DGApp.setVistualColumnsVisible(view.grid, false);
                }
            });
        }

        if(m_handlerViewLayoutHideVirtual === null) {
            m_handlerViewLayoutHideVirtual = grok.events.onViewLayoutApplied.subscribe((layout) => {
                    const view = layout.view;
                    const view1 = view.c$;

                    if (view1.grid !== undefined)
                        DGApp.setVistualColumnsVisible(view1.grid, false);
                }
            );
        }*/
    }

    return true;
}


//name: toEntityCentricGrid
//description: Loads CDF meta data and opens MedChem Grid
//input: object view
//input: object fnFilterOpener
//input: object fnFilterCloser
//output: bool result
export async function toEntityCentricGrid(view, fnFilterOpener, fnFilterCloser) {
   //alert("Start Loading: " + fnFilterOpener);
    grok.events.onViewerAdded.subscribe((args) => {

        const viewer = args.args.viewer;

        console.log("Here");
        //DGApp.setVistualColumnsVisible(view.grid, true);

    });
    if(fnFilterOpener !== null && fnFilterOpener !== undefined)
        FiltersUtils.setFilterOpener(fnFilterOpener);

    if(fnFilterCloser !== null && fnFilterCloser !== undefined)
        FiltersUtils.setFilterCloser(fnFilterCloser);

    const strViewName = view.name;

    const dframe = view.grid.dataFrame;
    const configSemTypes = new SemTypeConfig();
    //configSemTypes.add(CpdSemType);
    configSemTypes.add(QNumDateSemType);
    configSemTypes.add(TextDateSemType);
    configSemTypes.add(ActivityConcentrationSemType, "../images/dose_response.png");
    configSemTypes.add(IdeaSemType);
    configSemTypes.add(ObjectSemType);

    const arFilterClasses = [CpdFilter, ActivityConcentrationFilter, IdeaFilter, QNumDateFilter, TextDateFilter, RangeSliderFilter,
        RecentTimeFilter, ComboBoxFilter, CheckBoxesFilter, MultilineTextFilter, ObjectFilter];

    /*
    const arFilterClasses = [ComboBoxFilter, CheckBoxesFilter, RecentTimeFilter, RangeSliderFilter,
        MultilineTextFilter, CpdFilter, IdeaFilter, ActivityConcentrationFilter, QNumDateFilter, TextDateFilter, ObjectFilter];
*/
    //const arFilterClasses = [ComboBoxFilter, CheckBoxesFilter, RecentTimeFilter, RangeSliderFilter, MultilineTextFilter, CpdFilter, IdeaFilter, ActivityConcentrationFilter, QNumDateFilter, TextDateFilter, ObjectFilter];

    const configRenderers = new GridRenderersConfig();

    //const mapGridCellRenderers = new ClassMap();

    //const rendererRowHeader= await RDKitDGGridRowHeaderRenderer.create();//await RDKitCpdGridRowHeaderRenderer.create();
    //configRenderers.setRowHeaderRenderer(rendererRowHeader);
    //const rendererRowHeader= new OCLCpdGridRowHeaderRenderer();

    const rendererCpdCell = await RDKitDGCpdGridCellRenderer.create();//await RDKitCpdGridCellRenderer.create();
    //const rendererCpdCell = new OCLCpdGridCellRenderer();

    configRenderers.setCellRenderer(SemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(ObjectSemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(CpdSemType, rendererCpdCell);
    configRenderers.setCellRenderer(InVivoSemType, new InVivoGridCellRenderer());
    configRenderers.setCellRenderer(InSilicoSemType, new InSilicoGridCellRenderer());
    configRenderers.setCellRenderer(IC50SemType, new IC50GridCellRenderer());
    configRenderers.setCellRenderer(IdeaSemType, new IdeaGridCellRenderer());
    configRenderers.setCellRenderer(ActivityConcentrationSemType, new ActivityConcentrationGridCellRenderer());
    configRenderers.setCellRenderer(QNumDateSemType, new QNumDateGridCellRenderer());
    configRenderers.setCellRenderer(TextDateSemType, new TextDateGridCellRenderer());
    configRenderers.setCellRenderer(PrimitiveSemType, new PrimitiveGridCellRenderer());
    configRenderers.setCellRenderer(CpdPrimitiveSemType, new CpdPrimitiveGridCellRenderer());

    //const mapGridHeaderRenderers = new ClassMap();
    configRenderers.setColHeaderRenderer(SemType, new MetaDataGridColHeaderRenderer());//new ButtonGridColumnHeaderRenderer());
    //configRenderers.setColHeaderRenderer(CpdSemType, new CpdColumnHeaderRenderer());
    configRenderers.setColHeaderRenderer(NullSemType, new EmptyButtonGridColumnRenderer());

    configRenderers.setDeffaultColHeaderRenderer(new ButtonGridColumnHeaderRenderer());

    configRenderers.setColRowHeaderRenderer(new CpdColumnHeaderRenderer());

    const bHasMetaData = DGApp.hasVirtualColumns(dframe);
    if(!bHasMetaData)
    {
        const dialProgress = new ProgressDialog("Opening MedChem Grid...", true);
        await dialProgress.show();

        const arDescriptors = await AnalysisLoader.populateMetaData(dframe);
        await IntuenceDiscovery.processAnalysisData(dframe, configSemTypes, _package, dialProgress);

         const lstTViews = grok.shell.tableViews;
         for(let tview of lstTViews)
         {
             DGApp.setVistualColumnsVisible(tview.grid, false);
         }

        if(m_handlerViewAddedHideVirtual === null) {
            m_handlerViewAddedHideVirtual =  grok.events.onViewAdded.subscribe((view) => {
                if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
                    DGApp.setVistualColumnsVisible(view.grid, false);
                }
            });
        }

        if(m_handlerViewLayoutHideVirtual === null) {
            m_handlerViewLayoutHideVirtual = grok.events.onViewLayoutApplied.subscribe((layout) => {
                    const view = layout.view;
                    const view1 = view.c$;

                    if (view1.grid !== undefined)
                        DGApp.setVistualColumnsVisible(view1.grid, false);
                }
            );
        }

        await dialProgress.close();
    }
    await IntuenceDiscovery.open(view, configRenderers, arFilterClasses, _package, fnFilterOpener, fnFilterCloser);


    return true;
}



//tags: appp
//name: Super Primitives App Demo
export async function startSuperPrimitivesApp() {
    let tv = grok.shell.addTableView(grok.data.demo.demog());

    const viewer = tv.filters(
        {filters: [
            {type: DG.FILTER_TYPE.HISTOGRAM, column: 'subj'},
            {type: 'gnfentity:compoundFilter', column: 'site' },
        ]}
    );

    let aaa = 0;
/*
    const filter = compoundFilter();
    filter.attach(viewer.dataFrame);
    FiltersUtils.addFilter(filter, viewer);
*/

    //FiltersUtils.showFilters(false);

}


//tags: app
//name: Primitives App Demo
export async function startPrimitivesApp() {

    //const pc = new PinnedColumn(null);

    const nRowCount = 100;
    const nColCount = 500;
    const dframe = grok.data.demo.randomWalk(nRowCount, nColCount);

    const colDate = dframe.columns.addNewDateTime('datetime');
    const colText = dframe.columns.addNewString('Text');

    let col = dframe.columns.byIndex(0);
    col.setTag("TestTag", "TestValue");
    col.set(0, 1.0);
    col.set(1, 2.0);
    col.set(2, 3.0);
    col.set(3, 4.0);
    col.set(4, null);

    col = dframe.columns.byIndex(2);
    col.set(0, null);
    col.set(1, NaN);

    const fNaN = col.get(1);

    colText.set(4, "aaaaa");

    col = dframe.columns.byIndex(0);
    let nTime = null;
    for(let nRow=5; nRow<nRowCount; ++nRow)
    {
        col.set(nRow, 2.2);

        nTime = RandUtils.generateRandomDate(2005, 2019);
        colDate.set(nRow, DG.DateTime.fromMillisecondsSinceEpoch(nTime));

        colText.set(nRow, "bbbbb" + nRow);
    }


    const configRenderers = new GridRenderersConfig();
    const rendererRowHeader= new ButtonGridRowHeaderRenderer();
    configRenderers.setRowHeaderRenderer(rendererRowHeader);
    //const rendererRowHeader= new OCLCpdGridRowHeaderRenderer();

    configRenderers.setCellRenderer(SemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(ObjectSemType, new ObjectGridCellRenderer());
    configRenderers.setCellRenderer(PrimitiveSemType, new PrimitiveGridCellRenderer());

    //const mapGridHeaderRenderers = new ClassMap();
    configRenderers.setColHeaderRenderer(SemType, new ButtonGridColumnHeaderRenderer());
    configRenderers.setColHeaderRenderer(NullSemType, new EmptyButtonGridColumnRenderer());

    configRenderers.setDeffaultColHeaderRenderer(new ButtonGridColumnHeaderRenderer());

    const arFilterClasses = [ComboBoxFilter, CheckBoxesFilter, RecentTimeFilter, RangeSliderFilter, MultilineTextFilter];
    let viewSS = null;
    try{viewSS = grok.shell.addTableView(dframe);}
    catch(e)
    {
        throw new e;
    }

    const opts = viewSS.grid.getOptions(true);
    const nHColHeader = opts.look.colHeaderHeight;
    const nHRow = opts.look.rowHeight;

    viewSS.grid.setOptions({
        colHeaderHeight: 30,
        rowHeight: 25
    });


    let colGridHtml = viewSS.grid.columns.byIndex(3);
    colGridHtml.cellType = 'html';
    colGridHtml.format = null;

    viewSS.grid.onCellPrepare(function (gc) {
        if (gc.isTableCell && gc.gridColumn.name == colGridHtml.name) {
            gc.style.element = ui.divText(gc.cell.value === null ? "" : gc.cell.value.toString());
        }
    });


    const col5 = viewSS.grid.columns.byIndex(5);
    col5.format = '<div>Test</div>';
    col5.cellType = 'html';

    //const lstCs = viewSS.grid.root.getElementsByTagName("canvas");
    const eCanvasTable = viewSS.grid.overlay;
    grok.events.onContextMenu.subscribe((args) => {

   // rxjs.fromEvent(lstCs[2], 'contextmenu').subscribe((e) => {
        const grid = args.args.context;
        const e = args.causedBy;
        const cell = viewSS.grid.hitTest(e.offsetX, e.offsetY);
        if(cell === undefined || cell.cellType === null)//bug in DG , top left cell
            return;

        const elem = document.elementFromPoint(e.clientX, e.clientY);
        const colGrid = cell.gridColumn;
        if(!PC.isPinnedColumn(colGrid) && (cell.isTableCell || cell.isColHeader) && colGrid.cellType !== "html" && (elem === viewSS.grid.canvas || elem === viewSS.grid.overlay))
        {
            const menu = args.args.menu;
            menu.item("Pin Column", () => {

                const headerRows = new GridRowHeaderNew(colGrid);
                //new PinnedColumn(colGrid);
            });
            e.preventDefault();
            e.stopPropagation();
            return;
        }
    });


    //TableViewLayoutManager.install(new FlowTileLayoutManager());

    //if(true)
      //  return;

    viewSS.grid.addProperty("Test", DG.TYPE.STRING, 'aaa');
    const arPs = viewSS.grid._properties;
    for(let n=0; n<arPs.length; ++n) {

    let p = DG.toJs(arPs[n].dart);
    let v = arPs[n].defaultValue;
    if(v === "aaa")
    {
        console.log();
    }

}
    /*
    const viewerSP = viewSS.scatterPlot();
    const manager = viewerSP.view.dockManager;
    const nodeSP = manager.findNode(viewerSP.root);
    manager.close(viewerSP.root);
    nodeSP.container.destroy();
    manager.dock(viewerSP.root, DG.DOCK_TYPE.RIGHT, manager.rootNode, "SP", 0.11);
*/

    //my changes await DGApp.open(_package, viewSS, configRenderers, arFilterClasses);

    TableViewLayoutManager.install(new FlowTileLayoutManager());



    grok.events.onViewerAdded.subscribe((args) => {
        const viewer = args.args.viewer;

        if(viewer.type == "Grid") {
            viewer.setOptions({
                colHeaderHeight: 30,
                rowHeight: 25
            });

            //const settings = viewer.columns.byIndex(1).settings;
            //viewer.columns.byIndex(1).settings = {a: "Test"};
       }

        let fff = 0;
    });






    viewSS.grid.columns.byIndex(0).visible = false;

    const dlg = ui.dialog("Layout");
    dlg.addButton("Save", async function (){

        const layout = viewSS.saveLayout();

        let rrr = layout.id;
        layout.id = "13b27150-1890-11ec-a7b9-eff30b5796ed";
        layout.setUserDataValue("myvalue", "99");

        const json = layout.toJson();
        await grok.dapi.layouts.save(layout);

    });

    let viewerFilters = null;
    dlg.addButton("Filters", async function () {

        viewerFilters = viewSS.filters(/*{
            filters: [
                {type: 'compoundFilter', column: '#2'},
                {type: DG.FILTER_TYPE.FREE_TEXT},
            ]
        }*/);

        addCompoundFilter(viewerFilters);
    });

    dlg.addButton("Remove Cpd", async function () {

        removeCompoundFilter(viewerFilters);
    });


    dlg.addButton("Load", async function (){
        //ViewerLayoutManager.uninstall();
        const lu = await grok.dapi.layouts.find("13b27150-1890-11ec-a7b9-eff30b5796ed");
        viewSS.loadLayout(lu);
        let val = (await lu).getUserDataValue("myvalue");
        let gfhfh = 0;

    });


    dlg.show();


    grok.events.onViewLayoutGenerated.subscribe((layout) => {
        const strJSON = layout.toJson();
        const obJSON = JSON.parse(strJSON);
        TextUtils.printValues(obJSON);

        const viewTmp = layout.view;
        const viewNew = DG.toJs(viewTmp);
        const viewTmpSh = grok.shell.v;
        const nodeRoot = grok.shell.dockManager.rootNode;//viewTmpSh.dockNode;
        const lst = document.querySelectorAll ('.d4-viewer');
        const set = new Set();
        DGUtils.iterateDockNodes(nodeRoot, set, 0);
        let e = null;
        const arViewers = Array.from(viewTmpSh.viewers);
        for(let v of viewTmpSh.viewers) {

            const parent = v.root.parentNode;
            const parentOffset = v.root.offsetParent;
            const nW = v.root.style.width;
            let nE=0;
            for(; nE<lst.length; ++nE)
            {
                e = lst.item(nE);
                if(e === v.root)
                {
                    break;
                }
                let bbb = nE;
            }

            let aaa = 0;

                    //if(v.root.offsetParent === null)
                      //  v.close();
                }



        let aaa = 0;
    });

    grok.events.onViewLayoutApplying.subscribe((layout) =>  {
        let strJSON = layout.toJson();
       // strJSONstrJSON.replaceAll()
        //const map = await layout.getProperties();
        //layout.viewState
    let aaa = 0;
    });

    grok.events.onViewLayoutApplied.subscribe((layout) => {
        const itViewers = layout.view.viewers;
        const arViewers = Array.from(itViewers);

        let viewer = null;
        for(let n=0; n<arViewers.length; ++n) {
            viewer = arViewers[n];
            if(viewer.type !== "Grid")
                continue;

            viewer.setOptions({
                colHeaderHeight: 30,
                rowHeight: 25
            });

            const colG = viewer.columns.byIndex(5);
            try {
                colG.format = '<div></div>Test<div>';
                colG.cellType = 'html';
            }
            catch(e) {
                throw e;
            }

            //PC.installPinnedColumns(viewer);

/*
            let colGrid = null;
            let settings = null;
            const lstCols = viewer.columns;
            const arPinnedCols = new Array();

            for(let nCol=0;nCol<lstCols.length; ++nCol) {
                colGrid = lstCols.byIndex(nCol);
                settings = colGrid.settings;
                if(settings !== null && settings !== undefined && settings.isPinned) {
                    arPinnedCols.push(colGrid);
               }
            }

            arPinnedCols.sort((colOne, colTwo) => {
               if(colOne.settings.idxPinned === colTwo.settings.idxPinned)
                   throw new Error("Pinned indices cannot be equal for different columns");

               return colOne.settings.idxPinned < colTwo.settings.idxPinned ? -1 : 1;
            });

            for(let n=0;n<arPinnedCols.length; ++n) {
                colGrid = arPinnedCols[n];
                let header = new GridRowHeaderNew(colGrid);
            }
*/
            let aaa = 0;
        }
        //ViewerLayoutManager.install(new FlowTileLayoutManager());
        //const viewTmp = layout.view;
/*
        for(let v of viewTmp.viewers) {
            if(v.root.offsetParent === null)
                v.close();
        }*/
    });

}

async function Testing() {
    //const nA = a();
    //b();
    //c();


}







//tags: app
//name: GNF Test App
export async function test(s) {

    const view = grok.shell.addTableView(grok.data.demo.demog());
    let col = view.grid.columns.byName('disease');
    view.grid.setOptions({'rowHeight': 100});
    col.width = 200;
   //col.format = `<div style="display:flex; flex-direction:column; padding: 5px"><div><div>sex: <b>#{t.row[sex]}</b></div><div> height: <span style="background-color:#{t.color(height)}"><b>#{t.row[height]}</b></span> </div><button class="ui-btn ui-btn-ok ui-btn-raised">CONTACT</button></div></div>`;
    col.format = '<div class="intuence-discovery-edu-links-format"><a href="http://novartis.com" target="_blank">#{t.row[height]}</a></div>';
    col.cellType = 'html';

    if(true)
        return;

    //const df = grok.data.demo.molecules();//grok.shell.addTableView(grok.data.demo.molecules());
    //const view = grok.shell.addTableView(df);
    //const col0 = view.grid.columns.byIndex(0);
    //col0.visible = false;

    if(true)
        return;

    const idLayout = "13b27150-1890-11ec-a7b9-eff30b5796ed";

    const dlg = ui.dialog("Layout");
    dlg.addButton("Save", async () => {

        const layout = view.saveLayout();
        layout.id = idLayout;
        await grok.dapi.layouts.save(layout);
    });

    dlg.addButton("Load", async () => {
        const layout = await grok.dapi.layouts.find(idLayout);
        view.loadLayout(layout);
    });

    dlg.show();

    if(true)
    return;

    const options = view.grid.getOptions(true);
    const nRH =  options.look.rowHeight;
    const nCH = options.look.colHeaderHeight;

    console.log("Row Height: " + nH);

    if(true)
        return;

    const lstCs = view.grid.root.getElementsByTagName("canvas");
    for(let nC=0; nC<lstCs.length; ++nC){
            lstCs[nC].style.left = "100px";
    }


    const eCanvasRowHeader = ui.canvas(100, 1500);
    lstCs[1].parentNode.insertBefore(eCanvasRowHeader, lstCs[1]);

    rxjs.fromEvent(view.root, 'contextmenu').subscribe((e) => {
        let menu = DG.Menu.popup();

        /*e.screenX
        e.clientX
        e.offsetX
        e.x*/
        const elem = document.elementFromPoint(e.clientX, e.offsetY);
        const b = elem === eCanvasRowHeader;
        let rc = eCanvasRowHeader.getBoundingClientRect();
        if(b) {
            menu = menu.item("Test Bug", (str) => {
                console.log("clicked");
            });

            menu.show();
        }
    });


    if(true)
        return;


    const observable = new Observable(function subscribe(subscriber) {
        // Keep track of the interval resource

        let aaa = 0;

        const intervalId = setInterval(() => {
            subscriber.next('hi');
        }, 1000);

        // Provide a way of canceling and disposing the interval resource
        return function unsubscribe() {
            clearInterval(intervalId);
        };
    });

    const subscription1 = observable.subscribe(x => console.log("from 1 " + x));
    const subscription2 = observable.subscribe(x => console.log("from 2 " + x));


    //await Testing();

    let rrrrrr = 0;

    if(true)
        return;


    let dframe = grok.data.demo.randomWalk(100, 10);
    let griddd = grok.shell.addTableView(dframe);


    //let col = null;
    try{
        col = dframe.columns.addNewVirtual('Virt Column', function (nRow)
        {
            if(nRow === undefined || nRow === null)
            {
                let ffhfhf = 0;
            }

            const bV = dframe.columns.byName("Virt Column").isVirtual;

            return 3;
        });


        col = dframe.columns.addNewVirtual('Virt Column 1', function (nRow)
        {
            if(nRow === undefined || nRow === null)
            {
                let ffhfhf = 0;
            }


            return 5;
        });

    }
    catch(e)
    {
        let err = e.message;
    }

    const bVirt = col.isVirtual;



    let ffhfhg = 0;

    return;


    //const dial = ui.dialog("Analysis Loader");
    //dial.show();
       /*
    const strURL = "https://cdf-ew.nibr.novartis.net/services/v1/selectors/metadata/ComparableFactSet/7955.json";
    const obSON = await WebUtils.makeGetRequest(strURL);
    await MetaDataLoader.parseJSON(obSON);

    return;
         */
    let nRowCount = 100;
    let nColCount = 10;
    const dframeRand = grok.data.demo.randomWalk(100, 5);

    let gridd = grok.shell.addTableView(grok.data.demo.randomWalk(100, 100)).grid;
    let dial = ui.dialog("BUG Demo");

    dial.addButton("Click", () => {
        gridd.setRowOrder([1, 56, 3, 6, 4]);
     });

    dial.show();
    return;

    dframeRand.columns.byIndex(0).set(0, null);
    dframeRand.columns.byIndex(0).set(1, null);
    dframeRand.columns.byIndex(0).set(2, null);
    dframeRand.columns.byIndex(0).set(3, null);
    dframeRand.columns.byIndex(0).set(4, null);


    const colDate = dframeRand.columns.addNewDateTime('datetime');
    const colFloat = dframeRand.columns.addNewFloat('float');
    let nTime = -1;
    for(var n=5; n<nRowCount; ++n)
    {
        colFloat.set(n, 2.2);

        nTime = RandUtils.generateRandomDate(2005, 2019);
        colDate.set(n, DG.DateTime.fromMillisecondsSinceEpoch(nTime));
    }


    let grid = grok.shell.addTableView(dframeRand).grid;

    return;

    grid.onRowsResized.subscribe((ev) => {
        const options = grid.getOptions();
        const nH =  options.look.rowHeight;
        console.log("Row Height: " + nH);
    });

    let nCallCount = 0;
    grok.events.onViewerAdded.subscribe((args) => {
        ++nCallCount;
       console.log("Call Count: " + nCallCount);
    });

    grok.events.onViewerClosed.subscribe((args) => {
        console.log("Never gets called.");
    });


     return;
    let typeCompose = CDFComposites.getCompositeType("ActivityConcentration");
    let typePrimi = typeCompose.getPrimitiveType("activity");
    let typeParent = typePrimi.getParentType();
    let bSame = typeCompose === typeParent;

    //const strQuery = "{[\"structureRef,Label:Structure\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaId'].tValue,Label:CAST Idea ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'EJ3JyFsegb3dEVKbX4NwVSlgCPc')].compositeFieldInstances['9KsNABsTiq_smUmzLavqSFSnvbc'].activity,Aggregator:AVG,Label:IL17AA LZV623\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'EJ3JyFsegb3dEVKbX4NwVSlgCPc')].compositeFieldInstances['3XX9i9bjaIvD_CBDqikn9-lbudI'].dValue,Aggregator:MAX,Label:IL17AA LZV623 exp date\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'EJ3JyFsegb3dEVKbX4NwVSlgCPc')].compositeFieldInstances['X66Jkio5VXbEftpci6LA7LY2gjo'].imageId,Label:IL17AA LZV623 curve\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZoUwo87Ne896wGBnOb_XE-gtq7k')].compositeFieldInstances['6XpYwIht4MoHhd5OccXW7KTIVGE'].activity,Aggregator:AVG,Label:IL17AA IL17RA v2\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZoUwo87Ne896wGBnOb_XE-gtq7k')].compositeFieldInstances['qJo7i9V7-shoZylhm3vIqIQigjU'].imageId,Label:IL17AA IL17RA v2 curve\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'AIBBf7UwXJSvlsGB69HM2Dy78zU')].compositeFieldInstances['Ln82S1KvDZAiUqp5080krglYndg'].activity,Aggregator:AVG,Label:IL17AF IL17RA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ghCxfuDXaFA8RRlCVNVW8XiQFPo')].compositeFieldInstances['GGgFRl54EKEGLqB-wUCDLJXW0kw'].activity,Aggregator:AVG,Label:IL17 cell IL6\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'jbttmnOwP7eaHm64JvXe8RkChhM')].compositeFieldInstances['oPRlls9CnNzsKrDCge0Ai_YcbgE'].activity,Aggregator:AVG,Label:IL6 cell PP shift IL6\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rpyQkQjkGFF5Nc1XdmgymOCQJXw')].compositeFieldInstances['Sw-_YO8ZvRrn5WEpabhNPa_3MlQ'].activity,Aggregator:AVG,Label:IL17 cell rat\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'BcLJqvuV5zAxvRw8SQd3gNHpGgY')].compositeFieldInstances['_yWmg6nSl9bYOKE9QMLtxgF3D7k'].nValue,Aggregator:AVG,Label:IL17A SPR KD\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'BcLJqvuV5zAxvRw8SQd3gNHpGgY')].compositeFieldInstances['aLgA4eidTm1xvqc78KnQ-v9VvQE'].nValue,Aggregator:AVG,Label:IL17A SPR Kon\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'BcLJqvuV5zAxvRw8SQd3gNHpGgY')].compositeFieldInstances['ZN5gRV2H6vKHvPrJwHnFX2N7xYM'].nValue,Aggregator:AVG,Label:IL17A SPR Koff\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'x-OoKi2Vow7clPR2ES-KinMqsm8')].compositeFieldInstances['smfAvailability'].amount,Label:SMF Availability Availability Amount\",\"smallMoleculeSampleRefs.labhead,Label:Lab Head 521\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '7sTLZPVAcxsj7wzjnpNn8WKCEZo')].compositeFieldInstances['DJAy7t-40IIlp68XSfdM4wa8TW8'].tValue,Label:MAP in vivo PK DOG Route Admin\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '7sTLZPVAcxsj7wzjnpNn8WKCEZo')].compositeFieldInstances['K7mar7oX4FgddC0YfyfLwsznjwc'].nValue,Label:MAP in vivo PK DOG Dose\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '7sTLZPVAcxsj7wzjnpNn8WKCEZo')].compositeFieldInstances['H9qhyrc44YnwGNVhy6I3jyw1q2Y'].nValue,Label:MAP in vivo PK DOG CL\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['PZgPepgGKGP8YIoA4S9xTgjCRpM'].nValue,Aggregator:AVG,Label:Passive Permeability with MDCK-LE V2 PappAtoB\",\"smallMoleculeSampleRefs.preferredName,Label:Sample ID\",\"smallMoleculeSampleRefs.registrationDate,Label:Registration Date\",\"smallMoleculeSampleRefs.registrationDate,Aggregator:MIN,Label:First Registered\",\"smallMoleculeSampleRefs.chemist,Aggregator:MOST_COMMON,Label:Chemist\",\"smallMoleculeSampleRefs.labNotebookId,Label:Lab Notebook\",\"stereoCategory,Label:Stereo Category\",\"stereoDescriptor,Label:Stereo Descriptor\",\"smallMoleculeSampleRefs.projectCode,Aggregator:UNIQUE,Label:Project Code\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yo_8k8A73kiguu9TfV8SAqx3xn4')].compositeFieldInstances['Rah7fnilYoxxSfgrSBt5IJh21G0'].nValue,Aggregator:GEO,Label:HT Solubility pH6.8 [mM]\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yo_8k8A73kiguu9TfV8SAqx3xn4')].compositeFieldInstances['gK8vUaZLRgOwoTYA9SgyYWMjYHY'].nValue,Aggregator:GEO,Label:FaSSIF\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yo_8k8A73kiguu9TfV8SAqx3xn4')].compositeFieldInstances['3xfR8KjMIPd_8iKGHo3Jmsqtw_s'].dValue,Aggregator:MAX,Label:HT Solubility Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Sylv79kBc3leD0U-cR3A_63NuEc')].compositeFieldInstances['3TR05dUgcyFmyVvWmqzVYVjxLS0'].nValue,Aggregator:GEO,Label:SF-Solubility pH6.8 [mM]\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Sylv79kBc3leD0U-cR3A_63NuEc')].compositeFieldInstances['KIblChoT9LSMupVkStE0zpvPmgY'].tValue,Aggregator:MOST_COMMON,Label:SF-Solubility Center\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Sylv79kBc3leD0U-cR3A_63NuEc')].compositeFieldInstances['iOWQ8xL5eQdJ7dlMa69oUVnx4JM'].dValue,Aggregator:MAX,Label:SF-Solubility Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rV32B58AH_LPdsULrc435vhynKA')].compositeFieldInstances['0hvfnRo8FbmiHkqu_DHrR2WE7h8'].nValue,Aggregator:AVG,Label:logPAMPA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rV32B58AH_LPdsULrc435vhynKA')].compositeFieldInstances['vV7V9E8Ok6KN215nf6UGwebrOCA'].dValue,Aggregator:MAX,Label:logPAMPA Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['qkZ_DQyzhsRf3ECRCPOVBI9_74M'].nValue,Aggregator:GEO,Label:LE-MDCK PappAtoB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['bBVic8BCg-xxYUflnYTZMFklCm0'].activity,Aggregator:AVG,Label:LE-MDCK Recovery\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['u650cjE9QpyHLkRQkRzxXyezDic'].tValue,Aggregator:MOST_COMMON,Label:LE-MDCK FA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['1SpmIY9A9F5SA6lV2g0aYUCFzbE'].dValue,Aggregator:MAX,Label:LE-MDCK Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['PZgPepgGKGP8YIoA4S9xTgjCRpM'].nValue,Aggregator:GEO,Label:LE-MDCK v2 PappAtoB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['x6ZLFG75FELMp41BRumoi_WyqjQ'].activity,Aggregator:AVG,Label:LE-MDCK v2  Recovery\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['hRca1-gc2nwBGAaKGvGo9ZoHm1c'].activity,Aggregator:AVG,Label:LE-MDCK v2 FA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['bS76DvBdPdICdE1dda211D02b60'].dValue,Aggregator:MAX,Label:LE-MDCK v2 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['A--0FCuOBBnNz0D5iIxJLPHj3PE'].nValue,Aggregator:GEO,Label:MDCK-MDR1 ER\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['BpfvavyJXK0obRR_3kR07dCXHgU'].nValue,Aggregator:AVG,Label:MDCK-MDR1 A-B\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['jU0XfZmFjh1Mb7GaqbOHUcJXsg8'].nValue,Aggregator:AVG,Label:MDCK-MDR1 B-A\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['Sni9_j1Qo0HrQcLGfKTbvctwGbo'].dValue,Aggregator:MAX,Label:MDCK-MDR1 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'cGEAi42quRk6vF554y41YFSfQVw')].compositeFieldInstances['MzZ3k5P-wW0cYMLFsDvcMJBEPVw'].tValue,Aggregator:MOST_COMMON,Label:OATP1B1 Substrate Assay\",\"length(most_common(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'cGEAi42quRk6vF554y41YFSfQVw')].compositeFieldInstances['MzZ3k5P-wW0cYMLFsDvcMJBEPVw'].tValue)),Label:OATP1B1 Substrate\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'cGEAi42quRk6vF554y41YFSfQVw')].compositeFieldInstances['HQAB_B4nA-VOh-lZe4cOVkY2xhU'].dValue,Aggregator:MAX,Label:OATP1B1 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'd0KLAfcOM8Vr5yzDaVSA-obfmi0')].compositeFieldInstances['s-ySb1q6vPopbIyEWnDRDM03JuY'].nValue,Aggregator:AVG,Label:HT LogP\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'd0KLAfcOM8Vr5yzDaVSA-obfmi0')].compositeFieldInstances['0YiXvrJFA9Lg2P0SbfdGXIFFKsk'].dValue,Aggregator:MAX,Label:HT LogP Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['nFO-P99OG_Ct6D7yJcZR2MOwPGA'].nValue,Aggregator:AVG,Label:Direct LogP\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['l38JNOlK9YR1XcHRQZgedBL074Y'].nValue,Aggregator:AVG,Label:Direct LogD\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['1wlD76Y3geUcXarB36a8rTWHYeM'].dValue,Aggregator:MAX,Label:Direct LogPD Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['B1Kx8dAmAESo-Wtc69gWCxUCp98'].nValue,Aggregator:MAX,Label:pKa1\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['yqGqf8hPiyf3M4QUwRXg2TIWiGI'].nValue,Aggregator:MAX,Label:pKa2\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['T2o8Lm5MtXZRvkvhSKKG6yHWV6U'].nValue,Aggregator:MAX,Label:pKa3\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['3_5biEAay1ekL8C9Rpf4wgfMpYM'].nValue,Aggregator:MAX,Label:pKa4\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['6ORyt6EWxGIpdHozCN6e0uPq12o'].tValue,Aggregator:MOST_COMMON,Label:pKa1_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['s2Q4yAR7ipaCO-8BR6ZOgpOa0jE'].tValue,Aggregator:MOST_COMMON,Label:pKa2_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['cdBoHQs26Hrzt3Hi7poPXh7VMKE'].tValue,Aggregator:MOST_COMMON,Label:pKa3_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['eTtiNX_ZYRIQSG2O1ZDPZazqqfE'].tValue,Aggregator:MOST_COMMON,Label:pKa4_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'SW0kTdONJDEcF_KnycpS2J-18yI')].compositeFieldInstances['Zs46VVlLT4ziyZdQEeQgfi2Sew4'].nValue,Aggregator:GEO,Label:rLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'SW0kTdONJDEcF_KnycpS2J-18yI')].compositeFieldInstances['mHLIaq5yQGqSfi-rkC7ByMdeLD8'].dValue,Aggregator:MAX,Label:rLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'kQkCHBBEasvjH8Y8xNVI8dAQNnA')].compositeFieldInstances['QC2COEqmw2WVQdgtNhTNhL4U4Dw'].nValue,Aggregator:GEO,Label:mLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'kQkCHBBEasvjH8Y8xNVI8dAQNnA')].compositeFieldInstances['nyUdZpolM9dQD75Q2TiZ2Vncc58'].dValue,Aggregator:MAX,Label:mLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '-vddsB-WtoW7jlwQQTBPTGZgajs')].compositeFieldInstances['bMYHFdNUbjUgyc8p5bSpI4LjiOE'].nValue,Aggregator:GEO,Label:hLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '-vddsB-WtoW7jlwQQTBPTGZgajs')].compositeFieldInstances['AtpVqbOtagB8vxmy8JWK89ALDHs'].dValue,Aggregator:MAX,Label:hLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'hH4hg1T_m42_ZCJDq8mqXe20O0E')].compositeFieldInstances['ezEHQQMyhsECazauG8LtU77Dt0I'].nValue,Aggregator:GEO,Label:dLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'hH4hg1T_m42_ZCJDq8mqXe20O0E')].compositeFieldInstances['BpJmlYK02pFrIjhbZ1J86yETD7o'].dValue,Aggregator:MAX,Label:dLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Gq5S0jGMIhj7dGpdRcrDvb6giII')].compositeFieldInstances['yfLuBonQy0Ny8k25UEsd3Tul6MQ'].nValue,Aggregator:GEO,Label:cynoLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'e3ZPOcb6VhVEWZqtHhAEppTYnj0')].compositeFieldInstances['MhNRI6-NppT1nfK--Vcve_umwjA'].nValue,Aggregator:GEO,Label:minipigLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'm6kM9DPXPpgZs5lHgcBPlVUnjm4')].compositeFieldInstances['OPoAXu6ckjvDTF8JdqlrpQan1co'].nValue,Aggregator:GEO,Label:rHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LMl5PyWaQXt8CW1_7DrzTKV357E')].compositeFieldInstances['XiODJf-kNWqB4QgJjWlO2LMUWYU'].nValue,Aggregator:GEO,Label:mHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Ed-roIcf615pfStMEHvW6n903_M')].compositeFieldInstances['Yf6hOGHcEgnPyyU2YOXln4zzbrY'].nValue,Aggregator:GEO,Label:hHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yv_0ZfUg0PkZvmOTaaSKnQmQX4M')].compositeFieldInstances['dTaypxGOk-4qkFsGzz5s0ra3Uqk'].nValue,Aggregator:GEO,Label:dHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'nrfwl7yWsZ6XMkSVVH-_kgxRc3Y')].compositeFieldInstances['wHzQhi_zYRJ-vDMHMicKN_fAVeI'].nValue,Aggregator:GEO,Label:cynoHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'pxgyV6XiF6sKzutBPnAtWYMYoKQ')].compositeFieldInstances['zwfMDuceMoZkWpfavBfC_610sPI'].nValue,Aggregator:GEO,Label:minipigHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'mVivAeeHG0gQSlqXK8XU6FzJFxU')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:rPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '0t8LmVLFGFeq0iKTtUOxslgglGI')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:mPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'nWPGpBxxnEXWP2MMeo4sLbiG158')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:hPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'sTfm8emCiu6f92PfW5NjF7qY-Cs')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:dPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rV5INo9qP628SaP9FBPcN8GqGg8')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:cynoPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 't8kq0fsl7JA2OUguGwr144LuDMk')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:minipigPPB PB\",\"divide(remove_qual(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'eRvF0tyNCvVnvCLeUArIn2SOKaU')].compositeFieldInstances['iNK51tXhvRVP1M-Lys_0hqviAyQ'].nValue)),100),Label:rfu(mic)\",\"divide(remove_qual(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'sXzxu6x4_DxkM27u7Yd4yXyLHes')].compositeFieldInstances['iNK51tXhvRVP1M-Lys_0hqviAyQ'].nValue)),100),Label:hfu(mic)\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['dL_0aWPbTKonqnGqZQvOHWw26yk'].nValue,Aggregator:GEO,Label:CYP3A4 TDI kObs\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['_ib14_aeIFm_ByUaaLTsGLWVyEA'].tValue,Aggregator:MOST_COMMON,Label:CYP3A4 TDI Inh-Rev\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['AXW4CGHn41mJg5nfaJgWVZiuEVM'].tValue,Aggregator:MOST_COMMON,Label:CYP3A4 TDI Flag\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['7GamtjRP9RigQ6oBFjZMC9IcxhI'].dValue,Aggregator:MAX,Label:CYP3A4 TDI Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zCjqunJldyDxOKs6xxOEzchLOO8')].compositeFieldInstances['NORxqb6pWHCmLHs62f6Q57ChKGE'].activity,Aggregator:GEO,Label:CYP3A4 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zCjqunJldyDxOKs6xxOEzchLOO8')].compositeFieldInstances['5LnswVHMIVO98Q1ZTG7dcPzdiA0'].dValue,Aggregator:MAX,Label:CYP3A4 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'E3UyRf8YXKf_ysiAjo6KImmtsaE')].compositeFieldInstances['vpyt3dghj1VQ9thRte6GcHwG0MY'].activity,Aggregator:GEO,Label:CYP2D6 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'l-zKqbBDyJv2Yoian60ZOVRDVvc')].compositeFieldInstances['rIBFbovjkIUZq_SOXFSzbREW-0g'].activity,Aggregator:GEO,Label:CYP2C9 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['izC7ons3rj0jpiiW_Q-ptyt7D4w'].activity,Aggregator:GEO,Label:hERG Binding IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['oTttI7IzFIAabLCS8NS425jHTJg'].activity,Aggregator:AVG,Label:hERG Binding inh@1uM\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['S74KXMvFXsal2nN8GNZifKeQCsM'].activity,Aggregator:AVG,Label:hERG Binding inhib@30uM\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['izC7ons3rj0jpiiW_Q-ptyt7D4w'].imageId,Label:hERG Binding Curve ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['T4gbZBctaeUMKgLYIK0Zrpko_z4'].dValue,Aggregator:MAX,Label:hERG Binding Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'TTyaVjvh-yhQDBlokdjq_DVymtU')].compositeFieldInstances['1y5Lk8vxAlXAGJCd0cEVszsRGpw'].activity,Aggregator:GEO,Label:hERG QPatch IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'TTyaVjvh-yhQDBlokdjq_DVymtU')].compositeFieldInstances['1y5Lk8vxAlXAGJCd0cEVszsRGpw'].hill,Aggregator:MIN,Label:hERG QPatch Hill slope\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'TTyaVjvh-yhQDBlokdjq_DVymtU')].compositeFieldInstances['YSI1gn-QJHiaGU-MdfHhr37axyM'].dValue,Aggregator:MAX,Label:hERG QPatch Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZGU-0ne7HW7_9_DdJZxlL8zGc-o')].compositeFieldInstances['ODcF86fJwqifZ9ChnTXCNQHeSWQ'].activity,Aggregator:GEO,Label:BSEP IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZGU-0ne7HW7_9_DdJZxlL8zGc-o')].compositeFieldInstances['Au5cnkk1XysekBtjLZ2IuoKrcRY'].dValue,Aggregator:MAX,Label:BSEP Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '9YZ73Ke-sVcQuT29Wg1X9Stw1Rs')].compositeFieldInstances['JSj7zuEnJ2nrXCWdDu-533QqWaQ'].activity,Aggregator:GEO,Label:Nav1.5 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '9YZ73Ke-sVcQuT29Wg1X9Stw1Rs')].compositeFieldInstances['JSj7zuEnJ2nrXCWdDu-533QqWaQ'].imageId,Label:Nav1.5 Curve ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '9YZ73Ke-sVcQuT29Wg1X9Stw1Rs')].compositeFieldInstances['Mx2kdAbsw22_o-tjYbeujYiW1lc'].dValue,Aggregator:MIN,Label:Nav1.5 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'RI0c7xNdC0J2askGjmQeJuDoUHc')].compositeFieldInstances['mi4-u0Oq3M6FWOTx_MtCNSK56qg'].activity,Aggregator:GEO,Label:Cav1.2 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'RI0c7xNdC0J2askGjmQeJuDoUHc')].compositeFieldInstances['OqzK0_SCHVo4sGpkJiLa1vg5OZA'].dValue,Aggregator:MAX,Label:Cav1.2 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ze8J1qI_fscCkAewFBUUjXgJjxU')].compositeFieldInstances['WZJY9Rn7klxeEkUpTzmk4sPBZCc'].activity,Aggregator:GEO,Label:KCNQ1 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ze8J1qI_fscCkAewFBUUjXgJjxU')].compositeFieldInstances['WZJY9Rn7klxeEkUpTzmk4sPBZCc'].imageId,Label:KCNQ1 Curve ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ze8J1qI_fscCkAewFBUUjXgJjxU')].compositeFieldInstances['GFJ_lAiuTRTlpzi8aZQNBeMt9DA'].dValue,Aggregator:MAX,Label:KCNQ1 Date-out\",\"averageMass,Label:Average Mass\",\"tpsa,Label:TPSA\",\"nRotBonds,Label:Num Rotatable Bonds\",\"nHeavyAtoms,Label:Num Heavy Atoms\",\"factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logp'].numericPrediction,Aggregator:AVG,Label:NIBR logP\",\"abs(subtract(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['nFO-P99OG_Ct6D7yJcZR2MOwPGA'].nValue),factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logp'].numericPrediction)),Label:Abs(Direct LogP - NIBR logP)\",\"factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logd'].numericPrediction,Aggregator:AVG,Label:NIBR logD\",\"abs(subtract(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['l38JNOlK9YR1XcHRQZgedBL074Y'].nValue),factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logd'].numericPrediction)),Label:Abs(Direct LogD - NIBR logD)\",\"factInstances[?(@.comparableFactSet.consolidationId == 'nJoy4ZmTecPINHzHxZFHEszWxxI')].compositeFieldInstances['pampa'].classPrediction,Label:pH6.8 HT Solubility Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'fRqZ0aL_lPZkWoyuzU9ALT_yltg')].compositeFieldInstances['pampa'].classPrediction,Label:PAMPA Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == '_8fg12KqRyqttsahD_fNxH1RX1k')].compositeFieldInstances['le_mdck'].classPrediction,Label:LE-MDCK Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'TZBnOQlYEL-5sz6GBe7t9eFBTUU')].compositeFieldInstances['rlm_clint'].classPrediction,Label:RLM CLint Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'TqzGxU_xWQiQ_dq3SWGF0bJzC1g')].compositeFieldInstances['hlm_clint'].classPrediction,Label:HLM CLint Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'Xugvrv0hpo5wqsZfE_g4yfQM96s')].compositeFieldInstances['herg'].classPrediction,Label:hERG Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castBoardName'].tValue,Label:Board Name\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaName'].tValue,Label:Idea Name\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaStatus'].tValue,Label:Idea Status\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaAuthor'].tValue,Label:Idea Author\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaSynthesizer'].tValue,Label:Idea Synthesizer\",\"factInstances[?(@.comparableFactSet.consolidationId == 'aXP-s3siyN7ouxg3aKjCeClHx_Y')].compositeFieldInstances['castIdeaCommentText'].tValue,Label:Idea Comment\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaNVP'].tValue,Label:Idea NVP\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaTags'].tValue,Label:Idea Tags\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaCsfSeries'].tValue,Label:Idea CFS Series\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaAvgRating'].nValue,Label:Idea Priority\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaCreated'].dValue,Label:Idea Created\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaLastActive'].dValue,Label:Idea Last Active\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castBoardId'].tValue,Label:Board ID\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaId'].tValue,Label:Idea ID\"]}";
    const strArraySlectors = "[\"structureRef,Label:Structure\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaId'].tValue,Label:CAST Idea ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'EJ3JyFsegb3dEVKbX4NwVSlgCPc')].compositeFieldInstances['9KsNABsTiq_smUmzLavqSFSnvbc'].activity,Aggregator:AVG,Label:IL17AA LZV623\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'EJ3JyFsegb3dEVKbX4NwVSlgCPc')].compositeFieldInstances['3XX9i9bjaIvD_CBDqikn9-lbudI'].dValue,Aggregator:MAX,Label:IL17AA LZV623 exp date\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'EJ3JyFsegb3dEVKbX4NwVSlgCPc')].compositeFieldInstances['X66Jkio5VXbEftpci6LA7LY2gjo'].imageId,Label:IL17AA LZV623 curve\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZoUwo87Ne896wGBnOb_XE-gtq7k')].compositeFieldInstances['6XpYwIht4MoHhd5OccXW7KTIVGE'].activity,Aggregator:AVG,Label:IL17AA IL17RA v2\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZoUwo87Ne896wGBnOb_XE-gtq7k')].compositeFieldInstances['qJo7i9V7-shoZylhm3vIqIQigjU'].imageId,Label:IL17AA IL17RA v2 curve\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'AIBBf7UwXJSvlsGB69HM2Dy78zU')].compositeFieldInstances['Ln82S1KvDZAiUqp5080krglYndg'].activity,Aggregator:AVG,Label:IL17AF IL17RA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ghCxfuDXaFA8RRlCVNVW8XiQFPo')].compositeFieldInstances['GGgFRl54EKEGLqB-wUCDLJXW0kw'].activity,Aggregator:AVG,Label:IL17 cell IL6\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'jbttmnOwP7eaHm64JvXe8RkChhM')].compositeFieldInstances['oPRlls9CnNzsKrDCge0Ai_YcbgE'].activity,Aggregator:AVG,Label:IL6 cell PP shift IL6\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rpyQkQjkGFF5Nc1XdmgymOCQJXw')].compositeFieldInstances['Sw-_YO8ZvRrn5WEpabhNPa_3MlQ'].activity,Aggregator:AVG,Label:IL17 cell rat\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'BcLJqvuV5zAxvRw8SQd3gNHpGgY')].compositeFieldInstances['_yWmg6nSl9bYOKE9QMLtxgF3D7k'].nValue,Aggregator:AVG,Label:IL17A SPR KD\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'BcLJqvuV5zAxvRw8SQd3gNHpGgY')].compositeFieldInstances['aLgA4eidTm1xvqc78KnQ-v9VvQE'].nValue,Aggregator:AVG,Label:IL17A SPR Kon\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'BcLJqvuV5zAxvRw8SQd3gNHpGgY')].compositeFieldInstances['ZN5gRV2H6vKHvPrJwHnFX2N7xYM'].nValue,Aggregator:AVG,Label:IL17A SPR Koff\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'x-OoKi2Vow7clPR2ES-KinMqsm8')].compositeFieldInstances['smfAvailability'].amount,Label:SMF Availability Availability Amount\",\"smallMoleculeSampleRefs.labhead,Label:Lab Head 521\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '7sTLZPVAcxsj7wzjnpNn8WKCEZo')].compositeFieldInstances['DJAy7t-40IIlp68XSfdM4wa8TW8'].tValue,Label:MAP in vivo PK DOG Route Admin\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '7sTLZPVAcxsj7wzjnpNn8WKCEZo')].compositeFieldInstances['K7mar7oX4FgddC0YfyfLwsznjwc'].nValue,Label:MAP in vivo PK DOG Dose\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '7sTLZPVAcxsj7wzjnpNn8WKCEZo')].compositeFieldInstances['H9qhyrc44YnwGNVhy6I3jyw1q2Y'].nValue,Label:MAP in vivo PK DOG CL\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['PZgPepgGKGP8YIoA4S9xTgjCRpM'].nValue,Aggregator:AVG,Label:Passive Permeability with MDCK-LE V2 PappAtoB\",\"smallMoleculeSampleRefs.preferredName,Label:Sample ID\",\"smallMoleculeSampleRefs.registrationDate,Label:Registration Date\",\"smallMoleculeSampleRefs.registrationDate,Aggregator:MIN,Label:First Registered\",\"smallMoleculeSampleRefs.chemist,Aggregator:MOST_COMMON,Label:Chemist\",\"smallMoleculeSampleRefs.labNotebookId,Label:Lab Notebook\",\"stereoCategory,Label:Stereo Category\",\"stereoDescriptor,Label:Stereo Descriptor\",\"smallMoleculeSampleRefs.projectCode,Aggregator:UNIQUE,Label:Project Code\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yo_8k8A73kiguu9TfV8SAqx3xn4')].compositeFieldInstances['Rah7fnilYoxxSfgrSBt5IJh21G0'].nValue,Aggregator:GEO,Label:HT Solubility pH6.8 [mM]\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yo_8k8A73kiguu9TfV8SAqx3xn4')].compositeFieldInstances['gK8vUaZLRgOwoTYA9SgyYWMjYHY'].nValue,Aggregator:GEO,Label:FaSSIF\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yo_8k8A73kiguu9TfV8SAqx3xn4')].compositeFieldInstances['3xfR8KjMIPd_8iKGHo3Jmsqtw_s'].dValue,Aggregator:MAX,Label:HT Solubility Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Sylv79kBc3leD0U-cR3A_63NuEc')].compositeFieldInstances['3TR05dUgcyFmyVvWmqzVYVjxLS0'].nValue,Aggregator:GEO,Label:SF-Solubility pH6.8 [mM]\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Sylv79kBc3leD0U-cR3A_63NuEc')].compositeFieldInstances['KIblChoT9LSMupVkStE0zpvPmgY'].tValue,Aggregator:MOST_COMMON,Label:SF-Solubility Center\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Sylv79kBc3leD0U-cR3A_63NuEc')].compositeFieldInstances['iOWQ8xL5eQdJ7dlMa69oUVnx4JM'].dValue,Aggregator:MAX,Label:SF-Solubility Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rV32B58AH_LPdsULrc435vhynKA')].compositeFieldInstances['0hvfnRo8FbmiHkqu_DHrR2WE7h8'].nValue,Aggregator:AVG,Label:logPAMPA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rV32B58AH_LPdsULrc435vhynKA')].compositeFieldInstances['vV7V9E8Ok6KN215nf6UGwebrOCA'].dValue,Aggregator:MAX,Label:logPAMPA Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['qkZ_DQyzhsRf3ECRCPOVBI9_74M'].nValue,Aggregator:GEO,Label:LE-MDCK PappAtoB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['bBVic8BCg-xxYUflnYTZMFklCm0'].activity,Aggregator:AVG,Label:LE-MDCK Recovery\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['u650cjE9QpyHLkRQkRzxXyezDic'].tValue,Aggregator:MOST_COMMON,Label:LE-MDCK FA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ORimVhasUcpP1-KytXn2ZMifzuo')].compositeFieldInstances['1SpmIY9A9F5SA6lV2g0aYUCFzbE'].dValue,Aggregator:MAX,Label:LE-MDCK Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['PZgPepgGKGP8YIoA4S9xTgjCRpM'].nValue,Aggregator:GEO,Label:LE-MDCK v2 PappAtoB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['x6ZLFG75FELMp41BRumoi_WyqjQ'].activity,Aggregator:AVG,Label:LE-MDCK v2  Recovery\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['hRca1-gc2nwBGAaKGvGo9ZoHm1c'].activity,Aggregator:AVG,Label:LE-MDCK v2 FA\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'D5nIm-c7Y1Rs1jq4esw_9vhpy24')].compositeFieldInstances['bS76DvBdPdICdE1dda211D02b60'].dValue,Aggregator:MAX,Label:LE-MDCK v2 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['A--0FCuOBBnNz0D5iIxJLPHj3PE'].nValue,Aggregator:GEO,Label:MDCK-MDR1 ER\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['BpfvavyJXK0obRR_3kR07dCXHgU'].nValue,Aggregator:AVG,Label:MDCK-MDR1 A-B\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['jU0XfZmFjh1Mb7GaqbOHUcJXsg8'].nValue,Aggregator:AVG,Label:MDCK-MDR1 B-A\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'CgFKzqPy0ed6ibJZdURmQ99QFDc')].compositeFieldInstances['Sni9_j1Qo0HrQcLGfKTbvctwGbo'].dValue,Aggregator:MAX,Label:MDCK-MDR1 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'cGEAi42quRk6vF554y41YFSfQVw')].compositeFieldInstances['MzZ3k5P-wW0cYMLFsDvcMJBEPVw'].tValue,Aggregator:MOST_COMMON,Label:OATP1B1 Substrate Assay\",\"length(most_common(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'cGEAi42quRk6vF554y41YFSfQVw')].compositeFieldInstances['MzZ3k5P-wW0cYMLFsDvcMJBEPVw'].tValue)),Label:OATP1B1 Substrate\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'cGEAi42quRk6vF554y41YFSfQVw')].compositeFieldInstances['HQAB_B4nA-VOh-lZe4cOVkY2xhU'].dValue,Aggregator:MAX,Label:OATP1B1 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'd0KLAfcOM8Vr5yzDaVSA-obfmi0')].compositeFieldInstances['s-ySb1q6vPopbIyEWnDRDM03JuY'].nValue,Aggregator:AVG,Label:HT LogP\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'd0KLAfcOM8Vr5yzDaVSA-obfmi0')].compositeFieldInstances['0YiXvrJFA9Lg2P0SbfdGXIFFKsk'].dValue,Aggregator:MAX,Label:HT LogP Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['nFO-P99OG_Ct6D7yJcZR2MOwPGA'].nValue,Aggregator:AVG,Label:Direct LogP\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['l38JNOlK9YR1XcHRQZgedBL074Y'].nValue,Aggregator:AVG,Label:Direct LogD\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['1wlD76Y3geUcXarB36a8rTWHYeM'].dValue,Aggregator:MAX,Label:Direct LogPD Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['B1Kx8dAmAESo-Wtc69gWCxUCp98'].nValue,Aggregator:MAX,Label:pKa1\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['yqGqf8hPiyf3M4QUwRXg2TIWiGI'].nValue,Aggregator:MAX,Label:pKa2\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['T2o8Lm5MtXZRvkvhSKKG6yHWV6U'].nValue,Aggregator:MAX,Label:pKa3\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['3_5biEAay1ekL8C9Rpf4wgfMpYM'].nValue,Aggregator:MAX,Label:pKa4\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['6ORyt6EWxGIpdHozCN6e0uPq12o'].tValue,Aggregator:MOST_COMMON,Label:pKa1_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['s2Q4yAR7ipaCO-8BR6ZOgpOa0jE'].tValue,Aggregator:MOST_COMMON,Label:pKa2_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['cdBoHQs26Hrzt3Hi7poPXh7VMKE'].tValue,Aggregator:MOST_COMMON,Label:pKa3_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LJONnCmRKBH0nZsJGN-1XsWD3lI')].compositeFieldInstances['eTtiNX_ZYRIQSG2O1ZDPZazqqfE'].tValue,Aggregator:MOST_COMMON,Label:pKa4_type\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'SW0kTdONJDEcF_KnycpS2J-18yI')].compositeFieldInstances['Zs46VVlLT4ziyZdQEeQgfi2Sew4'].nValue,Aggregator:GEO,Label:rLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'SW0kTdONJDEcF_KnycpS2J-18yI')].compositeFieldInstances['mHLIaq5yQGqSfi-rkC7ByMdeLD8'].dValue,Aggregator:MAX,Label:rLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'kQkCHBBEasvjH8Y8xNVI8dAQNnA')].compositeFieldInstances['QC2COEqmw2WVQdgtNhTNhL4U4Dw'].nValue,Aggregator:GEO,Label:mLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'kQkCHBBEasvjH8Y8xNVI8dAQNnA')].compositeFieldInstances['nyUdZpolM9dQD75Q2TiZ2Vncc58'].dValue,Aggregator:MAX,Label:mLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '-vddsB-WtoW7jlwQQTBPTGZgajs')].compositeFieldInstances['bMYHFdNUbjUgyc8p5bSpI4LjiOE'].nValue,Aggregator:GEO,Label:hLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '-vddsB-WtoW7jlwQQTBPTGZgajs')].compositeFieldInstances['AtpVqbOtagB8vxmy8JWK89ALDHs'].dValue,Aggregator:MAX,Label:hLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'hH4hg1T_m42_ZCJDq8mqXe20O0E')].compositeFieldInstances['ezEHQQMyhsECazauG8LtU77Dt0I'].nValue,Aggregator:GEO,Label:dLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'hH4hg1T_m42_ZCJDq8mqXe20O0E')].compositeFieldInstances['BpJmlYK02pFrIjhbZ1J86yETD7o'].dValue,Aggregator:MAX,Label:dLM Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Gq5S0jGMIhj7dGpdRcrDvb6giII')].compositeFieldInstances['yfLuBonQy0Ny8k25UEsd3Tul6MQ'].nValue,Aggregator:GEO,Label:cynoLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'e3ZPOcb6VhVEWZqtHhAEppTYnj0')].compositeFieldInstances['MhNRI6-NppT1nfK--Vcve_umwjA'].nValue,Aggregator:GEO,Label:minipigLM CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'm6kM9DPXPpgZs5lHgcBPlVUnjm4')].compositeFieldInstances['OPoAXu6ckjvDTF8JdqlrpQan1co'].nValue,Aggregator:GEO,Label:rHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'LMl5PyWaQXt8CW1_7DrzTKV357E')].compositeFieldInstances['XiODJf-kNWqB4QgJjWlO2LMUWYU'].nValue,Aggregator:GEO,Label:mHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Ed-roIcf615pfStMEHvW6n903_M')].compositeFieldInstances['Yf6hOGHcEgnPyyU2YOXln4zzbrY'].nValue,Aggregator:GEO,Label:hHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'Yv_0ZfUg0PkZvmOTaaSKnQmQX4M')].compositeFieldInstances['dTaypxGOk-4qkFsGzz5s0ra3Uqk'].nValue,Aggregator:GEO,Label:dHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'nrfwl7yWsZ6XMkSVVH-_kgxRc3Y')].compositeFieldInstances['wHzQhi_zYRJ-vDMHMicKN_fAVeI'].nValue,Aggregator:GEO,Label:cynoHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'pxgyV6XiF6sKzutBPnAtWYMYoKQ')].compositeFieldInstances['zwfMDuceMoZkWpfavBfC_610sPI'].nValue,Aggregator:GEO,Label:minipigHEP CLint\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'mVivAeeHG0gQSlqXK8XU6FzJFxU')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:rPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '0t8LmVLFGFeq0iKTtUOxslgglGI')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:mPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'nWPGpBxxnEXWP2MMeo4sLbiG158')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:hPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'sTfm8emCiu6f92PfW5NjF7qY-Cs')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:dPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'rV5INo9qP628SaP9FBPcN8GqGg8')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:cynoPPB PB\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 't8kq0fsl7JA2OUguGwr144LuDMk')].compositeFieldInstances['imdQIUrX4hEs7w-R2_jvQdZQolk'].nValue,Aggregator:GEO,Label:minipigPPB PB\",\"divide(remove_qual(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'eRvF0tyNCvVnvCLeUArIn2SOKaU')].compositeFieldInstances['iNK51tXhvRVP1M-Lys_0hqviAyQ'].nValue)),100),Label:rfu(mic)\",\"divide(remove_qual(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'sXzxu6x4_DxkM27u7Yd4yXyLHes')].compositeFieldInstances['iNK51tXhvRVP1M-Lys_0hqviAyQ'].nValue)),100),Label:hfu(mic)\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['dL_0aWPbTKonqnGqZQvOHWw26yk'].nValue,Aggregator:GEO,Label:CYP3A4 TDI kObs\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['_ib14_aeIFm_ByUaaLTsGLWVyEA'].tValue,Aggregator:MOST_COMMON,Label:CYP3A4 TDI Inh-Rev\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['AXW4CGHn41mJg5nfaJgWVZiuEVM'].tValue,Aggregator:MOST_COMMON,Label:CYP3A4 TDI Flag\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'xI-aJodsMENg41O_ffqRliSCWYQ')].compositeFieldInstances['7GamtjRP9RigQ6oBFjZMC9IcxhI'].dValue,Aggregator:MAX,Label:CYP3A4 TDI Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zCjqunJldyDxOKs6xxOEzchLOO8')].compositeFieldInstances['NORxqb6pWHCmLHs62f6Q57ChKGE'].activity,Aggregator:GEO,Label:CYP3A4 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zCjqunJldyDxOKs6xxOEzchLOO8')].compositeFieldInstances['5LnswVHMIVO98Q1ZTG7dcPzdiA0'].dValue,Aggregator:MAX,Label:CYP3A4 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'E3UyRf8YXKf_ysiAjo6KImmtsaE')].compositeFieldInstances['vpyt3dghj1VQ9thRte6GcHwG0MY'].activity,Aggregator:GEO,Label:CYP2D6 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'l-zKqbBDyJv2Yoian60ZOVRDVvc')].compositeFieldInstances['rIBFbovjkIUZq_SOXFSzbREW-0g'].activity,Aggregator:GEO,Label:CYP2C9 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['izC7ons3rj0jpiiW_Q-ptyt7D4w'].activity,Aggregator:GEO,Label:hERG Binding IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['oTttI7IzFIAabLCS8NS425jHTJg'].activity,Aggregator:AVG,Label:hERG Binding inh@1uM\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['S74KXMvFXsal2nN8GNZifKeQCsM'].activity,Aggregator:AVG,Label:hERG Binding inhib@30uM\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['izC7ons3rj0jpiiW_Q-ptyt7D4w'].imageId,Label:hERG Binding Curve ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'zw82I12TukMCA5usE7BaHRkpiXE')].compositeFieldInstances['T4gbZBctaeUMKgLYIK0Zrpko_z4'].dValue,Aggregator:MAX,Label:hERG Binding Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'TTyaVjvh-yhQDBlokdjq_DVymtU')].compositeFieldInstances['1y5Lk8vxAlXAGJCd0cEVszsRGpw'].activity,Aggregator:GEO,Label:hERG QPatch IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'TTyaVjvh-yhQDBlokdjq_DVymtU')].compositeFieldInstances['1y5Lk8vxAlXAGJCd0cEVszsRGpw'].hill,Aggregator:MIN,Label:hERG QPatch Hill slope\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'TTyaVjvh-yhQDBlokdjq_DVymtU')].compositeFieldInstances['YSI1gn-QJHiaGU-MdfHhr37axyM'].dValue,Aggregator:MAX,Label:hERG QPatch Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZGU-0ne7HW7_9_DdJZxlL8zGc-o')].compositeFieldInstances['ODcF86fJwqifZ9ChnTXCNQHeSWQ'].activity,Aggregator:GEO,Label:BSEP IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ZGU-0ne7HW7_9_DdJZxlL8zGc-o')].compositeFieldInstances['Au5cnkk1XysekBtjLZ2IuoKrcRY'].dValue,Aggregator:MAX,Label:BSEP Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '9YZ73Ke-sVcQuT29Wg1X9Stw1Rs')].compositeFieldInstances['JSj7zuEnJ2nrXCWdDu-533QqWaQ'].activity,Aggregator:GEO,Label:Nav1.5 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '9YZ73Ke-sVcQuT29Wg1X9Stw1Rs')].compositeFieldInstances['JSj7zuEnJ2nrXCWdDu-533QqWaQ'].imageId,Label:Nav1.5 Curve ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == '9YZ73Ke-sVcQuT29Wg1X9Stw1Rs')].compositeFieldInstances['Mx2kdAbsw22_o-tjYbeujYiW1lc'].dValue,Aggregator:MIN,Label:Nav1.5 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'RI0c7xNdC0J2askGjmQeJuDoUHc')].compositeFieldInstances['mi4-u0Oq3M6FWOTx_MtCNSK56qg'].activity,Aggregator:GEO,Label:Cav1.2 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'RI0c7xNdC0J2askGjmQeJuDoUHc')].compositeFieldInstances['OqzK0_SCHVo4sGpkJiLa1vg5OZA'].dValue,Aggregator:MAX,Label:Cav1.2 Date-out\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ze8J1qI_fscCkAewFBUUjXgJjxU')].compositeFieldInstances['WZJY9Rn7klxeEkUpTzmk4sPBZCc'].activity,Aggregator:GEO,Label:KCNQ1 IC50\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ze8J1qI_fscCkAewFBUUjXgJjxU')].compositeFieldInstances['WZJY9Rn7klxeEkUpTzmk4sPBZCc'].imageId,Label:KCNQ1 Curve ID\",\"smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'ze8J1qI_fscCkAewFBUUjXgJjxU')].compositeFieldInstances['GFJ_lAiuTRTlpzi8aZQNBeMt9DA'].dValue,Aggregator:MAX,Label:KCNQ1 Date-out\",\"averageMass,Label:Average Mass\",\"tpsa,Label:TPSA\",\"nRotBonds,Label:Num Rotatable Bonds\",\"nHeavyAtoms,Label:Num Heavy Atoms\",\"factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logp'].numericPrediction,Aggregator:AVG,Label:NIBR logP\",\"abs(subtract(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['nFO-P99OG_Ct6D7yJcZR2MOwPGA'].nValue),factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logp'].numericPrediction)),Label:Abs(Direct LogP - NIBR logP)\",\"factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logd'].numericPrediction,Aggregator:AVG,Label:NIBR logD\",\"abs(subtract(avg(smallMoleculeSampleRefs.factInstances[?(@.comparableFactSet.consolidationId == 'MFpNOKKTRL2OdBdypvRONh3rT1s')].compositeFieldInstances['l38JNOlK9YR1XcHRQZgedBL074Y'].nValue),factInstances[?(@.comparableFactSet.consolidationId == 'Vv9y0eBDvMets8gTkx26EAi-YKY')].compositeFieldInstances['nibr_logd'].numericPrediction)),Label:Abs(Direct LogD - NIBR logD)\",\"factInstances[?(@.comparableFactSet.consolidationId == 'nJoy4ZmTecPINHzHxZFHEszWxxI')].compositeFieldInstances['pampa'].classPrediction,Label:pH6.8 HT Solubility Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'fRqZ0aL_lPZkWoyuzU9ALT_yltg')].compositeFieldInstances['pampa'].classPrediction,Label:PAMPA Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == '_8fg12KqRyqttsahD_fNxH1RX1k')].compositeFieldInstances['le_mdck'].classPrediction,Label:LE-MDCK Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'TZBnOQlYEL-5sz6GBe7t9eFBTUU')].compositeFieldInstances['rlm_clint'].classPrediction,Label:RLM CLint Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'TqzGxU_xWQiQ_dq3SWGF0bJzC1g')].compositeFieldInstances['hlm_clint'].classPrediction,Label:HLM CLint Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'Xugvrv0hpo5wqsZfE_g4yfQM96s')].compositeFieldInstances['herg'].classPrediction,Label:hERG Classification\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castBoardName'].tValue,Label:Board Name\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaName'].tValue,Label:Idea Name\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaStatus'].tValue,Label:Idea Status\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaAuthor'].tValue,Label:Idea Author\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaSynthesizer'].tValue,Label:Idea Synthesizer\",\"factInstances[?(@.comparableFactSet.consolidationId == 'aXP-s3siyN7ouxg3aKjCeClHx_Y')].compositeFieldInstances['castIdeaCommentText'].tValue,Label:Idea Comment\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaNVP'].tValue,Label:Idea NVP\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaTags'].tValue,Label:Idea Tags\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaCsfSeries'].tValue,Label:Idea CFS Series\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaAvgRating'].nValue,Label:Idea Priority\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaCreated'].dValue,Label:Idea Created\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaLastActive'].dValue,Label:Idea Last Active\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castBoardId'].tValue,Label:Board ID\",\"factInstances[?(@.comparableFactSet.consolidationId == 'AZMmTED35E_ygpCklyc8QKxvdTc')].compositeFieldInstances['castIdeaId'].tValue,Label:Idea ID\"]";
    const arArraySlectors  = JSON.parse(strArraySlectors);
   // let arDescriptors = await MetaDataLoader.load(arArraySlectors);

}


//tags: appp
//name: MolService Demo
export async function molSeiviceDemo(s)
{
    const nMolCount = 100000;
    const arSmiles = new Array(nMolCount);
    ChemUtils.fillRandomSmiles(arSmiles);

    let bCancel = false;
    const dialog = new ProgressDialog("MolService", true, () => {
     bCancel = true;
    });

    await dialog.show();
    await dialog.setUpperText("Creating Molecules...");

    let service = null;
    try{service = await MolServiceNew.createMolService(RDKitEngineNew, arSmiles, false, async (nPrgMolCount, argsCancel) =>
    {
        const nPct = Math.floor(nPrgMolCount/nMolCount*100);
        if(bCancel)//nPct > 50){
            argsCancel.cancel();
        //}

        await dialog.setProgress(nPct);
        await dialog.setBottomText(nPrgMolCount + " out of " + nMolCount);
    });}
    catch(e)
    {
        let fghfhh= 0;
        throw e;
    }

    //await dialog.setProgress(100);
    await dialog.setUpperText(bCancel ? "Cancelled" : "Success...");
    await dialog.setBottomText(bCancel ? "Process Stopped" : "Loading Complete");
    //////////////await dialog.close();

    /*
const nRecordCount = 10000;
const dframe = await GNFUtils.generateRandomDataFrame(nRecordCount, _package);
const nCol = GridUtils.findColumnBySemType(dframe, CpdSemType);
const colCpd = dframe.columns.byIndex(nCol);
const service = colCpd.getTag("MOL_SERVICE");*/
    const bitmap =  await service.createMolImage(1, 200, 200);

    const eCanvasCpd = ui.canvas();
    eCanvasCpd.id     = "idDrawAreaCpd";
    eCanvasCpd.width  = 200;
    eCanvasCpd.height = 200;
    let dial = ui.dialog("Canvas");
    dial.add(eCanvasCpd);
    try{dial.show(true);}
    catch(e)
    {
        throw e;
    }
    const g = eCanvasCpd.getContext("2d");

    /*
  const gradient = g.createLinearGradient(0,0, 0,200);
  gradient.addColorStop(0, 'white');
  gradient.addColorStop(.5, 'LightGray');
  gradient.addColorStop(1, 'white');
  g.fillStyle = gradient;
  g.fillRect(0, 0, 200, 200);
     */


    g.drawImage(bitmap, 0,0,200,200);

    await service.dispose();

    //await DGApp.open(_package, dframe);
}


//tags: appp
//name: MolServiceNew Demo
export async function molSeiviceDemoNew(s)
{
    const nMolCount = 100000;
    const arSmiles = new Array(nMolCount);
    ChemUtils.fillRandomSmiles(arSmiles);


/*
    let service = null;
    try{service = await MolServiceNew.createMolServiceNew(RDKitEngineNew, arSmiles, false, async (nPrgMolCount, argsCancel) =>
    {
        const nPct = Math.floor(nPrgMolCount/nMolCount*100);
        if(bCancel)//nPct > 50){
            argsCancel.cancel();
        //}

        await dialog.setProgress(nPct);
        await dialog.setBottomText(nPrgMolCount + " out of " + nMolCount);
    });}
    catch(e)
    {
       throw e;
    }
*/

    async function OnProgress(nPrgMolCount)
    {
        const nPct = Math.floor(nPrgMolCount/nMolCount*100);
        await dialog.setProgress(nPct);
        await dialog.setBottomText(nPrgMolCount + " out of " + nMolCount);
    }


    const service = new MolServiceNew(RDKitEngineNew);
    let task = service.addMolecules(arSmiles, false);
    const subscrPrg = task.onUpdate.subscribe(OnProgress);// = OnProgress;
    const subscrFin = task.onFinished.subscribe((obResult) => {

        let aaa = obResult;
    });


    let bCancel = false;
    const dialog = new ProgressDialog("MolService", true, () => {
        bCancel = true;
        task.cancel();
    });

    await dialog.show();
    await dialog.setUpperText("Creating Molecules...");

    task.run();

    await task.awaitResult();
    const nProcMolCount = task.getResult();
    const status = task.getStatus();

    //await dialog.setProgress(100);
    await dialog.setUpperText(status.getTitle() + "...");
    await dialog.setBottomText(status !== TaskStatus.RanSuccesfully ? "Process Stopped" : "Loading Complete");
    //////////////await dialog.close();

    subscrPrg.unsubscribe();


    task =  service.createMolImage(1, 200, 200);
    task.run();

    await task.awaitResult();
    const bitmap = task.getResult();




    const eCanvasCpd = ui.canvas();
    eCanvasCpd.id     = "idDrawAreaCpd";
    eCanvasCpd.width  = 200;
    eCanvasCpd.height = 200;
    let dial = ui.dialog("Canvas");
    dial.add(eCanvasCpd);
    try{dial.show(true);}
    catch(e)
    {
        throw e;
    }
    const g = eCanvasCpd.getContext("2d");
    g.drawImage(bitmap, 0,0,200,200);


    //Structure Search
    const arFlags = new Array(nMolCount);
    task =  service.structureSearch(arFlags, arSmiles[2]);
    task.run();

    await task.awaitResult();
    const aResFlags = task.getResult();


    //Tanimoto Scores
    const arScores = new Array(nMolCount);
    task =  service.tanimotoScores(arScores, arSmiles[2]);
    task.run();

    await task.awaitResult();
    const aResScoress = task.getResult();




    task = service.dispose();
    task.run();

    await task.awaitResult();
    let aaa = 0;
}



//tags: appp
//name: Progress Indicator Demo
export async function progressDemo(s)
{
    let dialog = new ProgressDialog("Progress Indicator...");
    await dialog.show();
    await dialog.setUpperText("Please Wait...");
    let eCanvas = dialog.getCanvas();

    let nProgress = 0;
    eCanvas.addEventListener("click",async function(e){
        nProgress += 5;
        await dialog.setProgress(nProgress);
        await dialog.setBottomText("Setting Super Long Long Bottom Text on " + nProgress);
    });
}




