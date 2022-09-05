import * as DG from "datagrok-api/dg";
import {SemType} from "../entity/SemType";
import {CpdSemType} from "../entity/gnf/cpd/CpdSemType";
import {CpdEntity} from "../entity/gnf/cpd/CpdEntity";
import {DGApp} from "../entity/app/DGApp";
import {GridUtils} from "../utils/GridUtils";
import {RDKitEngine} from "../chem/rdkit/RDKitEngine";
import {AnalysisLoader} from "../service/nx/analysis/AnalysisLoader";
import {CDFEntities} from "../service/nx/cdf/entity/CDFEntities";
import {DataSetTagGrouping} from "./DataSetTagGrouping";
import {MolService} from "../chem/MolService";
import {CpdPrimitiveSemType} from "../entity/gnf/cpd/CpdPrimitiveSemType";
import {MathUtils} from "../utils/MathUtils";
import {FiltersUtils} from "../ui/filters/FiltersUtils";
import {FilterPanel} from "../ui/filters/FilterPanel";

export class IntuenceDiscovery {
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
 }

IntuenceDiscovery.isCpdColumn = function(obTagDescriptor, strColName) {

    if(obTagDescriptor !== null && obTagDescriptor !== undefined) {

        const typeParent = obTagDescriptor.getPrimitiveType().getParentType();
        const b = typeParent === CDFEntities.SmallMoleculeConcept || typeParent === CDFEntities.SmallMoleculeSample;
        return b;
    }

   const b = strColName === "Id" || strColName === "Sample ID" || strColName === "Structure" || strColName === "Last Experimental Date" || strColName === "Last Published Date";
    return b;
}

IntuenceDiscovery.processAnalysisData = async function(dframe, configSemTypes, obPackage, dialProgress) {

    const DEBUG = false;

    if (dialProgress !== undefined) {

        await dialProgress.setUpperText("Please Wait...");
        await dialProgress.setBottomText("Creating Entities...");
        await dialProgress.setProgress(0);
    }

    let bFinished = false;

    const KEY_CDF_CFS_ID = ".CDF_CFS_ID";

    let arCpdColNames = [];

    const mapAssayCols = new Map();
    let arColNames1 = null;

    let lstCols = null;
    let nColCount = -1;
    let col = null;
    let strTag = null;
    let strTagColId = null;
    let obTagDescriptor = null;
    let typeSem = null;

    lstCols = dframe.columns.toList();
    nColCount = lstCols.length;
    strTag = null;
    strTagColId = null;
    obTagDescriptor = null;
    for (let nCol = 0; nCol < nColCount; ++nCol) {
        col = lstCols[nCol];
        typeSem = GridUtils.assignSemType(col);
        typeSem.init(obPackage);
        //GridUtils.setSemType(col, PrimitiveSemType.Instance);

        strTag = col.getTag(KEY_CDF_CFS_ID);
        strTagColId = col.getTag(AnalysisLoader.KEY_TAG_SELECTOR);
        obTagDescriptor = col.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);

        if (IntuenceDiscovery.isCpdColumn(obTagDescriptor, col.name)) {
            arCpdColNames.push(col.name);
            GridUtils.setSemType(col, CpdPrimitiveSemType.Instance);
            continue;
        }

        console.log(nCol + " " + col.name + " " + strTag);
        if (strTag !== null && strTag !== undefined) {
            arColNames1 = mapAssayCols.get(strTag);
            if (arColNames1 === null || arColNames1 === undefined) {
                arColNames1 = [];
                mapAssayCols.set(strTag, arColNames1);
            }

            arColNames1.push(col.name);//{name: col.name, property: SemType.UnknownProp});
        }
    }


    if (arCpdColNames.length > 0) {
        console.log("Compound: ");
        console.log("\t" + arCpdColNames);

        const arCpdColsGroup = GridUtils.names2Cols(dframe, arCpdColNames);

        let colVirt;
        try {
            colVirt = dframe.columns.addNewVirtual("Compound_V", function (nRow) {

                if(!bFinished)
                    return null;

                const arTmpColNames = colVirt.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);

                let colPri = dframe.getCol("Id");
                const strConceptId = colPri.get(nRow);
                colPri = dframe.getCol("Sample ID");
                const strSampleId = colPri.get(nRow);
                colPri = dframe.getCol("Structure");
                const strSmiles = colPri.get(nRow);

                let nTime = NaN;

                if (arTmpColNames.includes("Last Experimental Date")) {
                    colPri = dframe.getCol("Last Experimental Date");
                    nTime = colPri.get(nRow);
                    if(MathUtils.isNullValue(nTime))
                        nTime = NaN;
                    else if (typeof ob !== "number")// && !MathUtils.isNullValue(nTime))
                        nTime = nTime.a;
                } else {

                    if (arTmpColNames.includes("Last Published Date")) {
                        colPri = dframe.getCol("Last Published Date");
                        nTime = colPri.get(nRow);
                        if(MathUtils.isNullValue(nTime))
                            nTime = NaN;
                        else if (typeof ob !== "number")// && !MathUtils.isNullValue(nTime))
                            nTime = nTime.a;
                   }
                }


                //colPri = dframe.getCol("Last Experimental Date");
                //let t = colPri.type;
                //let dt = colPri.get(nRow);
                //let nTime = dt === undefined || dt === null ? -1 : dt.a;
                //nTime = new Date().getTime();
                const cpd = new CpdEntity(strSmiles, true, strSampleId, strConceptId, nTime);

                return cpd;//arCpdColNames.toString();
            });
        } catch (e) {
            throw e;
        }

        const typeSemCpd = new CpdSemType();
        GridUtils.setSemType(colVirt, typeSemCpd);
        typeSemCpd.init(obPackage, arCpdColsGroup);
        colVirt.setTag(DGApp.PRIMITIVE_COLS_TAG_NAME, arCpdColNames);

        const colStructure = dframe.getCol("Structure");
        const nRecordCount = colStructure.length;
        const arMolFiles = new Array(nRecordCount);
        let strMolFile = null;
        for (var nR = 0; nR < nRecordCount; ++nR) {
            strMolFile = colStructure.get(nR);
            arMolFiles[nR] = strMolFile;
        }


        if (dialProgress !== undefined) {
            await dialProgress.setUpperText("Creating Molecules...");
            await dialProgress.setBottomText("0 out of " + nRecordCount);
        }
        let service = null;
        try {
            service = await MolService.createMolService(RDKitEngine, arMolFiles, true, async (nPrgMolCount) => {

                if (dialProgress === undefined)
                    return;

                const nPct = Math.floor(nPrgMolCount / nRecordCount * 100);

                await dialProgress.setProgress(nPct === 100 ? 99 : nPct);
                await dialProgress.setBottomText(nPrgMolCount + " out of " + nRecordCount);


            });
        }//Progress Handler
        catch (e) {
            throw e;
        }


        typeSemCpd.setMolService(service);

        if (dialProgress !== undefined) {
            await dialProgress.setUpperText("Success...");
            await dialProgress.setBottomText("Performing Final Steps... ");
            await dialProgress.setProgress(99);
        }
    }


    if (dialProgress !== undefined) {
        await dialProgress.setUpperText("Creating Virtual Columns...");
    }

    let itKeys = mapAssayCols.keys();
    let nKey = 0;
    let strColName = null;
    let arColNames = null;
    let arColsGroup = null;

    for (let strKey of itKeys) {

        if((nKey % 10) === 0 && dialProgress !== undefined)
            await dialProgress.setBottomText((nKey + 1) + " out of " + mapAssayCols.size);

        arColNames = mapAssayCols.get(strKey);
        if(arColNames === null || arColNames === undefined)
            throw Error("Array of columns cannot be null or undefined for key:  " + strKey);

        arColsGroup = GridUtils.names2Cols(dframe, arColNames);

        try {typeSem = configSemTypes.findSemType(arColsGroup);}
        catch (e) {
            configSemTypes.findSemType(arColsGroup);
            throw e;
        }

        if (typeSem === null) {//left  for debufgging
            for (var nC = 0; nC < arColsGroup.length; ++nC) {
                col = arColsGroup[nC];
                obTagDescriptor = col.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);
                let typeCDFPrimitive = obTagDescriptor.getPrimitiveType();
                let typeCol = col.type;
                let tagGrouuping = col.getTag(DataSetTagGrouping.KEY);
                let ccc = 0;
            }
        }


        col = dframe.columns.byName(arColNames[0]);
        obTagDescriptor = col.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);
        if (obTagDescriptor !== null && obTagDescriptor !== undefined && obTagDescriptor.getTitle() !== null)
            strColName = obTagDescriptor.getTitle() + " (V)";
        else
            strColName = arColNames[0] + "_V";

        if (strColName === null)
            console.log("Null Assay Title");

        let colExisting = dframe.columns.byName(strColName);
        if(colExisting !== null)
        {
            let nVersion = 1;
            for(;nVersion<1000; ++nVersion)
            {
                colExisting = dframe.columns.byName(strColName + " " + nVersion.toString());
                if(colExisting === null)
                  break;
            }

            strColName += (" " + nVersion.toString());
        }



        if(DEBUG) {
            console.log("Tag: " + strKey);
            console.log("\t" + strColName);
            console.log("\t PRI COLS: " + arColNames);
        }

        try {
            let colVirt;
            colVirt = dframe.columns.addNewVirtual(strColName, function (nRow) {

                if(!bFinished || nRow === undefined || nRow === null)//bug in DG
                    return null;

                const typeSemVirt = GridUtils.getSemType(colVirt);
                if (typeSemVirt instanceof SemType)
                {
                    const arTmpColNames = colVirt.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
                    if(arTmpColNames === null || arTmpColNames === undefined)
                        throw Error("Array of columns cannot be null or undefined for virt col " + colVirt.name);

                    const arTmpCols = GridUtils.names2Cols(dframe, arTmpColNames);
                    let entity = null;
                    try {entity = typeSemVirt.createEntity(arTmpCols, nRow);}
                    catch(e)
                    {

                        try {entity = typeSemVirt.createEntity(arTmpCols, nRow);}
                        catch(eee)
                        {

                        }

                       const msg = e.message;
                       throw e;
                    }
                    return entity;
                }
                //return colVirt.d.primitives.toString();
            });

            GridUtils.setSemType(colVirt, typeSem);
            await typeSem.init(obPackage, arColsGroup);
            colVirt.setTag(DGApp.PRIMITIVE_COLS_TAG_NAME, arColNames);
            const strAltColName = typeSem.createColName(arColsGroup);
            if(strAltColName !== null)
             colVirt.setTag(DGApp.VIRTUAL_COL_ALT_NAME_TAG_NAME, strAltColName);

        } catch (e) {

            if (dialProgress !== undefined) {
                await dialProgress.setUpperText("ERROR...");
                await dialProgress.setBottomText("Failure on Column  " + strColName);
                console.log("Failure on column: " + strColName)
            }
        }
        ++nKey;
     }

    bFinished = true;

    if(dialProgress !== undefined)
        await dialProgress.setProgress(99);
}

IntuenceDiscovery.open = async function(viewSS, configGridRenderers, arFilterClasses, obPackage, fnFilterOpener, fnFilterCloser) {

        const grid = viewSS.grid;

        let b = true;
       if(fnFilterCloser !== null && fnFilterCloser !== undefined) {
           try{fnFilterCloser();}
           catch(e)
           {
               b = false;
               console.log("ERROR: Failed to close Filters");
           }
       }
        await DGApp.open(obPackage, viewSS, configGridRenderers, arFilterClasses);

       //if(b) {
           if (fnFilterOpener !== null && fnFilterOpener !== undefined) {
               try{fnFilterOpener();}
               catch(e) {
                   console.log("ERROR: Failed to open Filters");
               }
           }
       //}

        //let viewGrid = grok.shell.addTableView(dframe);
        const listGridCols = grid.columns;
        const nGridColCount = GridUtils.getColumnCount(grid);
        const arColOrderNames = [];

        let colGrid = null;
        let nChildCount = -1;
        let colGridPri = null;
        let arPrimiColNames = null;
        let obTagDescriptor = null;
        for (var nGCol = 1; nGCol < nGridColCount; ++nGCol)
        {
            colGrid = listGridCols.byIndex(nGCol);
            if(GridUtils.getSemType(colGrid.column) instanceof SemType && DGApp.isVirtual(colGrid.column))
            {
                colGrid.visible = true;

                arColOrderNames.push(colGrid.column.name);

                arPrimiColNames = colGrid.column.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
                nChildCount = arPrimiColNames.length;
                for(var nPCol = 0; nPCol <nChildCount ; ++nPCol)
                {
                    colGridPri = listGridCols.byName(arPrimiColNames[nPCol]);
                    colGridPri.visible = false;

                    colGridPri.column.setTag(DGApp.VIRTUAL_PARENT_TAG_NAME, colGrid.column.name);
                    arColOrderNames.push(arPrimiColNames[nPCol]);


                    obTagDescriptor = colGridPri.column.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);
                    if (obTagDescriptor !== null && obTagDescriptor !== undefined) {
                      const typeCDFPrimitive = obTagDescriptor.getPrimitiveType();
                        colGridPri.column.setTag(DGApp.PRIMITIVE_COL_LABEL_TAG_NAME, typeCDFPrimitive.getKey());
                    }

                }
            }
        }

        try{listGridCols.setOrder(arColOrderNames);}
        catch(e)
        {
         let rrr = 0;
         }


    const nColCpd = GridUtils.findGridColumnBySemType(grid, CpdSemType);
    const colGridCpd = grid.columns.byIndex(nColCpd);  //my changes 1
    const strColNameCpd = colGridCpd.name;
   //my changes  colGcpd.visible = false;

    FiltersUtils.showFilters(false);

    const viewerFilters = FiltersUtils.getFiltersViewer();
    if(viewerFilters !== null) {
        const panelFilters = new FilterPanel(arFilterClasses);

        viewerFilters.root.appendChild(panelFilters.root);
        panelFilters.attach(viewerFilters.dataFrame);
        viewerFilters.dart.m_panelFilters = panelFilters;

        FiltersUtils.openFilter(colGridCpd.column);
    }
}
