import {CpdSemType} from "./cpd/CpdSemType";
import {CpdEntity} from "./cpd/CpdEntity";
import {InVivoEntity} from "./invivo/InVivoEntity";
import {InVivoSemType} from "./invivo/InVivoSemType";
import {InSilicoEntity} from "./insilico/InSilicoEntity";
import {InSilicoSemType} from "./insilico/InSilicoSemType";
import {IC50Entity} from "./ic50/IC50Entity";
import {IC50SemType} from "./ic50/IC50SemType";
import {RDKitEngine} from "../../chem/rdkit/RDKitEngine";
//import {OCLEngine} from "../../chem/ocl/OCLEngine";
import {ProgressDialog} from "../../ui/progress/ProgressDialog";
import * as ui from "datagrok-api/ui";
import {MolService} from "../../chem/MolService";
import {GridUtils} from "../../utils/GridUtils";
import {RandUtils} from "../../utils/RandUtils";

export class GNFUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}

GNFUtils.generateRandomDataFrame = async function(nRecordCount, obPackage) {

    if(nRecordCount === undefined)
        throw new Error("THe number of records cannot be undefined.");

    let nIC50ColCount = 50;
    let arCols = new Array(nIC50ColCount+3);
    arCols[0] = DG.Column.fromType(DG.TYPE.OBJECT, "Cpd", nRecordCount);
    const typeSemCpd = new CpdSemType();
    GridUtils.setSemType(arCols[0], typeSemCpd);

    arCols[1] = DG.Column.fromType(DG.TYPE.OBJECT, "InVivo", nRecordCount);//DG.TYPE.FLOAT, "InVivo", nRecordCount);
    GridUtils.setSemType(arCols[1], new InVivoSemType());
    arCols[2] = DG.Column.fromType(DG.TYPE.OBJECT, "InSilico", nRecordCount);
    GridUtils.setSemType(arCols[2], new InSilicoSemType());

    var nCol=0;
    for(; nCol<nIC50ColCount; ++nCol)
    {
        arCols[nCol+3] = DG.Column.fromType(DG.TYPE.OBJECT, "A super very long IC50 assay name imaginable " + nCol.toString(), nRecordCount);
        GridUtils.setSemType(arCols[nCol+3], IC50SemType.Instance);//new IC50SemType();
    }



    //colIC50.semType = "ic50";
    //colIC50.setTag("quality", "ic50");
    //let colFake = new ThisColumn(colF.d, "NameFake", 5);

    let strSmiles = "";
    let strSampleId = "";
    let nTime = -1;
    for(var n=0; n<nRecordCount; ++n)
    {
        strSmiles = CpdEntity.generateRandomSmiles();

        strSampleId = "NVP-ASF-" + RandUtils.generateRandomInteger(10000, 90000);
        nTime = RandUtils.generateRandomDate(2007, 2021);
        try{arCols[0].set(n, new CpdEntity(strSmiles, false, strSampleId, strSampleId, nTime));}
        catch(e)
        {
            let rdfg = 0;
        }

        arCols[1].set(n, new InVivoEntity());
        arCols[2].set(n, new InSilicoEntity());

        for(nCol=0; nCol<nIC50ColCount; ++nCol)
        {
            if(Math.random() >= 0.2)
                arCols[nCol + 3].set(n, new IC50Entity());
        }
    }

    // arCols[arCols.length -1] = DG.Column.fromType(DG.TYPE.FLOAT, "Primitive", nRecordCount);

    let arSmiles = new Array(nRecordCount);
    for(var nR=0; nR<nRecordCount; ++nR)
    {
        arSmiles[nR] = arCols[0].get(nR).getSmiles();
    }

    for(nCol=0; nCol<arCols.length; ++nCol)
    {
        await GridUtils.getSemType(arCols[nCol]).init(obPackage);
    }

    //Progress Dialog

    const dialog = new ProgressDialog("MolService", true);
    await dialog.show();
    await dialog.setUpperText("Creating Molecules...");

    let service = null;
    try{service = await MolService.createMolService(RDKitEngine, arSmiles, false, async (nPrgMolCount) =>
    {
        let nPct = Math.floor(nPrgMolCount/arSmiles.length*100);

        await dialog.setProgress(nPct);
        await dialog.setBottomText(nPrgMolCount + " out of " + arSmiles.length);

    });}
    catch(e)
    {
        throw e;
    }

    await dialog.setProgress(100);
    await dialog.setUpperText("Success...");
    await dialog.setBottomText("Loading Complete");
    await dialog.close();

    //await service.dispose();
    //let arFlags = new Array(arSmiles.length);
    //await service.similaritySearch(arFlags, arSmiles[0], 0.5);


    let dframe =  null;
    try{dframe = DG.DataFrame.fromColumns(arCols);}
    catch(e)
    {
        throw e;
    }

    typeSemCpd.setMolService(service);
    return dframe;
}
