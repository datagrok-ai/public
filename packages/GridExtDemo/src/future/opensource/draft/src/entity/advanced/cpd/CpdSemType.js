import {SemType} from "../../SemType";
import {SemSorter} from "../../SemSorter";
import {MathUtils} from "../../../utils/MathUtils";
import {GridUtils} from "../../../utils/GridUtils";

export class CpdSemType extends SemType {
    constructor() {
        super(undefined, [new CpdSemType.CpsLastTestedSorter(), new CpdSemType.SimilaritySorter(),
            new CpdSemType.GNFIDSorter(), new CpdSemType.SIDSorter(), new CpdSemType.ConceptIDSorter(),
            new CpdSemType.SampleIDSorter()]);

        this.m_serviceMol = null;
    }

    getIconFilePath() {return "images/cpd.png";}
    getMolService() {return this.m_serviceMol;}
    setMolService(service)
    {
        this.m_serviceMol = service;
    }
}


CpdSemType.CpsLastTestedSorter = class extends SemSorter
{
    getName() {return "Sort by Last Tested";}

    isNullSortValue(ob)
    {
        if(super.isNullSortValue(ob))
            return true;

        const nTime = ob.getLastTested();
        return MathUtils.isNullValue(nTime);
    }

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const nTimeOne = obOne.getLastTested();
        const nTimeTwo = obTwo.getLastTested();

        const nResult =  nTimeOne < nTimeTwo ? 1 : nTimeOne > nTimeTwo ? -1 : 0;
        return nResult;
    }


}

CpdSemType.SimilaritySorter = class extends SemSorter
{
    getName() {return "Sort by Similarity";}
    supportsSortDirectionSwitch() {return false;}

    async fillSortValues(arValues, col, nRecordInitiator)
    {
        const typeSemCpd = GridUtils.getSemType(col);
        const service = typeSemCpd.getMolService();
        const strSmilesPattern = service.getSmilesOrMolFile(nRecordInitiator);

        await service.tanimotoScores(arValues, strSmilesPattern);
    }


    comparer(ctxOne, ctxTwo)
    {
        const fScoreOne = ctxOne.m_obValue;
        const fScoreTwo = ctxTwo.m_obValue;
        const nResult =  fScoreOne >= fScoreTwo ? -1 : fScoreOne < fScoreTwo ? 1 : 0;
        return nResult;
    }
}


CpdSemType.GNFIDSorter = class extends SemSorter
{
    getName() {return "Sort by GNF ID";}
    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        let str = obValue.getGnfRegId();
        return super.isNullSortValue(str);
    }

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const strOne = obOne.getGnfRegId();
        const strTwo = obTwo.getGnfRegId();
        const nResult =  strOne < strTwo ? -1 : strOne > strTwo ? 1 : 0;
        return nResult;
    }
}

CpdSemType.SIDSorter = class extends SemSorter
{
    getName() {return "Sort by Structure ID";}

    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        let str = obValue.getCpdSidId();
        return super.isNullSortValue(str);
    }

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const strOne = obOne.getCpdSidId();
        const strTwo = obTwo.getCpdSidId();
        const nResult =  strOne < strTwo ? -1 : strOne > strTwo ? 1 : 0;
        return nResult;
    }
}

CpdSemType.SampleIDSorter = class extends SemSorter
{
    getName() {return "Sort by Sample ID";}
    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        let str = obValue.getSampleId();
        return super.isNullSortValue(str);
    }

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const strOne = obOne.getSampleId();
        const strTwo = obTwo.getSampleId();
        const nResult =  strOne < strTwo ? -1 : strOne > strTwo ? 1 : 0;
        return nResult;
    }
}


CpdSemType.ConceptIDSorter = class extends SemSorter
{
    getName() {return "Sort by Concept ID";}
    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        let str = obValue.getConceptId();
        return super.isNullSortValue(str);
    }

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const strOne = obOne.getConceptId();
        const strTwo = obTwo.getConceptId();
        const nResult =  strOne < strTwo ? -1 : strOne > strTwo ? 1 : 0;
        return nResult;
    }
}
