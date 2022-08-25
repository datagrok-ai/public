import {SemType} from "../../SemType";
import {SemSorter} from "../../SemSorter";

export class IC50SemType extends SemType {

    constructor()
    {
        super(undefined, [new IC50SemType.MedianIC50Sorter(), new IC50SemType.LastTestedSorter()]);
    }
}


IC50SemType.MedianIC50Sorter = class extends SemSorter
{
    getName() {return "Median IC50";}
    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const valOne = obOne.getMedianIC50();
        const valTwo = obTwo.getMedianIC50();

        const nResult =  valOne < valTwo ? -1 : valOne > valTwo ? 1 : 0;
        return nResult;
    }
}


IC50SemType.LastTestedSorter = class extends SemSorter
{
    getName() {return "Last Tested";}
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

IC50SemType.Instance = new IC50SemType();
