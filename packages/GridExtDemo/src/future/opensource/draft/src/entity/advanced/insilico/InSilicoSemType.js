import {SemType} from "../../SemType";
import {SemSorter} from "../../SemSorter";

export class InSilicoSemType extends SemType {

    constructor()
    {
        super(undefined,  [new InSilicoSemType.LipinskiSorter(), new InSilicoSemType.MolWeightSorter(), new InSilicoSemType.CLogPSorter()]);
    }
}


InSilicoSemType.LipinskiSorter = class extends SemSorter
{
    getName() {return "Lipinski Violation";}
    comparer(ctxOne, ctxTwo)
    {
        let obOne = ctxOne.m_obValue;
        let obTwo = ctxTwo.m_obValue;
        let nLipOne = obOne.getLipinskiViolations();
        let nLipTwo = obTwo.getLipinskiViolations();

        let nResult =  nLipOne < nLipTwo ? -1 : nLipOne > nLipTwo ? 1 : 0;
        return nResult;
    }
}

InSilicoSemType.MolWeightSorter = class extends SemSorter
{
    getName() {return "MW";}
    comparer(ctxOne, ctxTwo)
    {
        let obOne = ctxOne.m_obValue;
        let obTwo = ctxTwo.m_obValue;
        let fOne = obOne.getMW();
        let fTwo = obTwo.getMW();

        let nResult =  fOne < fTwo ? -1 : fOne > fTwo ? 1 : 0;
        return nResult;
    }
}

InSilicoSemType.CLogPSorter = class extends SemSorter
{
    getName() {return "CLogP";}
    comparer(ctxOne, ctxTwo)
    {
        let obOne = ctxOne.m_obValue;
        let obTwo = ctxTwo.m_obValue;
        let fOne = obOne.getCLOGP();
        let fTwo = obTwo.getCLOGP();

        let nResult =  fOne < fTwo ? -1 : fOne > fTwo ? 1 : 0;
        return nResult;
    }
}