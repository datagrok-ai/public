import {SemType} from "../../SemType";
import {SemSorter} from "../../SemSorter";

export class InVivoSemType extends SemType {

    constructor()
    {
        super(undefined, [new InVivoSemType.InVivoSorter()]);
    }
}

InVivoSemType.InVivoSorter = class extends SemSorter
{
    getName() {return "InVivo";}
    comparer(ctxOne, ctxTwo)
    {
        let obOne = ctxOne.m_obValue;
        let obTwo = ctxTwo.m_obValue;
        let bOne = obOne.hasValue();
        let bTwo = obTwo.hasValue();

        let nResult = bOne && !bTwo ? -1 : !bOne && bTwo ? 1 : 0;
        return nResult;
    }
}
