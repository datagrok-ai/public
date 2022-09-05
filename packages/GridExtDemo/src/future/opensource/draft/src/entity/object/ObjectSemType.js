import {SemType} from "../SemType";
import {SemSorter} from "../SemSorter";
import {ObjectEntity} from "./ObjectEntity";
import {MathUtils} from "../../utils/MathUtils";

export class ObjectSemType extends SemType {

    constructor(strIconFilePath) {
        super(strIconFilePath, [new ObjectSemType.ObjectSorter()]);
    }

    getEntityClass()
    {
        return ObjectEntity;
    }

    getCompatibilityLevel(arCols)
    {
        return arCols.length > 0 ? 1 : 0;
   }

    createColName(arCols)
    {
        const col = arCols[0]; //0 must exist
        return col.name;
    }


    createEntity(arCols, nRow)
    {
        const col = arCols[0]; //0 must exist
        let ob = col.get(nRow);

        if(MathUtils.isNullValue(ob))
        {
           ob = null;
        }

        return new ObjectEntity(ob);
    }
}

ObjectSemType.ObjectSorter = class extends SemSorter
{
    getName() {return "Value";}
    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const valOne = obOne.getValue();
        const valTwo = obTwo.getValue();

        const nResult =  valOne < valTwo ? -1 : valOne > valTwo ? 1 : 0;
        return nResult;
    }
}





