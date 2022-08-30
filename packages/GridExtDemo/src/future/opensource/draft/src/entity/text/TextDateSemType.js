import {SemType} from "../SemType";
import {SemSorter} from "../SemSorter";
import {DateSemType} from "../date/DateSemType";
import {CDFGenericDateType} from "../../service/nx/cdf/composite/CDFGenericDateType";
import {CDFGenericTextType} from "../../service/nx/cdf/composite/CDFGenericTextType";
import {TextDateEntity} from "./TextDateEntity";
import {DateUtils} from "../../utils/DateUtils";

export class TextDateSemType extends DateSemType {

    constructor(strIconFilePath) {
        super(strIconFilePath, [new TextDateSemType.TextSorter(), new DateSemType.DateSorter()]);
    }

    getIconFilePath() {return "images/TTime.png";}

    getEntityClass()
    {
        return TextDateEntity;
    }

    fillHeaderValues(arHeaderValues, arCols)
    {
        SemType.fillHeaderValues(arHeaderValues, arCols, TextDateSemType.CDF_TYPES_HEADER_VALUES);
    }

    getCompatibilityLevel(arCols)
    {
        const arTypeKeys = TextDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arCols);
        const arHitColIdxs = new Array(arTypeKeys.length);
        let nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);

        return arHitColIdxs[0]  >= 0 ? nHitCount : 0;
   }

    createColName(arCols)
    {
        const arTypeKeys = TextDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arCols);
        const arHitColIdxs = new Array(arTypeKeys.length);
        const nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);
        if(arHitColIdxs[0] < 0)
            throw new Error("Qualified number was not identified.");

        const col = arCols[arHitColIdxs[0]]; //0 must exist
        return col.name;
    }

    createEntity(arCols, nRow)
    {
        const arTypeKeys = TextDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arCols);
        const arHitColIdxs = new Array(arTypeKeys.length);
        let nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);
        if(arHitColIdxs[0] < 0)
            throw new Error("Text value was not identified.");

        let col = arCols[arHitColIdxs[0]]; //0 must exist
        let str = col.get(nRow);

        if(str === undefined)
        {
           str = null;
        }

        let nTime = NaN;
        if(arHitColIdxs[1] >= 0)
        {
            col = arCols[arHitColIdxs[1]];
            const date = col.get(nRow);
            nTime = DateUtils.getTimeFromDGDate(date);
        }

        return new TextDateEntity(str, nTime);
    }
}


TextDateSemType.CDF_TYPES = [CDFGenericTextType.tValue, CDFGenericDateType.dValue];
TextDateSemType.CDF_TYPES_HEADER_VALUES = new Map([[TextDateSemType.CDF_TYPES[0], "Text"],
                                                           [TextDateSemType.CDF_TYPES[1],  "Date"]]);


TextDateSemType.TextSorter = class extends SemSorter
{
    getName() {return "Value";}
    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        const str = obValue.getText();
        return super.isNullSortValue(str);
    }

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const valOne = obOne.getText();
        const valTwo = obTwo.getText();

        const nResult =  valOne < valTwo ? -1 : valOne > valTwo ? 1 : 0;
        return nResult;
    }
}





