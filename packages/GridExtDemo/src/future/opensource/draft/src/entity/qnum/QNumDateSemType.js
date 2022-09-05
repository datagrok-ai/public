import {SemType} from "../SemType";
import {SemSorter} from "../SemSorter";
import {DateSemType} from "../date/DateSemType";
import {CDFGenericQualifiedNumberType} from "../../service/nx/cdf/composite/CDFGenericQualifiedNumberType";
import {CDFGenericDateType} from "../../service/nx/cdf/composite/CDFGenericDateType";
import {FZModifier} from "../../lang/FZModifier";
import {QNumDateEntity} from "./QNumDateEntity";
import {ErrorUtils} from "../../utils/ErrorUtils";
import {DateUtils} from "../../utils/DateUtils";

export class QNumDateSemType extends DateSemType {

    constructor(strIconFilePath) {
        super(strIconFilePath, [new QNumDateSemType.ValueSorter(), new DateSemType.DateSorter()]);
    }


    getIconFilePath() {return "images/QTime.png";}

    getEntityClass()
    {
        return QNumDateEntity;
    }


    fillHeaderValues(arHeaderValues, arCols)
    {
        SemType.fillHeaderValues(arHeaderValues, arCols, QNumDateSemType.CDF_TYPES_HEADER_VALUES);
    }

    getCompatibilityLevel(arCols)
    {

        const arTypeKeys = QNumDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arCols);
        const arHitColIdxs = new Array(arTypeKeys.length);
        const nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);
        if(nHitCount === 0)
        {
            arHitColIdxs.fill(-1);
            let col = 0;
            for(var nC=0; nC<arCols.length; ++nC)
            {
                col = arCols[nC];
                if(col.type === DG.COLUMN_TYPE.QNUM)
                {
                    arHitColIdxs[0] = nC;
                    return 0.1;
                }
            }
        }

        return arHitColIdxs[0]  >= 0 ? nHitCount : 0;
   }

    createColName(arCols)
    {
        const arTypeKeys = QNumDateSemType.CDF_TYPES;
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
        const arTypeKeys = QNumDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arCols);
        const arHitColIdxs = new Array(arTypeKeys.length);
        const nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);

        if(nHitCount === 0)
        {
            arHitColIdxs.fill(-1);
            let col = 0;
            for(var nC=0; nC<arCols.length; ++nC)
            {
                col = arCols[nC];
                if(col.type === DG.COLUMN_TYPE.QNUM)
                {
                    arHitColIdxs[0] = nC;
                    break;
                    //return 1;
                }
            }
        }

        if(arHitColIdxs[0] < 0) {

            SemType.fillCompatibilityLevelHits(arHitColIdxs, arKeys, arCols);
            throw new Error("Qualified number was not identified.");
        }
        let col = arCols[arHitColIdxs[0]]; //0 must exist
        const qnum = col.get(nRow);
        let fValue = null;
        let sign = null;
        if(qnum !== null && qnum !== undefined && qnum !== DG.FLOAT_NULL)
        {
            fValue = DG.Qnum.getValue(qnum);
            sign = DG.Qnum.getQ(qnum);
            sign = sign === DG.QNUM_LESS ? FZModifier.LESS_THAN : sign === DG.QNUM_GREATER ? FZModifier.GREATER_THAN : FZModifier.EXACT_EQUAL_TO;
        }

        let nTime = NaN;
        if(arHitColIdxs[1] >= 0)
        {
            col = arCols[arHitColIdxs[1]];
            const type = col.type;
            const date = col.get(nRow);
            nTime = DateUtils.getTimeFromDGDate(date);

           /*
            if(date === "")
            {
                nTime = NaN;
            }


            else if(date !== null && date !== undefined)
              nTime = date.a;

            if(typeof nTime !== Number.name.toLowerCase())
            {
                nTime = NaN;
            }*/
        }

        return new QNumDateEntity(fValue, sign, nTime);
    }
}


QNumDateSemType.CDF_TYPES = [CDFGenericQualifiedNumberType.nValue, CDFGenericDateType.dValue];
QNumDateSemType.CDF_TYPES_HEADER_VALUES = new Map([[QNumDateSemType.CDF_TYPES[0], "Value"],
                                                           [QNumDateSemType.CDF_TYPES[1],  "Date"]]);


QNumDateSemType.ValueSorter = class extends SemSorter
{
    getName() {return "Value";}

    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        const fValue = obValue.getValue();
        return super.isNullSortValue(fValue);
    }

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





