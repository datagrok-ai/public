import * as DG from "datagrok-api/dg";
import {AnalysisLoader} from "../service/nx/analysis/AnalysisLoader";
import {NotImplementedError} from "../lang/NotImplementedError";
import {ErrorUtils} from "../utils/ErrorUtils";

/**
 * SemType class defines a template for semantic types assigned to object columns in DataFrames.
 * Semantic type helps the application and its vusualizations interprete the column data in a meaningful way
 * either by rendering it, or be performing operations like Sort that make sense for the particular type of data.
 * Any semantic type is tightly coupled with the corresponding SemEntity object and can be considered as a "manager"
 * of the data represented by SemEntity.
 */
export class SemType {
    /**
     * Constructs a new instance of SemType class and initializes it with a displayable icon and sort options.
     * @param strIconFilePath a URL to the icon image that visually describes this semantic type.
     * @param arSorters an array of sorting options represeted by instances of the SemType class.
     */
    constructor(strIconFilePath, arSorters) {

        this.m_arSorters = arSorters;
        this.m_strUrlIcon = null;
        this.m_strIconFilePath = null;
        if(strIconFilePath !== undefined)
            this.m_strIconFilePath = strIconFilePath;
    }
    //api
    /**
     * Called by the application framework to perform initialization of this type instance.
     * @param obPackage a reference to the DataGrok package.
     * @param arCols an array containing child columns for this typee.
     */
    async init(obPackage, arCols) {
        const strPath = this.getIconFilePath();
        if(strPath !== null)
        {
            const strURLWR = obPackage.webRoot;
            this.m_strUrlIcon = strURLWR + strPath;
        }

        if(arCols !== undefined && arCols !== null)
        {
            const nColCount = arCols.length;
            const arHeaderVals = new Array(nColCount);
            this.fillHeaderValues(arHeaderVals, arCols);

            for(var nC=0; nC<nColCount; ++nC)
            {
                arCols[nC].setTag(SemType.KEY_TAG_CHILD_COLUMN_HEADEER_VALUE, arHeaderVals[nC]);
            }
        }
    }


    getIconURL() {return this.m_strUrlIcon;}
    getIconFilePath() {return this.m_strIconFilePath;}

    getEntityClass()
    {
        throw new NotImplementedError();
    }

    fillHeaderValues(arHeaderVals, arCols)
    {
        const nColCount = arCols.length;
        for(var nC=0; nC<nColCount; ++nC)
        {
         arHeaderVals[nC] = arCols[nC].name;
        }
    }

    getCompatibilityLevel(arCols)
    {
        throw new NotImplementedError();
    }

    createEntity(arCols, nRow)
    {
        throw new NotImplementedError();
    }

    createColName(arCols)
    {
        return null;
    }

    supportsSortDirectionSwitch(nOption)
    {
        return  this.m_arSorters[nOption].supportsSortDirectionSwitch();
    }

    /**
     * Returns the number sorting options supported by this type.
     * @returns {number|*}  the number sorting options.
     */
    getSortOptionCount()
    {
        return this.m_arSorters === undefined ? 0 : this.m_arSorters.length;
    }

    /**
     * Returns the name of the sorting option at a specified index.
     * @param nOption the specified option's index.
     * @returns {String} a string containing the name of the sorting option.
     */
    getSortOptionName(nOption)
    {
        return this.m_arSorters[nOption].getName();
    }

    /**
     * Determines whether the specified value for a sorting option must be treated as non-existing (or NULL) value.
     * @param obValue the value to be tested.
     * @param nOption the sorting option's index.
     * @returns {boolean|*} true if the value must be treated as NULL, otherwise false.
     */
    isNullSortOptionValue(obValue, nOption) {
        return this.m_arSorters[nOption].isNullSortValue(obValue);
    }


    /**
     * Returns the comparator function to be used for the specified sort option.
     * @param nOption the sorting option's index.
     * @returns the comparator function.
     */
    getSortOptionComparer(nOption)
    {
        return this.m_arSorters[nOption].comparer;
    }

    /**
     * Returns the insex of the default sort option.
     * @returns {number} by default return 0, or -1 of the number of sorting options is 0.
     */
    getDefaultSortOptionIndex()
    {
        let nCount = this.getSortOptionCount();
        return nCount === 0 ? -1 : 0;
    }



    /**
     * Called by the application framework to Fill the specified array with the data to be sorted using a specified sort option.
     * The default implementation fills with the data from the correspnding object column of the DataFrame.
     * @param arValues an array to be filled with the data to be sorted. THe data can be either JS primitive values
     * @param col the DataFrame's column on behalf of which the Sort operation is being performed (user clicked on the column header).
     * @param the index of the row in the DataFrame on behalf of which the Sort operation has been initiated (sort initiateb from the context menu
     * when use right clicked on a row). -1 if no row is involved.
     * @param nOption the sorting option's index.
     * @returns {Promise<void>} a JS promise
     */
    async fillSortOptionValues(arValues, col, nRowInitiator, nOption)
    {
        await this.m_arSorters[nOption].fillSortValues(arValues, col, nRowInitiator);
    }
}



SemType.KEY_TAG_CHILD_COLUMN_HEADEER_VALUE = "KEY_TAG_CHILD_COLUMN_HEADEER_VALUE";


SemType.SortCtx = class {
    constructor(nIdx, obValue)
    {
        this.m_nIdx = nIdx;
        this.m_obValue = obValue;
    }
}

SemType.getChildColHeaderValue = function(col)
{
    const valueHeader = col.getTag(SemType.KEY_TAG_CHILD_COLUMN_HEADEER_VALUE);
    return valueHeader === undefined ? null : valueHeader;
}


/**
 * Returns an array filled with the name of sorting options for a instance of semantic type.
 * @param typeSem the specified instance of semanic type.
 * @returns {String[]} a reference to and array object containing sorting options names.
 */
SemType.toSortOptionNamesArray = function(typeSem)
{
    let str = "";
    let nCount = typeSem.getSortOptionCount();
    let ar = new Array(nCount);
    for(var n=0; n<nCount; ++n)
    {
        str = typeSem.getSortOptionName(n);
        ar[n] = str;
    }

    return ar;
}


SemType.sortOptionIndexFromName = function(typeSem, strName)
{
    if(typeSem === undefined)
        throw new Error("Symantec type cannot be undefined.");

    if(strName === undefined)
        throw new Error("Sort option name cannot be undefined.");

    let str = "";
    let nCount = typeSem.getSortOptionCount();
    for(var n=0; n<nCount; ++n)
    {
        str = typeSem.getSortOptionName(n);
        if(strName === str)
            return n;
    }

    return -1;
}

/**
 * Called by the application framework to perform sorting operation for a particular sort option of a semantic type.
 * This function is intended for internal use only and is not a part of the public API.
 * @param arIndices an array of original indices corresponding to the records before sorting
 * @param col the colum on behalf of which the sorting is to be performed.
 * @param nRowInitiator the index of the row in the DataFrame on behalf of which the Sort operation has been initiated (sort initiateb from the context menu
 * when use right clicked on a row). -1 if no row is involved.
 * @param nOption the sorting option's index.
 * @param bAscend {boolean} indicates whether the sort direction is acsenting or descending.
 * @returns {Promise<void>} a JS promise.
 */
SemType.sort = async function(arIndices, col, nRowInitiator, nOption, bAscend) {
    //ErrorUtils.verifyType(arIndices, Array);
    ErrorUtils.verifyClass(col, DG.Column);

    const nRecordCount = col.length;
    if(nRecordCount !== arIndices.length)
        throw new Error("The return array's size (" + arIndices.length + ") is different from the column size (" + nRecordCount + ")");

    const arValues = new Array(nRecordCount);
    const typeSem = col.dart.m_typeSem;
    await typeSem.fillSortOptionValues(arValues, col, nRowInitiator, nOption);

    const arCtxs = [];
    const arNullIdxs = [];
    let nNullCount = 0;
    let obValue = null;
    let ctx = null;
    for(let nR=0; nR<nRecordCount; ++nR) {
        obValue = arValues[nR];//my changes col.get(nR);
        if(obValue !== null && !typeSem.isNullSortOptionValue(obValue, nOption))
        {
            ctx = new SemType.SortCtx(nR, obValue);
            arCtxs.push(ctx);
        }
        else
        {
            arNullIdxs.push(nR);
            ++nNullCount;
        }
    }

    const fnCompare = typeSem.getSortOptionComparer(nOption);
    if(fnCompare === null)
        arCtxs.sort();
    else
        arCtxs.sort(fnCompare);

    const nNonNullCount = nRecordCount - nNullCount;
    let n = 0;
    if(!bAscend) {
        let ctxTmp = null;

        let nNonNullCountHalf = Math.floor(nNonNullCount / 2);
        for (n = 0; n < nNonNullCountHalf; ++n) {
            ctxTmp = arCtxs[n];
            arCtxs[n] = arCtxs[nNonNullCount - n - 1];
            arCtxs[nNonNullCount - n - 1] = ctxTmp;
        }
    }

    for(n=0; n<nNonNullCount; ++n)
    {
        arIndices[n] = arCtxs[n].m_nIdx;
    }

    for(n=nNonNullCount; n<nNonNullCount + nNullCount; ++n)
    {
        arIndices[n] = arNullIdxs[n-nNonNullCount];
    }
}



SemType.col2CDFPrimTypes= function(arCols)
{
    const nColCount = arCols.length;
    let arCDFTypes = new Array(nColCount);
    let obTagDescriptor = null;
    let col = null;
    let typeCDFPrimitive = null;
    for(var nC=0; nC<nColCount; ++nC) {
        col = arCols[nC];

        obTagDescriptor = col.getTag(AnalysisLoader.TAG_COMPPOSE_DESCRIPTOR);
        if (obTagDescriptor === null || typeCDFPrimitive === undefined) {
            arCDFTypes[nC] = null;
            continue;
        }
        typeCDFPrimitive = obTagDescriptor.getPrimitiveType();
        arCDFTypes[nC] = typeCDFPrimitive;
    }

    return arCDFTypes;
}


SemType.fillCompatibilityLevelHits = function(arHitColIdxs, arTypeKeys, arColKeys)
{
    if(arTypeKeys.length !== arHitColIdxs.length)
        throw new Error("The length of the array containing keys must be the same as that of array to contain column indices, but they are different: " + arTypeKeys.length + " " + arColKeys.length)

    arHitColIdxs.fill(-1);

    let nHitCount = 0;
    let obColKey = null;
    let typeCol = null;

    for(let nC=0; nC<arColKeys.length; ++nC)
    {
        obColKey = arColKeys[nC];
        typeCol = obColKey === null ? null : obColKey.m_type;

        if(obColKey === null || obColKey === undefined)
            continue;

        for(var nKey=0; nKey<arTypeKeys.length; ++nKey)
        {
            if(arHitColIdxs[nKey] === -1 && obColKey === arTypeKeys[nKey]) {
                arHitColIdxs[nKey] = nC;
                ++nHitCount;

                if(nHitCount === arTypeKeys.length)
                    return nHitCount;

                break;
            }
        }

        if(nHitCount === arTypeKeys.length)
            return nHitCount;
    }

    return nHitCount;
}


SemType.fillHeaderValues = function(arHeaderValues, arCols, mapKeys2Titles)
{
    ErrorUtils.verifyNotUndefinedNotNull(arHeaderValues);
    ErrorUtils.verifyNotUndefinedNotNull(arCols);

    const arColKeys = SemType.col2CDFPrimTypes(arCols);
    let key = null;
    let str = null;
    const nColCount = arCols.length;
    for(let nC=0; nC<nColCount; ++nC)
    {
        key = arColKeys[nC];
        if(key === null || key === undefined)
        {
            arHeaderValues[nC] = arCols[nC].name;
            continue;
        }

        str = mapKeys2Titles.get(key);
        if(str === undefined || str === null)
        {
            arHeaderValues[nC] = arCols[nC].name;
            continue;
        }

        arHeaderValues[nC] =  "[" + str + "] " + arCols[nC].name;
    }
}


