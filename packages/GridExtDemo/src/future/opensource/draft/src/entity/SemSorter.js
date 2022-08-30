import {MathUtils} from "../utils/MathUtils";
import {NotImplementedError} from "../lang/NotImplementedError";

/**
 * The SemSorter class is an abstract class that defines an interface used by semantic types to perform Sort action.
 * Each instance of the SemSorter class represent one sorting option for the corresponding SemType instance. In the object-oriented paradigm,
 * one semantic type can have multiple sorting options. Each option can use a single or multiple properties from SemEntity objects or even the data
 * from other columns of the DataFrame.
 */
export class SemSorter {

    /**
     * Returns the name of the sorting option.
     */
    getName() {throw new NotImplementedError();}

    /**
     * Returns the comparator function to be used in the sort action.
     * @param ctxOne the first object to be compared
     * @param ctxTwo the secont object to be compared.
     * @returns {Number} a number whose sign is used to determine theresult of the compare opertation.
     */
    comparer(ctxOne, ctxTwo) {throw new NotImplementedError();}

    /**
     * Deteermines whether this sorting option supports switching between ascending and descending sort directions.
     * @returns {boolean} true if switching between directions is supportd; otherwise false.
     */
    supportsSortDirectionSwitch() {return true;}

    /**
     * Determines whether the specified value must be treated as non-existing (or NULL) value.
     * @param obValue the value to be tested.
     * @returns {boolean} true if the value must be treated as NULL, otherwise false.
     */
    isNullSortValue(obValue)
    {
        return MathUtils.isNullValue(obValue);
    }


    /**
     * Called by the corresponding SemType to Fill the specified array with the data to be sorted. THe default implementation
     * fills with the data from the correspnding object column of the DataFrame.
     * @param arValues an array to be filled with the data to be sorted. THe data can be either JS primitive values
     * like Number, String, Boolean, or can be complex objects. The filled values will be passed to the comparer function.
     * @param col the DataFrame's column on behalf of which the Sort operation is being performed (user clicked on the column header).
     * @param nRowInitiator the index of the row in the DataFrame on behalf of which the Sort operation has been initiated (sort initiateb from the context menu
     * when use right clicked on a row). -1 if no row is involved.
     * @returns {Promise<void>} a promise object,.
     */
    async fillSortValues(arValues, col, nRowInitiator)
    {
        let obValue = null;
        let nROwCount = col.length;
        for(var nR=0; nR<nROwCount; ++nR)
        {
            obValue = col.get(nR);
            arValues[nR] = obValue;
        }
    }
}