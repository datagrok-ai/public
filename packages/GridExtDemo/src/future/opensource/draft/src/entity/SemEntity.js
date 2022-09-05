import {NotImplementedError} from "../lang/NotImplementedError";

/**
 *  The SemEntity class defines a top-level interface for object-oriented data entities used
 *  by virtual columns of DataFrames. SemEntity object a tightly coupled with the corresponding DataFrame column type represented by
 *  SemType class. When a new column type introduced, both SemEntity and SemType should be implemented in a way that SemType could interpret
 *  the content of the corresponding SemEntity.
 */
export class SemEntity {

    /**
     * Returns the default value.
     */
    getDefaultValue()
    {
        throw new NotImplementedError();
    }
}