import {FZModifier} from "./FZModifier";

export class FZNumber
{
//Construction
    constructor(modifier)
    {
        //if(!isValidModifier(nModifier))
        if(modifier === null)
            throw new Error("Modifier cannot be null.");

        this.m_modifier = FZModifier.EXACT_EQUAL_TO;
        if(modifier !== undefined)
            this.m_modifier = modifier;
    }



//IMPLEMENTATION SECTION

    /**
     * Returns the value's modifier.
     * @return one of the following integer constants <code>APPROX_EQUAL_TO</code>,
     * <code>EXACT_EQUAL_TO</code>, <code>GREATER_THAN</code>,
     * <code>GREATER_THAN_OR_EQUAL_TO</code>, <code>LESS_THAN</code>,
     * <code>LESS_THAN_OR_EQUAL_TO</code>.
     */
    getModifier() {return this.m_modifier;}

    compareTo(o)
    {
        if(o === null)
            return 1;

        let nr = o;
        let nRTh = this.m_modifier.m_nSortRating; //SORTING_RATING[m_nModifier];
        let nRNr = nr.m_modifier.m_nSortRating;// SORTING_RATING[nr.m_nModifier];

        return nRTh > nRNr ? 1 : nRTh < nRNr ? -1 : 0;
    }

    toString()// {return m_nModifier == 0 ? "" : MODIFIERS[m_nModifier];}
    {
        return this.m_modifier.toString();
    }

//DATA SCTION


//public static final int EXACT_EQUAL_TO =0;
//public static final int APPROX_EQUAL_TO = 1;
//public static final int GREATER_THAN_OR_EQUAL_TO =2;
//public static final int GREATER_THAN =3;
//public static final int LESS_THAN_OR_EQUAL_TO = 4;
//public static final int LESS_THAN = 5;

//Modifier.EXACT_EQUAL_TO, Modifier.APPROX_EQUAL_TO,
//Modifier.GREATER_THAN_OR_EQUAL, Modifier.GREATER_THAN,
//Modifier.LESS_THAN_OR_EQUAL, Modifier.LESS_THAN};

//private static final String MODIFIERS [] = {"=", "~", ">=", ">", "<=", "<"};
//private static final int SORTING_RATING [] = {3, 2, 4, 5, 1, 0};//">", ">=", "=", "~",  "<=", "<"
//private int m_nModifier = 0;
}

FZNumber.EXCEPTION_INVALID_MODIFIER_VALUE  = "Invalid modifier value.";
FZNumber.EXCEPTION_INVALID_MODIFIER_STRING = "Invalid modifier string.";

FZNumber.MODIFIERS = [FZModifier.LESS_THAN, FZModifier.LESS_THAN_OR_EQUAL,
    FZModifier.APPROX_EQUAL_TO, FZModifier.EXACT_EQUAL_TO,
    FZModifier.GREATER_THAN_OR_EQUAL, FZModifier.GREATER_THAN];


FZNumber.modifierFromRating = function(nRating)
{
    if(nRating < 0 || nRating >= FZNumber.MODIFIERS.length)
        throw new Error("Rating is out of range: " + nRating);

    return FZNumber.MODIFIERS[nRating];
}


FZNumber.parseModifier = function(strFuzzy)
{
    if(strFuzzy === indefined)
        throw new Error("String containing fuzzy sign cannot be undefined.");

    if(strFuzzy === "")
        return FZModifier.EXACT_EQUAL_TO;

    let nModifierCount = FZNumber.MODIFIERS.length;
    for(let n=0; n<nModifierCount; ++n)
    {
        if(strFuzzy === FZNumber.MODIFIERS[n].m_strName)
            return FZNumber.MODIFIERS[n];
    }

    throw new Error(FZNumber.EXCEPTION_INVALID_MODIFIER_STRING + ": " + strFuzzy);
}

FZNumber.parseModifierFromFZNumber = function(strNumber)
{
    if(strNumber === null)
        return null;

    let strModifier = null;
    let nCharCount =0;
    let nModifierCount = MODIFIERS.length;
    for(var n=0; n<nModifierCount; ++n)
    {
        nCharCount = FZNumber.MODIFIERS[n].m_strName.length;
        try{strModifier = strNumber.substring(0, nCharCount);}
        catch(ex) {continue;}

        if(strModifier === "")
            return FZModifier.EXACT_EQUAL_TO;

        if(strModifier === FZNumber.MODIFIERS[n].m_strName)
            return FZNumber.MODIFIERS[n];
    }

    throw new Error(FZNumber.EXCEPTION_INVALID_MODIFIER_STRING + ": " + strModifier);
}

