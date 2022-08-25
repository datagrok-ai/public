export class FZModifier
{
    //Construction
    constructor(strName, nSortRating)
    {
        this.m_strName = strName;
        this.m_nSortRating = nSortRating;
    }

//IMPLEMENTATION SECTION
    getSortRating() {return this.m_nSortRating;}

    toString() {return this.m_strName;}

//DATA SECTION
}
FZModifier.SORT_RATING = 0;
FZModifier.LESS_THAN = new FZModifier("<", FZModifier.SORT_RATING++);
FZModifier.LESS_THAN_OR_EQUAL = new FZModifier("<=", FZModifier.SORT_RATING++);
FZModifier.APPROX_EQUAL_TO = new FZModifier("~", FZModifier.SORT_RATING++);
FZModifier.EXACT_EQUAL_TO = new FZModifier("=", FZModifier.SORT_RATING++);
FZModifier.GREATER_THAN_OR_EQUAL = new FZModifier(">=", FZModifier.SORT_RATING++);
FZModifier.GREATER_THAN = new FZModifier(">", FZModifier.SORT_RATING++);