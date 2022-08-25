import {SemEntity} from "../../SemEntity";
import {FZNumber} from "../../../lang/FZNumber";

export class IC50Entity extends SemEntity {
    constructor() {
        super();

        this.m_fMedianIC50 = Math.floor(Math.random() * 20);

        const nM = Math.floor(Math.random() * FZNumber.MODIFIERS.length);
        this.m_modifierFuzzy = FZNumber.MODIFIERS[nM];//FZNumber.parseModifier(value.toString());
        //this.m_nFuzzySortRating = modifier.getSortRating();

        this.m_fStdv  = Math.floor(Math.random() * 10)/10;
        this.m_fPctEff = Math.floor(Math.random() * 100) + (Math.floor(Math.random()*10)/10);

        const nMonth = 6;//Math.floor(Math.random() * 7);
        const nDay =  /*20 + */Math.floor(Math.random() * 29);
        this.m_nLastTested = new Date(2021, nMonth, nDay).getTime();
    }

    getDefaultValue()
    {
        return this.getMedianIC50();
    }

    getMedianIC50()
    {
        return this.m_fMedianIC50;
    }

    getFuzzyModifier()
    {
        return this.m_modifierFuzzy;
    }

    getPctEff()
    {
        return this.m_fPctEff;
    }

    getStdv()
    {
        return this.m_fStdv;

    }

    getLastTested() {return this.m_nLastTested;}
}
