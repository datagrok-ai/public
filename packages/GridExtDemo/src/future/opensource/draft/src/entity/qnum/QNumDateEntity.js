import {DateEntity} from "../date/DateEntity";

export class QNumDateEntity extends DateEntity {
    constructor(fValue, modifierFuzzy, nTime)
    {
        super(nTime);

        this.m_fValue = fValue;
        this.m_modifierFuzzy = modifierFuzzy;
    }

    getDefaultValue() {
        return this.getValue();
    }

    getValue()
    {
        return this.m_fValue;
    }

    getFuzzyModifier()
    {
        return this.m_modifierFuzzy;
    }
}