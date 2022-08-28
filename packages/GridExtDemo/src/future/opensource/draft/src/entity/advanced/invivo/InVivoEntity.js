import {SemEntity} from "../../SemEntity";

export class InVivoEntity extends SemEntity {
    constructor() {
        super();

        this.m_bHasValue = Math.random() >= 0.5;
    }

    getDefaultValue()
    {
        return this.hasValue();
    }

    hasValue() {
        return this.m_bHasValue;
    }
}
