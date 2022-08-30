import {SemEntity} from "../SemEntity";

export class DateEntity extends SemEntity
{
    constructor(nTime)
    {
        super();

        this.m_nTime = nTime;
    }

    getDefaultValue() {
        return getTime();
    }

    getTime() {return this.m_nTime;}
}
