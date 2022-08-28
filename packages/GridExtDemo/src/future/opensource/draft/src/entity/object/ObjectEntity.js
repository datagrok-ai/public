import {SemEntity} from "../SemEntity";

export class ObjectEntity extends SemEntity {
    constructor(ob)
    {
        super();

        this.m_ob = ob;
    }

    getDefaultValue() {
        return this.getValue();
    }

    getValue()
    {
        return this.m_ob;
    }
}