import {DateEntity} from "../date/DateEntity";

export class TextDateEntity extends DateEntity {
    constructor(strText, nTime)
    {
        super(nTime);

        this.m_strText = strText;
    }

    getDefaultValue() {
        return this.getText();
    }

    getText()
    {
        return this.m_strText;
    }
}