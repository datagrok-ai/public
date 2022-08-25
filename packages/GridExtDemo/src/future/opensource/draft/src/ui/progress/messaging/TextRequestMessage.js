export class TextRequestMessage
{
    constructor(strText, bUp)
    {
        this.m_strText = strText;
        this.m_bUp = bUp;
    }

    getText()
    {
        return this.m_strText;
    }

    isUp()
    {
        return this.m_bUp;
    }

}