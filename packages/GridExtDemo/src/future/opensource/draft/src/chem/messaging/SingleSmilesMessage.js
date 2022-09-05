export class SingleSmilesMessage
{
    constructor(strSmilesFrag)
    {
        this.m_strSmilesFrag = strSmilesFrag;

    }

    getFragmentSmiles()
    {
        return this.m_strSmilesFrag;
    }

}