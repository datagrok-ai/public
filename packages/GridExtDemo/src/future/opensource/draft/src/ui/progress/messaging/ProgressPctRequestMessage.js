export class ProgressPctRequestMessage
{
    constructor(nPct)
    {
        this.m_nPct = nPct;
    }

    getPct()
    {
        return this.m_nPct;
    }
}