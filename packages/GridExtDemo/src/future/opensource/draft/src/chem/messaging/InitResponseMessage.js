export class InitResponseMessage
{
    constructor(nMolCount, bCancelled)
    {
        this.m_nMolCount = nMolCount;
        this.m_bCancelled = bCancelled === undefined ? false : bCancelled;
    }

    getMolCount() {return this.m_nMolCount;}
    isCancelled() {
        return this.m_bCancelled;
    }
}
