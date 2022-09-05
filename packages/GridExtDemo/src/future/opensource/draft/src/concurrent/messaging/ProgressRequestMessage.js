export class ProgressRequestMessage
{
    constructor(nJobCount)
    {
        this.m_nJobCount = nJobCount;
    }

    getJobCount() {return this.m_nJobCount;}
}