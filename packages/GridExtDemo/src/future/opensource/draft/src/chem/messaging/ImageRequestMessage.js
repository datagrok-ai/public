export class ImageRequestMessage
{
    constructor(nIdxMol, nW, nH)
    {
        this.m_nIdxMol = nIdxMol;
        this.m_nW = nW;
        this.m_nH = nH;
    }

    getMolIndex()
    {
        return this.m_nIdxMol;
    }

    getW()
    {
        return this.m_nW;
    }

    getH()
    {
        return this.m_nH;
    }
}