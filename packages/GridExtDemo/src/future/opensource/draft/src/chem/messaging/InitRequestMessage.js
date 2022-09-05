export class InitRequestMessage
{
    constructor(strWebRoot, arSmilesOrMolFiles, nFr, nTo, bMolFIle, nCPU)
    {
        this.m_strWebRoot = strWebRoot;
        this.m_arSmilesOrMolFiles = arSmilesOrMolFiles;
        this.m_nFr = nFr;
        this.m_nTo = nTo;
        this.m_bMolFIle = bMolFIle;
        this.m_nCPU = nCPU;
    }

    getCPUId() {
        return this.m_nCPU;
    }

    getWebRoot()
    {
        return this.m_strWebRoot;
    }

    getSmilesOrMolFiles()
    {
        return this.m_arSmilesOrMolFiles;
    }

    isMolFile()
    {
        return this.m_bMolFIle;
    }

    getIdxFr() {return this.m_nFr;}
    getIdxTo() {return this.m_nTo;}

}
