import {createRDKIt} from "./RDKit_minimal_2021.03_02.js";
import {MolWorker} from "../MolWorker";

export class RDKitWorker extends MolWorker
{
    constructor(fnDispose)
    {
     super(fnDispose);

        this.m_arMols = null;
        this.m_viewFs = null;
        this.m_rdKitModule = null;
    }

    dispose()
    {
        const nMolCount = this.m_arMols.length;
        for(var nMol=0; nMol<nMolCount; ++nMol)
        {
            this.m_arMols[nMol].delete();
        }

        this.m_arMols = null;
        this.m_viewFs = null;
        this.m_rdKitModule = null;

        super.dispose();
    }

    async init(strWebRoot, arSmiles, bMolFile)
    {
        const strFullPath = strWebRoot + "dist/6976ad1747eb075e7a1d.module.wasm";//"dist/93e50067519e22ffe954.module.wasm";//"api/packages/published/files/Testpackage/admin/dist/93e50067519e22ffe954.module.wasm";
        this.m_rdKitModule = await createRDKIt(strFullPath);

        this.m_bMolFile = bMolFile;

        const nMolCount = arSmiles.length;
        this.m_arMols = new Array(nMolCount);

        const viewFingerprints = new Int32Array(nMolCount*RDKitWorker.nSizeFingerPrint);
        this.m_viewFs = viewFingerprints;
    }

    continueInit(arSmiles, nIdxSmiles, bMolFile)
    {
        if(arSmiles === null)
            throw new Error("Array containing smiles cannot be null at mol index " + nIdxSmiles);

        const PROGRESS_MOL_BATCH_COUNT = 5;
        const PROGRESS_MIN_MOL_COUNT_IN_BATCH = 100;

        let PROGRESS_MOL_COUNT_IN_BATCH = Math.floor(arSmiles.length/PROGRESS_MOL_BATCH_COUNT);
        if(PROGRESS_MOL_COUNT_IN_BATCH < PROGRESS_MIN_MOL_COUNT_IN_BATCH)
            PROGRESS_MOL_COUNT_IN_BATCH = PROGRESS_MIN_MOL_COUNT_IN_BATCH;

        let strSmiles = "";
        let molecule = null;
        let strIndex = null;

        const viewFingerprints = this.m_viewFs;

        const nMolCount = arSmiles.length;
        const arIndex = new Array(RDKitWorker.nSizeFingerPrint);

        let nMol=nIdxSmiles;

        try
        {
            for(let nMolTmp = 0; nMol<nMolCount && nMolTmp < PROGRESS_MOL_COUNT_IN_BATCH; ++nMolTmp, ++nMol)
            {
                strSmiles = arSmiles[nMol];
                molecule = this.m_rdKitModule.get_mol(strSmiles);

                this.m_arMols[nMol] = molecule;
                strIndex = molecule.get_morgan_fp(2, RDKitWorker.FP_SIZE);
                RDKitWorker.str2ArrayFingerprint(arIndex, strIndex);

                for(let n=0; n<RDKitWorker.nSizeFingerPrint; ++n)
                {
                    viewFingerprints[nMol*RDKitWorker.nSizeFingerPrint + n] = arIndex[n];
                }
            }

         return nMol;
        }
        catch(e)
        {
            throw new Error("Error creating a molecule: " + strSmiles);
        }
    }

    tanimotoScores(strSmilesFrag)
    {
        const moltFragment = this.m_rdKitModule.get_mol(strSmilesFrag);
        const strIndex = moltFragment.get_morgan_fp(2, RDKitWorker.FP_SIZE);

        const arIndexFrag = new Array(RDKitWorker.nSizeFingerPrint);
        RDKitWorker.str2ArrayFingerprint(arIndexFrag, strIndex);

        const nMolCount = this.m_arMols.length;
        const arScores = new Array(nMolCount);//new Float32Array(nMolCount);
        let fScore = 0.0;
        let b = false;

        const arIndex = new Array(RDKitWorker.nSizeFingerPrint);
        for(var nM=0; nM<nMolCount; ++nM)
        {
            for(var n=0; n<RDKitWorker.nSizeFingerPrint; ++n)
            {
                arIndex[n] = this.m_viewFs[RDKitWorker.nSizeFingerPrint*nM +n];
            }

            fScore = RDKitWorker.getSimilarityTanimoto(arIndexFrag, arIndex);
            arScores[nM] = fScore;
        }

        moltFragment.delete();
        return arScores;
    }

    structureSearch(strSmilesFrag)
    {
        const moltFragment = this.m_rdKitModule.get_mol(strSmilesFrag);

        let b = false;
        const nMolCount = this.m_arMols.length;
        const arFlags = new Array(nMolCount);
        let mol = null;
        let match = null;
        for(var nM=0; nM<nMolCount; ++nM)
        {
            mol = this.m_arMols[nM];
            match = mol.get_substruct_match(moltFragment);
            b = match !== "{}";
            arFlags[nM] = b ? 1 : 0;
        }

        moltFragment.delete();
        return arFlags;
   }


    paintStructure(offscreen, nMol, nX, nY, nW, nH)
    {
        const mol = this.m_arMols[nMol];

        const g = offscreen.getContext('2d');

        g.save();
        mol.draw_to_canvas_with_offset(g.canvas, nX, -nY, nW, nH);
        g.restore();
    }
}


RDKitWorker.FP_SIZE = 512;//128;
RDKitWorker.nSizeFingerPrint = RDKitWorker.FP_SIZE/32;

RDKitWorker.str2ArrayFingerprint = function (arIndex, strIndex)
{
    let strNum = null;
    let nNum = -1;

    for(var n=0; n<RDKitWorker.nSizeFingerPrint; ++n)
    {
        strNum = strIndex.substring(n*32, (n+1)*32);
        nNum = parseInt(strNum, 2);
        arIndex[n] = nNum;
    }
}

RDKitWorker.bitCount_1 = function (x_0){
    // jl.$clinit_Integer();
    x_0 -= x_0 >> 1 & 1431655765;
    x_0 = (x_0 >> 2 & 858993459) + (x_0 & 858993459);
    x_0 = (x_0 >> 4) + x_0 & 252645135;
    x_0 += x_0 >> 8;
    x_0 += x_0 >> 16;
    return x_0 & 63;
}

RDKitWorker.getSimilarityTanimoto = function (nIndexOne, nIndexTwo)
{
    let nSharedKeys = 0;
    let nAllKeys = 0;
    for(let i = 0; i < nIndexOne.length; i++)
    {
     nSharedKeys += RDKitWorker.bitCount_1(nIndexOne[i] & nIndexTwo[i]);
     nAllKeys += RDKitWorker.bitCount_1(nIndexOne[i] | nIndexTwo[i]);
    }
    return nSharedKeys / nAllKeys;
}
