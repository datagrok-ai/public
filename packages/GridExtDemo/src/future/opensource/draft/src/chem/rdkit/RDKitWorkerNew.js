import {createRDKIt} from "./RDKit_minimal_2021.03_02.js";
import {MolWorkerNew} from "../MolWorkerNew";
import {AsyncUtils} from "../../concurrent/AsyncUtils";

export class RDKitWorkerNew extends MolWorkerNew
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
            if(this.m_arMols[nMol] !== undefined && this.m_arMols[nMol] !== null)
             this.m_arMols[nMol].delete();
        }

        this.m_arMols = null;
        this.m_viewFs = null;
        this.m_rdKitModule = null;

        super.dispose();
    }

    async init(strWebRoot, arSmiles, bMolFile, scope)
    {
        const strFullPath = strWebRoot + "dist/6976ad1747eb075e7a1d.module.wasm";//"dist/93e50067519e22ffe954.module.wasm";//"api/packages/published/files/Testpackage/admin/dist/93e50067519e22ffe954.module.wasm";
        this.m_rdKitModule = await createRDKIt(strFullPath);

        this.m_bMolFile = bMolFile;

        const nMolCount = arSmiles.length;
        this.m_arMols = new Array(nMolCount);

        const viewFingerprints = new Int32Array(nMolCount*RDKitWorkerNew.nSizeFingerPrint);
        this.m_viewFs = viewFingerprints;

        this.notifyContinueInitProgress(0, scope);
        this.continueInit(this.m_arTmpSmiles, this.m_nIdxMol, this.m_bMolFile, scope);
    }

    async continueInit(arSmiles, nIdxSmiles, bMolFile, scope)
    {

        if(arSmiles === null)
            throw new Error("Array containing smiles cannot be null at mol index " + nIdxSmiles);

        console.log("Proc CPU: " + this.m_nCPU + " smiles: " + nIdxSmiles + " " + this.m_bCancelled);
        if(this.m_bCancelled)
        {
            console.log("Breaking on Cancel Proc CPU: " + this.m_nCPU + " smiles: " + nIdxSmiles + " " + this.m_bCancelled);
            this.notifyContinueInitProgress(nIdxSmiles, scope);
            this.notifyInitDone(nIdxSmiles, scope);
            return;// nIdxSmiles;
        }

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
        const arIndex = new Array(RDKitWorkerNew.nSizeFingerPrint);

        let nMol=nIdxSmiles;

        for(let nMolTmp = 0; nMol<nMolCount && nMolTmp < PROGRESS_MOL_COUNT_IN_BATCH; ++nMolTmp, ++nMol)
            {
                try
                {

                strSmiles = arSmiles[nMol];
                molecule = this.m_rdKitModule.get_mol(strSmiles);

                this.m_arMols[nMol] = molecule;
                strIndex = molecule.get_morgan_fp(2, RDKitWorkerNew.FP_SIZE);
                RDKitWorkerNew.str2ArrayFingerprint(arIndex, strIndex);

                for(let n=0; n<RDKitWorkerNew.nSizeFingerPrint; ++n)
                {
                    viewFingerprints[nMol*RDKitWorkerNew.nSizeFingerPrint + n] = arIndex[n];
                }
               }
               catch(e)
               {

                setTimeout(function() { throw new Error("Error creating a molecule: " + strSmiles + " " + nMol + "  " + nMolTmp + "  " + e.message) });
                throw new Error("Error creating a molecule: " + strSmiles + " " + nMol + "  " + nMolTmp);
               }
            }

        this.m_nIdxMol = nMol;
        if(nMol !== nMolCount)
        {
            this.notifyContinueInitProgress(nMol, scope);
            await AsyncUtils.sleep(10);
            this.continueInit(arSmiles, nMol, bMolFile, scope);
        }
        else {
            this.notifyContinueInitProgress(nMol, scope);
            this.notifyInitDone(nMol, scope);
        }
         return nMol;
    }

    tanimotoScores(strSmilesFrag)
    {
        const moltFragment = this.m_rdKitModule.get_mol(strSmilesFrag);
        const strIndex = moltFragment.get_morgan_fp(2, RDKitWorkerNew.FP_SIZE);

        const arIndexFrag = new Array(RDKitWorkerNew.nSizeFingerPrint);
        RDKitWorkerNew.str2ArrayFingerprint(arIndexFrag, strIndex);

        const nMolCount = this.m_nIdxMol;
        const arScores = new Array(nMolCount);//new Float32Array(nMolCount);
        let fScore = 0.0;
        let b = false;

        const arIndex = new Array(RDKitWorkerNew.nSizeFingerPrint);
        for(var nM=0; nM<nMolCount; ++nM)
        {
            for(var n=0; n<RDKitWorkerNew.nSizeFingerPrint; ++n)
            {
                arIndex[n] = this.m_viewFs[RDKitWorkerNew.nSizeFingerPrint*nM +n];
            }

            fScore = RDKitWorkerNew.getSimilarityTanimoto(arIndexFrag, arIndex);
            arScores[nM] = fScore;
        }

        moltFragment.delete();
        return arScores;
    }

    structureSearch(strSmilesFrag)
    {
        const moltFragment = this.m_rdKitModule.get_mol(strSmilesFrag);

        let b = false;
        const nMolCount = this.m_nIdxMol;
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


RDKitWorkerNew.FP_SIZE = 512;//128;
RDKitWorkerNew.nSizeFingerPrint = RDKitWorkerNew.FP_SIZE/32;

RDKitWorkerNew.str2ArrayFingerprint = function (arIndex, strIndex)
{
    let strNum = null;
    let nNum = -1;

    for(let n=0; n<RDKitWorkerNew.nSizeFingerPrint; ++n)
    {
        strNum = strIndex.substring(n*32, (n+1)*32);
        nNum = parseInt(strNum, 2);
        arIndex[n] = nNum;
    }
}

RDKitWorkerNew.bitCount_1 = function (x_0){
    // jl.$clinit_Integer();
    x_0 -= x_0 >> 1 & 1431655765;
    x_0 = (x_0 >> 2 & 858993459) + (x_0 & 858993459);
    x_0 = (x_0 >> 4) + x_0 & 252645135;
    x_0 += x_0 >> 8;
    x_0 += x_0 >> 16;
    return x_0 & 63;
}

RDKitWorkerNew.getSimilarityTanimoto = function (nIndexOne, nIndexTwo)
{
    let nSharedKeys = 0;
    let nAllKeys = 0;
    for(let i = 0; i < nIndexOne.length; i++)
    {
     nSharedKeys += RDKitWorkerNew.bitCount_1(nIndexOne[i] & nIndexTwo[i]);
     nAllKeys += RDKitWorkerNew.bitCount_1(nIndexOne[i] | nIndexTwo[i]);
    }
    return nSharedKeys / nAllKeys;
}
