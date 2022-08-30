import {initOCLModule} from "./openchemlib-full.pretty.js";
import {MolWorker} from "../MolWorker";

export class OCLWorker extends MolWorker {
    constructor(fnDispose) {
        super(fnDispose);

        this.m_arMols = null;
        this.m_viewFs = null;

    }

    dispose() {
        this.m_arMols = null;
        this.m_viewFs = null;

        super.dispose();
    }

    async init(strWebRoot, arSmiles, bMolFile, fnProgressCallback)
    {
        const nMolCount = arSmiles.length;

        this.m_arMols = new Array(nMolCount);

        //let strSmiles = "";
        //let molecule = null;
        //let arIndex = [];
        var viewFingerprints = new Int32Array(nMolCount * OCLWorker.nSizeFingerPrint);
        this.m_viewFs = viewFingerprints;
    }

    continueInit(arSmiles, nIdxSmiles, bMolFile)
    {
        const PROGRESS_MOL_BATCH_COUNT = 5;
        const PROGRESS_MIN_MOL_COUNT_IN_BATCH = 100;

        let PROGRESS_MOL_COUNT_IN_BATCH = Math.floor(arSmiles.length/PROGRESS_MOL_BATCH_COUNT);
        if(PROGRESS_MOL_COUNT_IN_BATCH < PROGRESS_MIN_MOL_COUNT_IN_BATCH)
            PROGRESS_MOL_COUNT_IN_BATCH = PROGRESS_MIN_MOL_COUNT_IN_BATCH;


        let strSmiles = "";
        let molecule = null;
        let arIndex = [];
        var viewFingerprints = this.m_viewFs;

        let nMol=nIdxSmiles;
        let nMolCount = arSmiles.length;
        try {
            for(var nMolTmp = 0; nMol<nMolCount && nMolTmp < PROGRESS_MOL_COUNT_IN_BATCH; ++nMolTmp, ++nMol)//for (var n = 0; n < nMolCount; ++n)
            {
              strSmiles = arSmiles[nMol];

                if (bMolFile)
                    molecule = OCL.Molecule.fromMolfile(strSmiles);
                else
                    molecule = OCL.Molecule.fromSmiles(strSmiles);

                this.m_arMols[nMol] = molecule;
                arIndex = new OCLM.SSSearcherWithIndex().createIndexMy(molecule);

                for (var i = 0; i < arIndex.length; ++i)
                {
                    viewFingerprints[nMol * OCLWorker.nSizeFingerPrint + i] = arIndex[i];
                }
            }

            return nMol;
        }
        catch (e)
        {

            throw new Error("Error creating a molecule: " + strSmiles + " " + e.message);
        }
    }



    tanimotoScores(strSmilesFrag) {
        let moltFragment = this.m_bMolFile ? OCL.Molecule.fromMolfile(strSmilesFrag) : OCL.Molecule.fromSmiles(strSmilesFrag);
        //moltFragment.setFragment(true);
        let search = new OCLM.SSSearcherWithIndex();
        let arIndexFrag = search.createIndexMy(moltFragment);

        const nMolCount = this.m_arMols.length;
        let arScores = new Array(nMolCount);
        let fScore = 0.0;
        let b = false;
        const nSizeFingerPrint = 16;
        let arIndex = new Array(nSizeFingerPrint);
        for (var nM = 0; nM < nMolCount; ++nM) {
            for (var n = 0; n < nSizeFingerPrint; ++n) {
                arIndex[n] = this.m_viewFs[nSizeFingerPrint * nM + n];
            }

            fScore = OCLM.SSSearcherWithIndex.getSimilarityTanimoto(arIndexFrag, arIndex);
            arScores[nM] = fScore;
        }
        return arScores;
    }

    structureSearch(strSmilesFrag) {
        let moltFragment = this.m_bMolFile ? OCL.Molecule.fromMolfile(strSmilesFrag) : OCL.Molecule.fromSmiles(strSmilesFrag);
        moltFragment.setFragment(true);

        let search = new OCL.SSSearcher();
        search.setFragment(moltFragment);

        let b = false;
        const nMolCount = this.m_arMols.length;
        let arFlags = new Array(nMolCount);

        for (var nM = 0; nM < nMolCount; ++nM) {
            search.setMolecule(this.m_arMols[nM]);
            b = search.isFragmentInMolecule();
            arFlags[nM] = b ? 1 : 0;
        }

        return arFlags;
    }

    paintStructure(offscreen, nMol, nX, nY, nW, nH) {
        let mol = this.m_arMols[nMol];

        let g = offscreen.getContext('2d');

        g.save();
        OCL.StructureView.drawMolecule(offscreen, mol, undefined, nX, nY, nW, nH, g);
        g.restore();
    }
}

OCLWorker.nSizeFingerPrint = 16;
