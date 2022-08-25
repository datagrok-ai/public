import {_package} from "../package.js";
import {InitSerializer} from "./messaging/InitSerializer";
import {InitRequestMessage} from "./messaging/InitRequestMessage";
import {ImageRequestMessage} from "./messaging/ImageRequestMessage";
import {ImageSerializer} from "./messaging/ImageSerializer";
import {TanimotoRequestMessage} from "./messaging/TanimotoRequestMessage";
import {TanimotoSerializer} from "./messaging/TanimotoSerializer";
import {StructSearchRequestMessage} from "./messaging/StructSearchRequestMessage";
import {StructSearchSerializer} from "./messaging/StructSearchSerializer";
import {ExecutorService} from "../concurrent/ExecutorService";

export class MolService
{
    constructor(engine)
    {
        if(engine === undefined)
            throw new Error("Engine object cannot be undefined.");

        this.m_engine = engine;
        this.m_service = new ExecutorService(engine.getWorkerClass(), false);
    }

    async init(arSmilesOrMolFIles, bMolFile, fnProgressCallback)
    {
        const engine = this.m_engine;
        let root = _package.webRoot;
        const nIdxFr = root.toString().indexOf("/api/");
        if(nIdxFr >= 0)
         root =  root.substring(nIdxFr+1);

        let bSuccess = true;


        this.m_arSmilesOrMolFiles = arSmilesOrMolFIles;
        this.m_bMolFile = bMolFile;

        const nMolCount = arSmilesOrMolFIles.length;
        this.m_nMolCount = nMolCount;

        this.m_service.init(nMolCount);

        //const service = this.m_service;
        function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
        {
            //let nFr = service.m_arFrs[nCPU];
            //let nTo = service.m_arTos[nCPU]
            const msg = new InitRequestMessage(root, arSmilesOrMolFIles, nJobFr, nJobTo, bMolFile);
            const arArgs = InitSerializer.Instance.toRequestArgs(msg);
            return arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = InitSerializer.Instance.toResponseMessage(args);
         let CCC = 0;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
            let CCC = 0;
        }

        const nTimeStart = new Date().getTime();

        await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError, fnProgressCallback);

        const nTimeEnd = new Date().getTime();
        console.log("Spent on Molecules Init: " + (nTimeEnd - nTimeStart)/1000);

        if(!bSuccess)
            throw new Error("Failed to initialize molecules.");
    }

     async dispose()
     {
         await this.m_service.dispose();
         this.m_service = null;
     }



    async createMolImage(nMol, nW, nH)
    {
        const service = this.m_service;
        let ctx = null;

        const nCPUMol = service.getCPU4Job(nMol);

       /*
        if (!service.m_arDoneFlags[nCPUMol]) {
         throw new Error("CPU " + nCPUMol + " is busy.")
        }*/


        function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
        {
         if(nCPU !== nCPUMol)
             return null;

            const nMolCPU = service.getJobIdxOnCPU4Job(nMol);// - this.m_arFrs[nCPUMol];
            const msg = new ImageRequestMessage(nMolCPU, nW, nH);
            const arArgs = ImageSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
          const args = obData;
          const message = ImageSerializer.Instance.toResponseMessage(args);

         ctx = message.getBitmapImage();
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
         let ccc = 0;
        }

        await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);

        return ctx;
    }


    async createMolImages(nW, nH) {

        const service = this;
        let ctx = null;


        function CreatePostArgs(nCPU, nCPUCount)
        {
            const msg = new ImageRequestMessage(-1, nW, nH);
            const arArgs = ImageSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = ImageSerializer.Instance.toResponseMessage(args);

            ctx = message.getBitmapImage();
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
            let ccc = 0;
        }

        await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);

        return ctx;
    }


    isMolFile() {return this.m_bMolFile;}
    getSmilesOrMolFile(nIdx)
    {
        return this.m_arSmilesOrMolFiles[nIdx];
    }


    async structureSearch(arFlags, strSmiles) {

        if(arFlags.length !== this.m_nMolCount)
            throw new Error("The size of the flags' array must be " + this.m_nMolCount + " but it is " + arFlags.length);

        arFlags.fill(false);

        function CreatePostArgs(nCPU, nCPUCount)
        {
            const msg = new StructSearchRequestMessage(strSmiles);
            const arArgs = StructSearchSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = StructSearchSerializer.Instance.toResponseMessage(args);
            const arFs = message.getFlags();
            const nMolCount = arFs.length;

            for(var n=nFr; n<=nTo; ++n)
            {
                arFlags[n] = arFs[n-nFr];
            }
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {

        }

        await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);
        let aasd = 0;
    }


    async similaritySearch(arFlags, strSmiles, fThresh) {

        if(arFlags.length !== this.m_nMolCount)
            throw new Error("The size of the flags' array must be " + this.m_nMolCount + " but it is " + arFlags.length);

        arFlags.fill(false);

        let arScores = new Array(arFlags.length);
        await this.tanimotoScores(arScores, strSmiles);

        for(var n=0; n<arFlags.length; ++n)
        {
            arFlags[n] = arScores[n] >= fThresh;
        }
    }

    async tanimotoScores(arScores, strSmiles)
    {

        if(arScores.length !== this.m_nMolCount)
            throw new Error("The size of the flags' array must be " + this.m_nMolCount + " but it is " + arFlags.length);

        arScores.fill(NaN);

        function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
        {
            const msg = new TanimotoRequestMessage(strSmiles);
            const arArgs = TanimotoSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = TanimotoSerializer.Instance.toResponseMessage(args);
            const arSs = message.getScores();
            const nMolCount = arSs.length;

            for(var n=nFr; n<=nTo; ++n) {
                arScores[n] = arSs[n-nFr];

                if(arScores[n] === null || 0.0 > arScores[n] || arScores[n] > 1.0)
                    throw new Error("Tanimoto score value " + arScores[n] + " is out of range [0,1] for record index " + n);
            }
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {

        }

        await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError);
    }
}




MolService.createMolService = async function(engine, arSmilesOrMolFIles, bMolFile, fnProgressCallback){
    const service = new MolService(engine);
    await service.init(arSmilesOrMolFIles, bMolFile, fnProgressCallback);
    return service;
}
