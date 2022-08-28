import {_package} from "../package.js";
import {InitSerializer} from "./messaging/InitSerializer";
import {InitRequestMessage} from "./messaging/InitRequestMessage";
import {ImageRequestMessage} from "./messaging/ImageRequestMessage";
import {ImageSerializer} from "./messaging/ImageSerializer";
import {TanimotoRequestMessage} from "./messaging/TanimotoRequestMessage";
import {TanimotoSerializer} from "./messaging/TanimotoSerializer";
import {StructSearchRequestMessage} from "./messaging/StructSearchRequestMessage";
import {StructSearchSerializer} from "./messaging/StructSearchSerializer";
import {ExecutorServiceNew} from "../concurrent/ExecutorServiceNew";
import {Task} from "../concurrent/Task";

export class MolServiceNew
{
    constructor(engine)
    {
        if(engine === undefined)
            throw new Error("Engine object cannot be undefined.");

        this.m_engine = engine;
        this.m_service = new ExecutorServiceNew(engine.getWorkerClass(), false);
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

        await this.m_service.init(nMolCount);

        //const service = this.m_service;
        function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
        {
            const msg = new InitRequestMessage(root, arSmilesOrMolFIles, nJobFr, nJobTo, bMolFile, nCPU);
            const arArgs = InitSerializer.Instance.toRequestArgs(msg);

            return arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = InitSerializer.Instance.toResponseMessage(args);

         return true;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
            return error;
        }

        const nTimeStart = new Date().getTime();
        const obResult = await this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError, fnProgressCallback);
        const nTimeEnd = new Date().getTime();
        console.log("Spent on Molecules Init: " + (nTimeEnd - nTimeStart)/1000);

        if(!bSuccess)
            throw new Error("Failed to initialize molecules.");

        return obResult;
    }


    addMolecules(arSmilesOrMolFIles, bMolFile)
    {
        //const service = this.m_service;
        function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
        {
            const msg = new InitRequestMessage(root, arSmilesOrMolFIles, nJobFr, nJobTo, bMolFile, nCPU);
            const arArgs = InitSerializer.Instance.toRequestArgs(msg);
            return arArgs;
        }

        function OnResultWorker(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = InitSerializer.Instance.toResponseMessage(args);
            const nCount = message.getMolCount();
            return nCount;
        }

        function OnResult(service, arResultsWorkers)
        {
            let nCountTotal = 0;
            for(let n=0; n<arResultsWorkers.length; ++n)
            {
                nCountTotal += arResultsWorkers[n];
            }

           return nCountTotal;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
            return error;
        }

       const task = new Task(this.m_service, CreatePostArgs, OnResultWorker, OnResult, OnError);

       const engine = this.m_engine;
       let root = _package.webRoot;
       const nIdxFr = root.toString().indexOf("/api/");
       if(nIdxFr >= 0)
          root =  root.substring(nIdxFr+1);

        this.m_arSmilesOrMolFiles = arSmilesOrMolFIles;
        this.m_bMolFile = bMolFile;

        const nMolCount = arSmilesOrMolFIles.length;
        this.m_nMolCount = nMolCount;

        this.m_service.init(nMolCount);
        //const nTimeStart = new Date().getTime();

        //const obResult = task.run();
        //this.m_service.executeOnWorkers(CreatePostArgs, OnResult, OnError, fnProgressCallback);
        //const nTimeEnd = new Date().getTime();
        //console.log("Spent on Molecules Init: " + (nTimeEnd - nTimeStart)/1000);

        //if(!bSuccess)
          //  throw new Error("Failed to initialize molecules.");

        return task;
    }


     dispose()
     {
         const serviceMol = this;
         function onResult(service, arResultsWorkers)
         {
             serviceMol.m_service = null;

             while(serviceMol.m_arSmilesOrMolFiles.length > 0)
             {
                 serviceMol.m_arSmilesOrMolFiles.pop();
             }

             serviceMol.m_arSmilesOrMolFiles = null;
         }

         const task = this.m_service.dispose(onResult);
         return task;
     }



    createMolImage(nMol, nW, nH)
    {
        const service = this.m_service;
        let ctx = null;

        const nCPUMol = service.getCPU4Job(nMol);

        function CreatePostArgs(nCPU, nCPUCount, nJobFr, nJobTo)
        {
         if(nCPU !== nCPUMol)
             return null;

            const nMolCPU = service.getJobIdxOnCPU4Job(nMol);// - this.m_arFrs[nCPUMol];
            const msg = new ImageRequestMessage(nMolCPU, nW, nH);
            const arArgs = ImageSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResultWorker(service, nCPU, nFr, nTo, obData)
        {
          const args = obData;
          const message = ImageSerializer.Instance.toResponseMessage(args);
          ctx = message.getBitmapImage();
          return ctx;
        }

        function OnResult(service)
        {
         return ctx;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
         return error;
        }

        const task = new Task(this.m_service, CreatePostArgs, OnResultWorker, OnResult, OnError);

        return task;
    }


    createMolImages(nW, nH)
    {
        //const service = this;
        const arImages = new Array(this.m_nMolCount);

        function CreatePostArgs(nCPU, nCPUCount)
        {
            const msg = new ImageRequestMessage(-1, nW, nH);
            const arArgs = ImageSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResultWorker(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = ImageSerializer.Instance.toResponseMessage(args);

            const ctx = message.getBitmapImage();
            for(var n=nFr; n<=nTo; ++n)
            {
                arImages[n] = ctx;
            }

            return null;
        }

        function OnResult(service, arResultsWorkers)
        {
          return arImages;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
           return error;
        }

        const task = new Task(service, CreatePostArgs, OnResultWorker, OnResult, OnError);
        return task;
    }


    isMolFile() {return this.m_bMolFile;}
    getSmilesOrMolFile(nIdx)
    {
        return this.m_arSmilesOrMolFiles[nIdx];
    }


    structureSearch(arFlags, strSmiles) {

        if(arFlags.length !== this.m_nMolCount)
            throw new Error("The size of the flags' array must be " + this.m_nMolCount + " but it is " + arFlags.length);

        arFlags.fill(false);

        function CreatePostArgs(nCPU, nCPUCount)
        {
            const msg = new StructSearchRequestMessage(strSmiles);
            const arArgs = StructSearchSerializer.Instance.toRequestArgs(msg);
            return  arArgs;
        }

        function OnResultWorker(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = StructSearchSerializer.Instance.toResponseMessage(args);
            const arFs = message.getFlags();
            const nMolCount = arFs.length;

            for(var n=nFr; n<=nTo; ++n)
            {
                arFlags[n] = arFs[n-nFr];
            }

            return null;
        }

        function OnResult(service, arResultsWorkers)
        {
         return arFlags;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
            return error;
        }

        const task = new Task(this.m_service, CreatePostArgs, OnResultWorker, OnResult, OnError);
        return task;
    }


    async similaritySearch(arFlags, strSmiles, fThresh)
    {
        if(arFlags.length !== this.m_nMolCount)
            throw new Error("The size of the flags' array must be " + this.m_nMolCount + " but it is " + arFlags.length);

        arFlags.fill(false);

        const arScores = new Array(arFlags.length);
        const task = this.tanimotoScores(arScores, strSmiles);
        task.run();
        await task.awaitResult();

        for(var n=0; n<arFlags.length; ++n)
        {
            arFlags[n] = arScores[n] >= fThresh;
        }
    }

    tanimotoScores(arScores, strSmiles)
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

        function OnResultWorker(service, nCPU, nFr, nTo, obData)
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

         return null;
        }

        function OnResult(service, arResultsWorkers)
        {
           return arScores;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
         return error;
        }

        const task = new Task(this.m_service, CreatePostArgs, OnResultWorker, OnResult, OnError);
        return task;
    }
}




MolServiceNew.createMolService = async function(engine, arSmilesOrMolFIles, bMolFile, fnProgressCallback)
{
    const service = new MolServiceNew(engine);
    await service.init(arSmilesOrMolFIles, bMolFile, fnProgressCallback);
    return service;
}

MolServiceNew.createMolServiceNew = async function(engine, arSmilesOrMolFIles, bMolFile, fnProgressCallback)
{
    const service = new MolServiceNew(engine);
    const task = service.addMolecules(arSmilesOrMolFIles, bMolFile);
    task.onUpdate = fnProgressCallback;

    task.run();

    await task.awaitResult();
    const obResult = task.getResult();

    return service;
}
