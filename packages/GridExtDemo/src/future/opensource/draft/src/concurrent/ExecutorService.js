import {ProgressSerializer} from "./messaging/ProgressSerializer";
import {ProgressResponseMessage} from "./messaging/ProgressResponseMessage";
import {DisposeRequestMessage} from "./messaging/DisposeRequestMessage";
import {DisposeSerializer} from "./messaging/DisposeSerializer";
import {AsyncUtils} from "./AsyncUtils";

export class CancelArgs {
    constructor() {
        this.m_bCancel = false;
    }

    isCancelled()
    {
        return this.m_bCancel;
    }

    cancel()
    {
        this.m_bCancel = true;
    }
}


/**
 * An istance of ExecutorService class provide means to manage execution of a heavy job by spltting it into
 * one or more asynchronous tasks that run in parallel on multiple system CPUs. THe service allows to track progress.
 */
export class ExecutorService
{
    constructor(classWorker, bSkipBusy)
    {
        if(classWorker === undefined)
            throw new Error("Worker class cannot be undefined.");

        if(bSkipBusy === undefined)
            bSkipBusy = false;

        this.m_classWorker = classWorker;
        this.m_bSkipBusy = bSkipBusy;
        this.m_arDoneFlags = null;
        this.m_bCancelled = false;
    }

    createWorker()
    {
     return new this.m_classWorker();
    }

    isCancelled() {
        return this.m_bCancelled;
    }

    isBusy(nCPU)
    {
        if(this.m_arWorkers === undefined)
            return true;

        if(nCPU === undefined) {
            const nCPUCountToUse = this.m_arWorkers.length;
            for(nCPU = 0; nCPU < nCPUCountToUse; ++nCPU)
            {
                if(!this.m_arDoneFlags[nCPU])
                    return true;
            }
        }
        else
        {
            if(!this.m_arDoneFlags[nCPU])
                return true;
        }

        return false;
    }

    getCPU4Job(nJob){
        const nCPUCount = this.m_arWorkers.length;
        for(var nCPU=0; nCPU<nCPUCount; ++nCPU)
        {
            if(this.m_arFrs[nCPU] <= nJob && nJob <= this.m_arTos[nCPU])
                return nCPU;
        }
        return -1;
    }


    getJobIdxOnCPU4Job(nJob){
        const nCPUCount = this.m_arWorkers.length;
        for(var nCPU=0; nCPU<nCPUCount; ++nCPU)
        {
            if(this.m_arFrs[nCPU] <= nJob && nJob <= this.m_arTos[nCPU])
                return nJob - this.m_arFrs[nCPU];
        }
        return -1;
    }

    async dispose()
    {
        function CreatePostArgs(nCPU, nCPUCount)
        {
            let msg = new DisposeRequestMessage();
            let arArgs = DisposeSerializer.Instance.toRequestArgs(msg);

            return arArgs;
        }

        function OnResult(service, nCPU, nFr, nTo, obData)
        {
            let args = obData;
            let message = DisposeSerializer.Instance.toResponseMessage(args);

            service.m_arWorkers[nCPU].terminate();
            service.m_arWorkers[nCPU] = null;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {

        }

        await this.executeOnWorkers(CreatePostArgs, OnResult, OnError);
        this.m_arWorkers = null;
    }

    init(nJobCount, fnProgressCallback)
    {

        let nCPUCount = self.navigator.hardwareConcurrency;
        let nCPUCountToUse = nCPUCount < 2 || nJobCount <= nCPUCount ? 1 : nCPUCount;
        if(nCPUCountToUse > 10)
            nCPUCountToUse = 20;

        this.m_arFrs = new Array(nCPUCountToUse);
        this.m_arTos = new Array(nCPUCountToUse);

        let nFr = 0;
        let nTo = 0;
        let nCountPerCPU = Math.floor(nJobCount/nCPUCountToUse);
        if(nCountPerCPU === 0)
            nCountPerCPU = 1;

        this.m_arDoneFlags = Array(nCPUCountToUse);

        let nCPU=0;
        for(; nCPU<nCPUCountToUse; ++nCPU)
        {
            this.m_arDoneFlags[nCPU] = true;
        }

        this.m_arWorkers = new Array(nCPUCountToUse);
        let worker = null;

        let nTimeStart = new Date().getTime();
        for(nCPU=0; nCPU<nCPUCountToUse; ++nCPU)
        {
            nFr = nCountPerCPU*nCPU;
            nTo = nCPU === nCPUCountToUse-1 ? nJobCount -1 : nCountPerCPU*(nCPU+1) -1;
            this.m_arFrs[nCPU] = nFr;
            this.m_arTos[nCPU] = nTo;

            worker = this.createWorker();
            this.m_arWorkers[nCPU] = worker;

            worker.m_nId = nCPU;
            worker.m_nFr = nFr;
            worker.m_nTo = nTo;
        }
    }

    async executeOnWorkers(fnCreatePostArgs, fnOnResult, fnOnError, fnPnProgress)
    {
        if(this.isBusy()) {
            if(this.m_bSkipBusy)
                return;

            throw new Error("The service is busy.");
        }
        const service = this;
        let bSuccess = true;

        let arArgs = null;
        let arrArgs = [];

        let nTotalJobCountPrgOld = -1;
        let nTotalJobCountPrg = 0;
        const nCPUCountToUse = this.m_arWorkers.length;
        const arJobCountWorkerPrgs = new Array(nCPUCountToUse);
        arJobCountWorkerPrgs.fill(0);

        let nCPU=0;
        for(; nCPU<nCPUCountToUse; ++nCPU)
        {
            arArgs = fnCreatePostArgs(nCPU, nCPUCountToUse, this.m_arFrs[nCPU], this.m_arTos[nCPU]);
            arrArgs.push(arArgs);
            if(arArgs === null)
                continue;

            this.m_arDoneFlags[nCPU] = false;
        }

        async function onMessageWorker(event) {
            const worker = event.currentTarget;
            const nIdWorker = worker.m_nId;
            const nFr = worker.m_nFr;
            const nTo = worker.m_nTo;

            const obData = event.data;
            const msgPrg = ProgressSerializer.Instance.toRequestMessage(obData);
            const nProgress = msgPrg !== null ? msgPrg.getJobCount() : nTo - nFr +1;

            arJobCountWorkerPrgs[nIdWorker] = nProgress;
            nTotalJobCountPrg = 0;
            for(let n=0; n<arJobCountWorkerPrgs.length; ++n)
            {
                nTotalJobCountPrg += arJobCountWorkerPrgs[n];
            }

            const args = new CancelArgs();
            nTotalJobCountPrg = Math.floor(nTotalJobCountPrg);
            if(nTotalJobCountPrg !== nTotalJobCountPrgOld && !service.m_bCancelled) {
                if(fnPnProgress !== undefined) {
                    console.log("Progress worker: " + nIdWorker + "  " + nTotalJobCountPrg);
                    await fnPnProgress(nTotalJobCountPrg, args);
                }
             nTotalJobCountPrgOld = nTotalJobCountPrg;
            }

            if(!service.m_bCancelled)
             service.m_bCancelled = args.isCancelled();

            if(msgPrg !== null)
            {
                const bCancel = args.isCancelled();
                if(bCancel) {
                    let aaa = 0;
                }

                const msgRes = new ProgressResponseMessage(service.m_bCancelled);
                const arArgs = ProgressSerializer.Instance.toResponseArgs(msgRes);
                worker.postMessage(arArgs);
                return;
            }
            //console.log("Done Worker: " + nIdWorker)
            const obResult = await fnOnResult(service, nIdWorker, nFr, nTo, event.data);

            worker.onmessage = null;
            worker.onerror = null;
            service.m_arDoneFlags[nIdWorker] = true;
        }

        async function onErrorWorker(error) {

            let worker = event.currentTarget;
            let nIdWorker = worker.m_nId;
            let nFr = worker.m_nFr;
            let nTo = worker.m_nTo;

            bSuccess = false;

            await fnOnError(service, nIdWorker, nFr, nTo, error);

            worker.onmessage = null;
            worker.onerror = null;
            this.m_arDoneFlags[nIdWorker] = true;
        }


        let worker = null;
        for(nCPU=0; nCPU<nCPUCountToUse; ++nCPU)
        {
            arArgs = arrArgs[nCPU];
            if(arArgs === null)
                continue;

            worker = this.m_arWorkers[nCPU];
            worker.onmessage = onMessageWorker;
            worker.onerror = onErrorWorker;

            if(!(arArgs instanceof Array))
                arArgs = [arArgs];

            try {worker.postMessage(...arArgs);}
            catch (e)
            {
                throw e;
            }
        }

        for(nCPU=0; nCPU<nCPUCountToUse; ++nCPU)
           {
            while (!this.m_arDoneFlags[nCPU]) {
                await AsyncUtils.sleep(50);
           }
        }

        return bSuccess;
    }
}

