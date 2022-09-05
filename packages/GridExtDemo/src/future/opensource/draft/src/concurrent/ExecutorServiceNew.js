import {ProgressSerializer} from "./messaging/ProgressSerializer";
import {ProgressResponseMessage} from "./messaging/ProgressResponseMessage";
import {DisposeRequestMessage} from "./messaging/DisposeRequestMessage";
import {DisposeSerializer} from "./messaging/DisposeSerializer";
import {AsyncUtils} from "./AsyncUtils";
import {Task} from "./Task";

/**
 * An istance of ExecutorService class provide means to manage execution of a heavy job by spltting it into
 * one or more asynchronous tasks that run in parallel on multiple system CPUs. THe service allows to track progress.
 */
export class ExecutorServiceNew
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
        this.m_bCancel = false;
    }

    createWorker()
    {
     return new this.m_classWorker();
    }

    isCancel() {
        return this.m_bCancel;
    }

    cancel()
    {
        if(this.m_bCancel)
            throw new Error("Operation has already been cancelled.");

        this.m_bCancel = true;
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

    dispose(fnResult)
    {
        function CreatePostArgs(nCPU, nCPUCount)
        {
            const msg = new DisposeRequestMessage();
            const arArgs = DisposeSerializer.Instance.toRequestArgs(msg);

            return arArgs;
        }

        function OnResultWorker(service, nCPU, nFr, nTo, obData)
        {
            const args = obData;
            const message = DisposeSerializer.Instance.toResponseMessage(args);

            service.m_arWorkers[nCPU].terminate();
            service.m_arWorkers[nCPU] = null;

            return true;
        }

        function OnResult(service, arResultsWorkers)
        {
            this.m_arWorkers = null;
            fnResult(service, arResultsWorkers);
            return true;
        }

        function OnError(service, nCPU, nFr, nTo, error)
        {
            return error;
        }

        const task = new Task(this, CreatePostArgs, OnResultWorker, OnResult, OnError);
        return task;
    }

    init(nJobCount)
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


    executeOnWorkersAsync(fnCreatePostArgs, fnOnResultWorker, fnOnResult, fnOnError, fnPnProgress, task)
    {
        if(this.isBusy()) {
            if(this.m_bSkipBusy)
                return;

            throw new Error("The service is busy.");
        }

        this.m_bCancel = false;
        let bCancelRequested = false;
        const serviceThis = this;

        let arArgs = null;
        let arrArgs = [];

        let nTotalJobCountPrgOld = -1;
        let nTotalJobCountPrg = 0;
        const nCPUCountToUse = this.m_arWorkers.length;
        const arJobCountWorkerPrgs = new Array(nCPUCountToUse);
        arJobCountWorkerPrgs.fill(0);

        const arResultsWorkers = new Array(nCPUCountToUse);

        let nCPU=0;
        for(; nCPU<nCPUCountToUse; ++nCPU)
        {
            arArgs = fnCreatePostArgs(nCPU, nCPUCountToUse, this.m_arFrs[nCPU], this.m_arTos[nCPU]);
            arrArgs.push(arArgs);
            if(arArgs === null)
                continue;

            this.m_arDoneFlags[nCPU] = false;
        }


        //async
        function onMessageWorker(event)
        {
            const worker = event.currentTarget;
            const nIdWorker = worker.m_nId;
            const nFr = worker.m_nFr;
            const nTo = worker.m_nTo;

            const obData = event.data;
            const msgPrg = ProgressSerializer.Instance.toRequestMessage(obData);
            const nProgress = msgPrg !== null ? msgPrg.getJobCount() : arJobCountWorkerPrgs[nIdWorker];//nTo - nFr +1;

            arJobCountWorkerPrgs[nIdWorker] = nProgress;
            nTotalJobCountPrg = 0;
            for(let n=0; n<arJobCountWorkerPrgs.length; ++n)
            {
                nTotalJobCountPrg += arJobCountWorkerPrgs[n];
            }

            nTotalJobCountPrg = Math.floor(nTotalJobCountPrg);
            if(nTotalJobCountPrg !== nTotalJobCountPrgOld /*&& !bCancelRequested*/)
            {
                if(fnPnProgress !== undefined && fnPnProgress !== null)
                {
                    console.log("Progress from worker: " + nIdWorker + "  " + nProgress);
                    /*await*/ //fnPnProgress(nTotalJobCountPrg, task);
                    fnPnProgress.next(nTotalJobCountPrg);
                }
                nTotalJobCountPrgOld = nTotalJobCountPrg;
            }


            if(msgPrg !== null) {
                const bCancel = serviceThis.isCancel();
                if(!bCancelRequested && bCancel)
                {
                    console.log("Sending Cancel: " + nIdWorker + " " + msgPrg.getJobCount());
                    let wrk = null;
                    for (let nW = 0; nW < serviceThis.m_arWorkers.length; ++nW)
                    {
                        wrk = serviceThis.m_arWorkers[nW];
                        const msgRes = new ProgressResponseMessage(bCancel);
                        const arArgs = ProgressSerializer.Instance.toResponseArgs(msgRes);
                        wrk.postMessage(arArgs);
                    }

                  bCancelRequested = bCancel;
                }

                return;
            }
            console.log("Done Worker: " + nIdWorker)
            const obResultWorker = /*await*/fnOnResultWorker(serviceThis, nIdWorker, nFr, nTo, event.data);
            arResultsWorkers[nIdWorker] = obResultWorker;

            serviceThis.m_arDoneFlags[nIdWorker] = true;

            if(!serviceThis.isBusy())
            {
                //await
                fnOnResult(serviceThis, arResultsWorkers);
            }

            worker.onmessage = null;
            worker.onerror = null;
        }

        async function onErrorWorker(event) {

            const worker = event.currentTarget;
            const nIdWorker = worker.m_nId;
            const nFr = worker.m_nFr;
            const nTo = worker.m_nTo;

            if(fnOnError !== undefined && fnOnError !== null)
             await fnOnError(serviceThis, nIdWorker, nFr, nTo, event);

            worker.onmessage = null;
            worker.onerror = null;
            serviceThis.m_arDoneFlags[nIdWorker] = true;
        }


        let worker = null;
        for(nCPU=0; nCPU<nCPUCountToUse; ++nCPU)
        {
            arArgs = arrArgs[nCPU];
            if(arArgs === null)
                continue;

            worker = this.m_arWorkers[nCPU]
            worker.onmessage = onMessageWorker;
            worker.onerror = onErrorWorker;

            if(!(arArgs instanceof Array))
                arArgs = [arArgs];

            try {worker.postMessage(...arArgs);}
            catch(e)
            {
                throw e;
            }
        }

    }


    async awaitResult() {

        const nCPUCountToUse =  this.m_arWorkers === null ? 0 : this.m_arWorkers.length;

        for(let nCPU=0; nCPU<nCPUCountToUse; ++nCPU)
        {
            while (!this.m_arDoneFlags[nCPU]) {
                await AsyncUtils.sleep(50);
            }
        }
    }
}


