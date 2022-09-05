import { Subject } from 'rxjs';
import {ErrorUtils} from "../utils/ErrorUtils";
import {TaskStatus} from "./TaskStatus";

export class Task
{
    constructor(service, fnPostArgs, fnResultWorker, fnResult, fnError)
    {
        ErrorUtils.verifyNotUndefinedNotNull(service);

        this.m_service = service;
        this.m_status = TaskStatus.Created;
        this.m_obResult = null;

        this.m_fnPostArgs = fnPostArgs;
        this.m_fnResultWorker = fnResultWorker;
        this.m_fnResult = fnResult;
        this.m_fnError = fnError;

        this.m_onFinished = new Subject();
        this.m_onError = new Subject();
        this.m_onUpdate = new Subject();
    }

    get onUpdate() {return this.m_onUpdate;}
    get onError() {return this.m_onError;}
    get onFinished() {return this.m_onFinished;}

    getStatus() {return this.m_status;}
    getResult() {return this.m_obResult;}

    cancel()
    {
        if(this.m_service.isCancel())
           return;

        this.m_status = TaskStatus.CancelRequested;
        this.m_service.cancel();
    }

    isCancelled()
    {
        return this.m_status == TaskStatus.Cancelled;
    }

    isError()
    {
        return this.m_status == TaskStatus.Faulted;
    }

    isSuccess()
    {
        return this.m_status == TaskStatus.RanSuccesfully;
    }


    run()
    {
        if(this.m_status !== TaskStatus.Created)
          throw new Error('Task has already been executed: ' + this.m_status.toString());

        const task = this;
        function onResult(service, arResultsWorkers)
        {
           const obResult = task.m_fnResult(service, arResultsWorkers);

           if(task.m_status === TaskStatus.CancelRequested)
               task.m_status = TaskStatus.Cancelled;
           else if(task.m_status == TaskStatus.Running)
            task.m_status = TaskStatus.RanSuccesfully;
           //else whatever ciurrent status is

           task.m_obResult = obResult;
           task.onFinished.next(obResult);
        }

        function OnError(serviceThis, nIdWorker, nFr, nTo, event)
        {
          task.m_status = TaskStatus.Faulted;
          task.m_fnError(serviceThis, nIdWorker, nFr, nTo, event);
          task.onError.next(event);
        }

        this.m_status = TaskStatus.Running;
        this.m_service.executeOnWorkersAsync(this.m_fnPostArgs, this.m_fnResultWorker, onResult, OnError, this.onUpdate, this);
    }


    async awaitResult()
    {
        await this.m_service.awaitResult();
    }

}














