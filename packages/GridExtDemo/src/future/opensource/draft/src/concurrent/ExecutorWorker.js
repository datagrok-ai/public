export class ExecutorWorker
{
    constructor(fnDispose)
    {
        this.m_fnDispose = fnDispose;
    }

    dispose()
    {
        this.m_fnDispose();
    }

    onMessage(scope, worker, event)
    {
        const strMsg = event.data.msg;
    }
}
