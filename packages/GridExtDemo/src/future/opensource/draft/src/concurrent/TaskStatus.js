export class TaskStatus {

    constructor(strName, strTitle) {
        this.m_strName = strName;
        this.m_strTitle = strTitle;
    }

    getTitle() {
        return this.m_strTitle;
    }

    toString() {
        return this.m_strName;
    }
}

TaskStatus.Created = new TaskStatus('Created', 'Task Created');
TaskStatus.CancelRequested = new TaskStatus('Cancel Requested', 'Cancel Requested');
TaskStatus.Cancelled = new TaskStatus('Cancelled','Cancelled');
TaskStatus.Faulted = new TaskStatus('Faulted', 'Error');
TaskStatus.Running = new TaskStatus('Running','Running');
TaskStatus.RanSuccesfully = new TaskStatus('RanSuccesfully', 'Success');
