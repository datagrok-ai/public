export class ProgressResponseMessage {
    constructor(bCancel) {
        this.m_bCancel = bCancel === undefined ? false : bCancel;
    }

    isCancel()
    {
        return this.m_bCancel;
    }
}
