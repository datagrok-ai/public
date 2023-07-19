export class LockedEntity<T extends any> {
    private _value: T;
    private _locked: boolean = false;
    constructor(value: T) {
        this._value = value;
    }

    public lock() {
        this._locked = true;
    }

    public release() {
        this._locked = false;
    }

    public get value(): T {
        return this._value;
    }

    public set value(value: T) {
        this._value = value;
    }

    public async unlockPromise() {
        return new Promise<void>((resolve) => {
            const check = () => {
                if (!this._locked) {
                    resolve();
                } else {
                    setTimeout(check, 0);
                }
            }
            check();
        });
    }


}
