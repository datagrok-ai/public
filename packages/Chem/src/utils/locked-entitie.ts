export class LockedEntity<T extends any> {
  private _value: T;
  private _lockPromise: Promise<void> | null = Promise.resolve();
  private _unlockFunction: () => void = () => { };
  constructor(value: T) {
    this._value = value;
  }

  public lock() {
    this._lockPromise = new Promise<void>((resolve) => {
      this._unlockFunction = resolve;
    });
  }

  public release() {
    this._unlockFunction();
  }

  public get value(): T {
    return this._value;
  }

  public set value(value: T) {
    this._value = value;
  }

  public async unlockPromise() {
    return this._lockPromise;
  }
}
