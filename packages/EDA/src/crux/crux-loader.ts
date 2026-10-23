const bundleUrl: string = (typeof document !== 'undefined' && document.currentScript) ?
  (document.currentScript as HTMLScriptElement).src : self.location.href;

const CRUX_BASE = new URL('crux/', bundleUrl).href;
const WORKER_URL = new URL('crux-module-worker.js', CRUX_BASE).href;

type Pending = {resolve: (value: any) => void, reject: (reason: Error) => void};

const MAX_LIVE_MODELS = 32;

export type CruxResource = {name: string, id: string, bytes: Uint8Array};

/** One Crux wasm module in its own module worker (the protocol of crux-core `crux-modules.mjs`) */
export class CruxModule {
  private worker: Worker | null = null;
  private ready: Promise<void> | null = null;
  private next = 0;
  private pending = new Map<number, Pending>();
  /** Resource name -> id of the bytes registered in the live worker */
  private readonly resources = new Map<string, string>();
  /** Crux memoises every registered model and never frees it: the worker restarts after this many */
  private registered = 0;

  constructor(readonly name: string) {}

  private init(): Promise<void> {
    if (this.ready !== null)
      return this.ready;

    const worker = new Worker(WORKER_URL, {type: 'module'});
    worker.onmessage = (e: MessageEvent) => {
      const {id, ok, result, error} = e.data;
      const p = this.pending.get(id);
      if (p === undefined)
        return;
      this.pending.delete(id);
      if (ok)
        p.resolve(result);
      else
        p.reject(new Error(error));
    };
    worker.onerror = (e: ErrorEvent) => {
      this.rejectAll(new Error(e.message || `Crux module "${this.name}": worker error`));
      this.terminate();
    };
    this.worker = worker;

    const base = new URL(`${this.name}/`, CRUX_BASE).href;
    this.ready = this.post('init', {js: new URL('crux_wasm.js', base).href}).then(() => {}, (e: Error) => {
      this.terminate();
      throw new Error(`Crux module "${this.name}" failed to load from ${base} (see crux/README.md): ${e.message}`);
    });
    return this.ready;
  }

  private post(op: string, args: object, transfer: Transferable[] = []): Promise<any> {
    const id = this.next++;
    return new Promise((resolve, reject) => {
      this.pending.set(id, {resolve, reject});
      this.worker!.postMessage({id, op, args}, transfer);
    });
  }

  private rejectAll(err: Error): void {
    for (const p of this.pending.values())
      p.reject(err);
    this.pending.clear();
  }

  /** Runs a Crux function on an Arrow IPC stream (transferred); `resource` is registered first unless live */
  async callFn(name: string, ipc: Uint8Array, params: object, resource?: CruxResource): Promise<Uint8Array> {
    if (this.registered >= MAX_LIVE_MODELS && this.pending.size === 0)
      this.terminate();
    await this.init();
    if (this.worker === null)
      throw new Error('Crux worker terminated');

    let reg;
    if (resource !== undefined && this.resources.get(resource.name) !== resource.id) {
      reg = this.post('fnResource', resource).catch((e: Error) => {
        if (this.resources.get(resource.name) === resource.id)
          this.resources.delete(resource.name);
        throw e;
      });
      this.resources.set(resource.name, resource.id);
      this.registered++;
    }
    const call = this.post('callFn', {name, ipc, params: JSON.stringify(params)}, [ipc.buffer]);

    const [out] = await Promise.all([call, reg]);
    return out;
  }

  async fnDescribe(name: string): Promise<string> {
    await this.init();
    return this.post('fnDescribe', {name});
  }

  terminate(): void {
    this.worker?.terminate();
    this.rejectAll(new Error('Crux worker terminated'));
    this.worker = null;
    this.ready = null;
    this.resources.clear();
    this.registered = 0;
  }
}

export const cruxMl = new CruxModule('eda');
