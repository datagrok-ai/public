const bundleUrl: string = (typeof document !== 'undefined' && document.currentScript) ?
  (document.currentScript as HTMLScriptElement).src : self.location.href;

const CRUX_BASE = new URL('crux/', bundleUrl).href;
const WORKER_URL = new URL('crux-module-worker.js', CRUX_BASE).href;

type Pending = {resolve: (value: any) => void, reject: (reason: Error) => void};

const MAX_LIVE_MODELS = 32;

export type CruxResource = {name: string, id: string, bytes: Uint8Array};

const POOL_SIZE = Math.max(1, Math.min(4, Math.floor((navigator.hardwareConcurrency ?? 2) / 2)));

/** One Crux wasm module in its own module worker (the protocol of crux-core `crux-modules.mjs`) */
export class CruxWorker {
  private worker: Worker | null = null;
  private ready: Promise<void> | null = null;
  private next = 0;
  private pending = new Map<number, Pending>();
  /** Resource name -> id of the bytes registered in the live worker */
  private readonly resources = new Map<string, string>();
  /** Crux memoises every registered model and never frees it: the worker restarts after this many */
  private registered = 0;
  /** Calls not settled yet, counted from the call itself: a call is posted only after the module loads */
  private calls = 0;

  constructor(readonly name: string) {}

  get pendingCount(): number {
    return this.calls;
  }

  hasResource(resource: CruxResource): boolean {
    return this.resources.get(resource.name) === resource.id;
  }

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
    if (this.registered >= MAX_LIVE_MODELS && this.calls === 0)
      this.terminate();
    ++this.calls;
    try {
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
    } finally {
      --this.calls;
    }
  }

  async fnDescribe(name: string): Promise<string> {
    ++this.calls;
    try {
      await this.init();
      return await this.post('fnDescribe', {name});
    } finally {
      --this.calls;
    }
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

/** A pool of up to `size` module workers (each its own wasm instance): a call goes to a worker that holds its
 * resource, else to an idle one, else to a new one while under the size, else to the least loaded */
export class CruxModule {
  private workers: CruxWorker[] = [];

  constructor(readonly name: string, readonly size = POOL_SIZE) {}

  get workerCount(): number {
    return this.workers.length;
  }

  private pick(resource?: CruxResource): CruxWorker {
    const byLoad = (list: CruxWorker[]) => [...list].sort((a, b) => a.pendingCount - b.pendingCount)[0];
    const holder = resource === undefined ? undefined : byLoad(this.workers.filter((w) => w.hasResource(resource)));
    if (holder !== undefined)
      return holder;
    const idle = this.workers.find((w) => w.pendingCount === 0);
    if (idle !== undefined)
      return idle;
    if (this.workers.length < this.size) {
      const worker = new CruxWorker(this.name);
      this.workers.push(worker);
      return worker;
    }
    return byLoad(this.workers);
  }

  callFn(name: string, ipc: Uint8Array, params: object, resource?: CruxResource): Promise<Uint8Array> {
    return this.pick(resource).callFn(name, ipc, params, resource);
  }

  fnDescribe(name: string): Promise<string> {
    return this.pick().fnDescribe(name);
  }

  terminate(): void {
    for (const w of this.workers)
      w.terminate();
    this.workers = [];
  }
}

export const cruxMl = new CruxModule('eda');
