/**
 * {@link CruxPool} — a pool of Web Workers that shards a dataset across all CPU
 * cores, parses each shard in parallel, then broadcasts every query to all
 * shards and merges the results. Built for the search-as-you-draw pattern: load
 * a dataset once, query it constantly with cancellation, progress and streaming.
 *
 *   await using pool = await CruxPool.create();
 *   await using ds   = await pool.load(smiles, { onProgress });
 *   const live = ds.session({ debounceMs: 50 });
 *   for await (const b of live.searchStream({ kind: "substructure", smarts }))
 *     render(b.indices);          // global input positions; map to your data
 *
 * Results are typed arrays of original input positions (decision: exact index),
 * transferred zero-copy from the workers. SMILES are not echoed back — map the
 * indices to the array you passed to {@link CruxPool.load}.
 */
import { initCrux } from "./wasm.js";
import { mergeSubstructure, mergeSimilarity } from "./merge.js";
import { planShards } from "./sharding.js";
import { InlineBackend, WorkerBackend } from "./backend.js";
// ------------------------------------------------------------------- helpers
function hardwareConcurrency() {
    const n = globalThis.navigator
        ?.hardwareConcurrency;
    return typeof n === "number" && n > 0 ? Math.floor(n) : 4;
}
function abortError() {
    return new DOMException("The operation was aborted.", "AbortError");
}
async function compileModule(wasm) {
    if (wasm === undefined)
        wasm = new URL("./wasm/crux_wasm_bg.wasm", import.meta.url);
    if (wasm instanceof WebAssembly.Module)
        return wasm;
    if (wasm instanceof Response)
        return WebAssembly.compileStreaming(wasm);
    if (wasm instanceof URL || typeof wasm === "string") {
        const url = wasm;
        try {
            return await WebAssembly.compileStreaming(fetch(url));
        }
        catch {
            return WebAssembly.compile(await (await fetch(url)).arrayBuffer());
        }
    }
    return WebAssembly.compile(wasm);
}
/** Resolve after `ms`, or reject with the signal's reason if it aborts first. */
function delay(ms, signal) {
    return new Promise((resolve, reject) => {
        if (signal.aborted)
            return reject(signal.reason ?? abortError());
        const onAbort = () => {
            clearTimeout(timer);
            reject(signal.reason ?? abortError());
        };
        const timer = setTimeout(() => {
            signal.removeEventListener("abort", onAbort);
            resolve();
        }, ms);
        signal.addEventListener("abort", onAbort, { once: true });
    });
}
/** A push/pull async queue. `break`ing the consumer invokes `onReturn`. */
class AsyncQueue {
    #values = [];
    #waiters = [];
    #ended = false;
    #error;
    #hasError = false;
    #onReturn;
    constructor(onReturn) {
        this.#onReturn = onReturn;
    }
    push(v) {
        if (this.#ended)
            return;
        const w = this.#waiters.shift();
        if (w)
            w.resolve({ value: v, done: false });
        else
            this.#values.push(v);
    }
    end() {
        if (this.#ended)
            return;
        this.#ended = true;
        for (const w of this.#waiters.splice(0))
            w.resolve({ value: undefined, done: true });
    }
    fail(err) {
        if (this.#ended)
            return;
        this.#ended = true;
        this.#hasError = true;
        this.#error = err;
        for (const w of this.#waiters.splice(0))
            w.reject(err);
    }
    [Symbol.asyncIterator]() {
        return {
            next: () => {
                if (this.#values.length)
                    return Promise.resolve({ value: this.#values.shift(), done: false });
                if (this.#hasError) {
                    this.#hasError = false;
                    return Promise.reject(this.#error);
                }
                if (this.#ended)
                    return Promise.resolve({ value: undefined, done: true });
                return new Promise((resolve, reject) => this.#waiters.push({ resolve, reject }));
            },
            return: (value) => {
                this.#onReturn?.();
                this.#ended = true;
                return Promise.resolve({ value: value, done: true });
            },
        };
    }
}
// --------------------------------------------------------------------- pool
export class CruxPool {
    #backend;
    #laneCount;
    #active = new Map();
    #nextQueryId = 1;
    #shardCount = 0;
    #disposed = false;
    #onReady;
    #onReadyError;
    #onLoadMessage;
    constructor(laneCount) {
        this.#laneCount = laneCount;
    }
    /** Number of live worker lanes. */
    get workerCount() {
        return this.#laneCount;
    }
    /**
     * Spawn the workers, compile the WASM once on the main thread and share the
     * compiled module to every worker. Resolves when all workers are ready. In a
     * Node / no-Worker environment, falls back to a single-thread inline backend.
     */
    static async create(opts = {}) {
        const workers = Math.max(1, Math.floor(opts.workers ?? hardwareConcurrency()));
        if (typeof Worker !== "undefined") {
            const module = await compileModule(opts.wasm);
            const pool = new CruxPool(workers);
            pool.#backend = new WorkerBackend(workers, opts.workerUrl, (id, msg) => pool.#handleLaneMessage(id, msg), (id, err) => pool.#handleLaneError(id, err));
            await pool.#initWorkers(module);
            return pool;
        }
        // Single-thread fallback: initialize WASM on this thread once.
        await initCrux(opts.wasm ?? new URL("./wasm/crux_wasm_bg.wasm", import.meta.url));
        const pool = new CruxPool(workers);
        pool.#backend = new InlineBackend(workers, (id, msg) => pool.#handleLaneMessage(id, msg));
        return pool;
    }
    #initWorkers(module) {
        const resolvers = new Map();
        const readys = [];
        return new Promise((resolve, reject) => {
            this.#onReadyError = (err) => reject(err);
            this.#onReady = (laneId) => resolvers.get(laneId)?.();
            for (let i = 0; i < this.#laneCount; i++) {
                readys.push(new Promise((res) => resolvers.set(i, res)));
                this.#backend.lanes[i].post({ kind: "init", module, laneId: i });
            }
            Promise.all(readys).then(() => resolve());
        }).finally(() => {
            this.#onReady = undefined;
            this.#onReadyError = undefined;
        });
    }
    /** Shard `smiles` across the workers and parse in parallel. */
    async load(smiles, opts = {}) {
        if (this.#disposed)
            throw new Error("CruxPool is disposed");
        if (this.#shardCount > 0)
            await this._dropDataset(); // one dataset at a time
        const total = smiles.length;
        const plan = planShards(total, this.#laneCount, opts);
        const shardCount = plan.shardCount;
        const buildIndex = opts.buildIndex ?? false;
        const chunkSize = opts.chunkSize ?? 5000;
        this.#shardCount = shardCount;
        const perShardDone = new Array(shardCount).fill(0);
        const shardSizes = new Array(shardCount).fill(0);
        let shardsLoaded = 0;
        let resolveLoad;
        let rejectLoad;
        const done = new Promise((res, rej) => {
            resolveLoad = res;
            rejectLoad = rej;
        });
        const report = () => {
            let d = 0;
            for (const x of perShardDone)
                d += x;
            opts.onProgress?.({ done: d, total, shardsDone: shardsLoaded, shardCount });
        };
        this.#onLoadMessage = (msg) => {
            if (msg.kind === "shardProgress") {
                perShardDone[msg.shardId] = msg.done;
                report();
            }
            else if (msg.kind === "shardLoaded") {
                perShardDone[msg.shardId] = msg.total;
                shardSizes[msg.shardId] = msg.size;
                shardsLoaded++;
                report();
                if (shardsLoaded === shardCount)
                    resolveLoad();
            }
            else if (msg.kind === "error") {
                rejectLoad(new Error(msg.message));
            }
        };
        const onAbort = () => rejectLoad(opts.signal?.reason ?? abortError());
        if (opts.signal) {
            if (opts.signal.aborted) {
                this.#onLoadMessage = undefined;
                this.#shardCount = 0;
                throw opts.signal.reason ?? abortError();
            }
            opts.signal.addEventListener("abort", onAbort, { once: true });
        }
        const enc = new TextEncoder();
        for (let s = 0; s < shardCount; s++) {
            const { start, end } = plan.shards[s];
            const laneId = s % this.#laneCount;
            const bytes = enc.encode(smiles.slice(start, end).join("\n"));
            const payload = bytes.buffer;
            this.#backend.lanes[laneId].post({ kind: "loadShard", shardId: s, shardOffset: start, payload, buildIndex, chunkSize }, [payload]);
        }
        try {
            await done;
        }
        catch (e) {
            this.#broadcast({ kind: "dropAll" });
            this.#shardCount = 0;
            throw e;
        }
        finally {
            if (opts.signal)
                opts.signal.removeEventListener("abort", onAbort);
            this.#onLoadMessage = undefined;
        }
        const size = shardSizes.reduce((a, b) => a + b, 0);
        return new CruxDataset(this, { size, failed: total - size, indexed: buildIndex, shardCount });
    }
    /** Terminate all workers and free their WASM memory. Idempotent. */
    async dispose() {
        if (this.#disposed)
            return;
        this.#disposed = true;
        this.#settleAll();
        this.#backend.dispose();
        this.#shardCount = 0;
    }
    [Symbol.asyncDispose]() {
        return this.dispose();
    }
    // ---- internal: used by CruxDataset (same module) ----
    /** @internal */
    _startSubstructure(smarts, opts) {
        const limit = opts.limit ?? 0;
        return this.#startSearch((queryId) => ({ kind: "searchSub", queryId, smarts, limit }), opts);
    }
    /** @internal */
    _startSimilarity(smiles, opts) {
        const limit = opts.limit ?? 0;
        const threshold = opts.threshold ?? 0;
        return this.#startSearch((queryId) => ({ kind: "searchSim", queryId, smiles, threshold, limit }), opts);
    }
    /** @internal Free the current dataset's shards and cancel its searches. */
    async _dropDataset() {
        this.#settleAll();
        this.#broadcast({ kind: "dropAll" });
        this.#shardCount = 0;
    }
    // ---- private orchestration ----
    #startSearch(makeMsg, opts) {
        if (this.#disposed)
            throw new Error("CruxPool is disposed");
        const queryId = this.#nextQueryId++;
        const op = {
            queryId,
            lanesOutstanding: new Set(Array.from({ length: this.#laneCount }, (_, i) => i)),
            shardsDone: 0,
            shardCount: this.#shardCount,
            hits: 0,
            queue: new AsyncQueue(() => this.#stopOp(op)),
            onProgress: opts.onProgress,
            signal: opts.signal,
            settled: false,
        };
        if (opts.signal?.aborted) {
            op.settled = true;
            op.queue.fail(opts.signal.reason ?? abortError());
            return op.queue;
        }
        this.#active.set(queryId, op);
        if (opts.signal) {
            op.abortListener = () => this.#abortOp(op, opts.signal.reason ?? abortError());
            opts.signal.addEventListener("abort", op.abortListener, { once: true });
        }
        this.#broadcast(makeMsg(queryId));
        return op.queue;
    }
    #broadcast(msg) {
        for (const lane of this.#backend.lanes)
            lane.post(msg);
    }
    #handleLaneMessage(laneId, msg) {
        switch (msg.kind) {
            case "ready":
                this.#onReady?.(laneId);
                return;
            case "shardProgress":
            case "shardLoaded":
                this.#onLoadMessage?.(msg);
                return;
            case "subPartial": {
                const op = this.#active.get(msg.queryId);
                if (!op)
                    return;
                op.shardsDone++;
                op.hits += msg.indices.length;
                op.onProgress?.({ shardsDone: op.shardsDone, shardCount: op.shardCount, hits: op.hits });
                op.queue.push({ indices: msg.indices, shardsDone: op.shardsDone, shardCount: op.shardCount });
                return;
            }
            case "simPartial": {
                const op = this.#active.get(msg.queryId);
                if (!op)
                    return;
                op.shardsDone++;
                op.hits += msg.indices.length;
                op.onProgress?.({ shardsDone: op.shardsDone, shardCount: op.shardCount, hits: op.hits });
                op.queue.push({
                    indices: msg.indices,
                    scores: msg.scores,
                    shardsDone: op.shardsDone,
                    shardCount: op.shardCount,
                });
                return;
            }
            case "searchDone": {
                const op = this.#active.get(msg.queryId);
                if (!op)
                    return;
                op.lanesOutstanding.delete(laneId);
                if (op.lanesOutstanding.size === 0)
                    this.#finishOp(op);
                return;
            }
            case "error": {
                if (msg.queryId !== undefined) {
                    const op = this.#active.get(msg.queryId);
                    if (op) {
                        this.#detach(op);
                        op.queue.fail(new Error(msg.message));
                    }
                }
                else {
                    // No query id: a load or init error. Surface to whichever is pending.
                    this.#onReadyError?.(new Error(msg.message));
                    this.#onLoadMessage?.(msg);
                }
                return;
            }
            case "disposed":
                return;
        }
    }
    #handleLaneError(_laneId, err) {
        // Fail-fast (v1): surface to any pending create, load, and searches.
        this.#onReadyError?.(err);
        this.#onLoadMessage?.({ kind: "error", message: err.message });
        for (const op of [...this.#active.values()]) {
            this.#detach(op);
            op.queue.fail(err);
        }
    }
    #finishOp(op) {
        if (op.settled)
            return;
        this.#detach(op);
        op.queue.end();
    }
    #abortOp(op, reason) {
        if (op.settled)
            return;
        this.#detach(op);
        this.#broadcast({ kind: "cancel", queryId: op.queryId });
        op.queue.fail(reason);
    }
    #stopOp(op) {
        // Consumer broke out of the stream: cancel workers, end gracefully.
        if (op.settled)
            return;
        this.#detach(op);
        this.#broadcast({ kind: "cancel", queryId: op.queryId });
        op.queue.end();
    }
    #detach(op) {
        op.settled = true;
        this.#active.delete(op.queryId);
        if (op.signal && op.abortListener)
            op.signal.removeEventListener("abort", op.abortListener);
    }
    #settleAll() {
        for (const op of [...this.#active.values()])
            this.#stopOp(op);
    }
}
// ------------------------------------------------------------------ dataset
export class CruxDataset {
    /** Molecules successfully loaded (parse failures excluded). */
    size;
    /** Input SMILES that failed to parse. */
    failed;
    indexed;
    shardCount;
    #pool;
    #disposed = false;
    /** @internal Created by {@link CruxPool.load}. */
    constructor(pool, meta) {
        this.#pool = pool;
        this.size = meta.size;
        this.failed = meta.failed;
        this.indexed = meta.indexed;
        this.shardCount = meta.shardCount;
    }
    /** Streaming substructure (SMARTS) search — one batch per shard as it finishes. */
    substructureSearchStream(smarts, opts = {}) {
        this.#check();
        return this.#pool._startSubstructure(smarts, opts);
    }
    /** Collected substructure search — merged, sorted by ascending index, capped at `limit`. */
    async substructureSearch(smarts, opts = {}) {
        const parts = [];
        for await (const b of this.substructureSearchStream(smarts, opts))
            parts.push(b.indices);
        return mergeSubstructure(parts, opts.limit ?? 0);
    }
    /** Streaming similarity search — one batch per shard as it finishes. */
    similaritySearchStream(querySmiles, opts = {}) {
        this.#check();
        return this.#pool._startSimilarity(querySmiles, opts);
    }
    /** Collected similarity search — globally merged + re-sorted desc, capped at `limit`. */
    async similaritySearch(querySmiles, opts = {}) {
        const parts = [];
        for await (const b of this.similaritySearchStream(querySmiles, opts))
            parts.push({ indices: b.indices, scores: b.scores });
        return mergeSimilarity(parts, opts.limit ?? 0);
    }
    /** Create a latest-wins session for search-as-you-draw (auto-cancels the prior query). */
    session(opts = {}) {
        this.#check();
        return new SearchSession(this, opts);
    }
    /** Free this dataset's shard collections (the pool stays alive). Idempotent. */
    async dispose() {
        if (this.#disposed)
            return;
        this.#disposed = true;
        await this.#pool._dropDataset();
    }
    [Symbol.asyncDispose]() {
        return this.dispose();
    }
    #check() {
        if (this.#disposed)
            throw new Error("CruxDataset is disposed");
    }
}
// ------------------------------------------------------------------ session
export class SearchSession {
    #dataset;
    #debounceMs;
    #controller;
    /** @internal Created by {@link CruxDataset.session}. */
    constructor(dataset, opts = {}) {
        this.#dataset = dataset;
        this.#debounceMs = opts.debounceMs ?? 0;
    }
    /** Run a query, auto-cancelling the previous one. Resolves to merged results. */
    async search(query) {
        const signal = this.#begin();
        if (this.#debounceMs > 0)
            await delay(this.#debounceMs, signal);
        if (query.kind === "substructure")
            return this.#dataset.substructureSearch(query.smarts, { ...query.options, signal });
        return this.#dataset.similaritySearch(query.smiles, { ...query.options, signal });
    }
    /** Stream a query, auto-cancelling the previous one. */
    searchStream(query) {
        const signal = this.#begin();
        const dataset = this.#dataset;
        const debounceMs = this.#debounceMs;
        return (async function* () {
            if (debounceMs > 0)
                await delay(debounceMs, signal);
            const stream = query.kind === "substructure"
                ? dataset.substructureSearchStream(query.smarts, { ...query.options, signal })
                : dataset.similaritySearchStream(query.smiles, { ...query.options, signal });
            yield* stream;
        })();
    }
    /** Cancel the current query without starting a new one. */
    cancel() {
        this.#controller?.abort();
        this.#controller = undefined;
    }
    dispose() {
        this.cancel();
    }
    [Symbol.dispose]() {
        this.dispose();
    }
    #begin() {
        this.#controller?.abort();
        const controller = new AbortController();
        this.#controller = controller;
        return controller.signal;
    }
}
//# sourceMappingURL=pool.js.map