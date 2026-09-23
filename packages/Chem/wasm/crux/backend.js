/**
 * The transport the pool talks to, abstracting over real Web Workers and a
 * single-thread inline fallback. The pool's orchestration (sharding, merge,
 * cancellation) is identical for both — only message delivery differs.
 *
 * - {@link WorkerBackend}: one `Worker` per lane (browser; true parallelism).
 * - {@link InlineBackend}: one in-process {@link WorkerCore} per lane (Node /
 *   no-Worker environments; correct but not parallel). Message delivery is
 *   emulated with macrotasks so the same between-shard cancellation checkpoints
 *   behave the same way.
 */
import { WorkerCore } from "./worker-core.js";
/** Run `fn` on a fresh macrotask (mirrors the async worker message boundary). */
function deliver(fn) {
    setTimeout(fn, 0);
}
/** Lanes backed by real Web Workers (browser). */
export class WorkerBackend {
    laneCount;
    lanes = [];
    #workers = [];
    constructor(laneCount, workerUrl, onMessage, onError) {
        this.laneCount = laneCount;
        for (let i = 0; i < laneCount; i++) {
            const laneId = i;
            // The literal `new URL("./pool.worker.js", import.meta.url)` form with a
            // STATIC options object lets Vite/webpack/Rollup find and bundle the
            // worker chunk (Vite rejects non-static worker options). The explicit
            // `workerUrl` is an escape hatch for unusual setups.
            const w = workerUrl
                ? new Worker(workerUrl, { type: "module", name: "crux-pool" })
                : new Worker(new URL("./pool.worker.js", import.meta.url), {
                    type: "module",
                    name: "crux-pool",
                });
            w.onmessage = (e) => onMessage(laneId, e.data);
            w.onerror = (e) => {
                const where = e.filename ? ` @ ${e.filename}:${e.lineno}:${e.colno}` : "";
                const detail = e.error instanceof Error ? `: ${e.error.stack ?? e.error.message}` : "";
                onError(laneId, new Error(`${e.message || "worker failed to load"}${where}${detail}`));
            };
            w.onmessageerror = () => onError(laneId, new Error("worker message deserialization failed"));
            this.#workers.push(w);
            this.lanes.push({ post: (msg, transfer) => w.postMessage(msg, transfer ?? []) });
        }
    }
    dispose() {
        for (const w of this.#workers)
            w.terminate();
        this.#workers = [];
    }
}
/** Lanes backed by in-process {@link WorkerCore}s (single thread). */
export class InlineBackend {
    laneCount;
    lanes = [];
    #cores = [];
    constructor(laneCount, onMessage) {
        this.laneCount = laneCount;
        for (let i = 0; i < laneCount; i++) {
            const laneId = i;
            const core = new WorkerCore((msg) => deliver(() => onMessage(laneId, msg)));
            this.#cores.push(core);
            // The transfer list is irrelevant inline (no structured clone); ignore it.
            this.lanes.push({ post: (msg) => deliver(() => core.handle(msg)) });
        }
    }
    dispose() {
        for (const lane of this.lanes)
            lane.post({ kind: "dispose" });
    }
}
//# sourceMappingURL=backend.js.map