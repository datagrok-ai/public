/// <reference lib="webworker" />
/**
 * Web Worker entry for {@link CruxPool}. Thin: it just wires a {@link WorkerCore}
 * to the worker's `postMessage` / `onmessage`. All real logic lives in
 * `worker-core.ts` so it is shared with the single-thread `InlineBackend`.
 */
import { WorkerCore } from "./worker-core.js";
const post = (msg) => self.postMessage(msg);
const core = new WorkerCore((msg, transfer) => self.postMessage(msg, transfer ?? []));
self.onmessage = (e) => core.handle(e.data);
// Surface otherwise-silent runtime / unhandled-rejection errors to the pool.
self.addEventListener("error", (e) => post({ kind: "error", message: e.message || String(e.error ?? "worker error") }));
self.addEventListener("unhandledrejection", (e) => post({ kind: "error", message: String(e.reason ?? "worker rejection") }));
export {};
//# sourceMappingURL=pool.worker.js.map