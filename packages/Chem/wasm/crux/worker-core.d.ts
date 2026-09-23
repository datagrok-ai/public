import type { MainToWorker, WorkerToMain } from "./protocol.js";
/** How a {@link WorkerCore} emits messages back to the main thread. */
export type PostFn = (msg: WorkerToMain, transfer?: Transferable[]) => void;
export declare class WorkerCore {
    #private;
    constructor(post: PostFn);
    handle(msg: MainToWorker): void;
}
//# sourceMappingURL=worker-core.d.ts.map