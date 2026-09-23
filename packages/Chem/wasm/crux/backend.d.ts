import type { MainToWorker, WorkerToMain } from "./protocol.js";
/** A single lane the pool can post to. */
export interface Lane {
    post(msg: MainToWorker, transfer?: Transferable[]): void;
}
export interface Backend {
    readonly laneCount: number;
    readonly lanes: Lane[];
    dispose(): void;
}
/** Lanes backed by real Web Workers (browser). */
export declare class WorkerBackend implements Backend {
    #private;
    readonly laneCount: number;
    readonly lanes: Lane[];
    constructor(laneCount: number, workerUrl: URL | string | undefined, onMessage: (laneId: number, msg: WorkerToMain) => void, onError: (laneId: number, err: Error) => void);
    dispose(): void;
}
/** Lanes backed by in-process {@link WorkerCore}s (single thread). */
export declare class InlineBackend implements Backend {
    #private;
    readonly laneCount: number;
    readonly lanes: Lane[];
    constructor(laneCount: number, onMessage: (laneId: number, msg: WorkerToMain) => void);
    dispose(): void;
}
//# sourceMappingURL=backend.d.ts.map