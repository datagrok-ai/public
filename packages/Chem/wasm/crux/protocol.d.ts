/**
 * The Web Worker message protocol for {@link CruxPool}. Raw `postMessage` over a
 * typed discriminated union (no Comlink — keeps the package dependency-free and
 * gives explicit control of the transferable list). All types are
 * `import type`-only; this module emits no runtime code.
 *
 * One worker (a "lane") owns several shards — each shard is a single WASM
 * `Collection` built from a contiguous slice of the dataset. Hit indices leaving
 * a worker are already rewritten to GLOBAL input positions (`shardOffset` added
 * inside the worker), so the main thread only merges.
 */
/** Monotonic id for a search; a `cancel` with this id supersedes/aborts it. */
export type QueryId = number;
/** Stable id of a shard (one WASM `Collection` in one lane). */
export type ShardId = number;
export type MainToWorker = {
    kind: "init";
    /** Precompiled module shared by every lane (structured-cloneable). */
    module: WebAssembly.Module;
    laneId: number;
} | {
    kind: "loadShard";
    shardId: ShardId;
    /** Global input position of this shard's first molecule. */
    shardOffset: number;
    /** Newline-joined SMILES, UTF-8 encoded; transferred (zero-copy). */
    payload: ArrayBuffer;
    buildIndex: boolean;
    chunkSize: number;
} | {
    kind: "searchSub";
    queryId: QueryId;
    smarts: string;
    limit: number;
} | {
    kind: "searchSim";
    queryId: QueryId;
    smiles: string;
    threshold: number;
    limit: number;
} | {
    kind: "cancel";
    queryId: QueryId;
} | {
    kind: "dropAll";
} | {
    kind: "dispose";
};
export type WorkerToMain = {
    kind: "ready";
    laneId: number;
} | {
    kind: "shardProgress";
    shardId: ShardId;
    done: number;
    total: number;
} | {
    kind: "shardLoaded";
    shardId: ShardId;
    /** Molecules successfully parsed in this shard. */
    size: number;
    /** Molecules handed to this shard (size + parse failures). */
    total: number;
    indexed: boolean;
} | {
    kind: "subPartial";
    queryId: QueryId;
    shardId: ShardId;
    /** Matching molecules' GLOBAL input positions. */
    indices: Uint32Array;
} | {
    kind: "simPartial";
    queryId: QueryId;
    shardId: ShardId;
    /** Hits' GLOBAL input positions, parallel to `scores`, sorted desc. */
    indices: Uint32Array;
    scores: Float32Array;
} | {
    kind: "searchDone";
    queryId: QueryId;
    cancelled: boolean;
} | {
    kind: "error";
    queryId?: QueryId;
    shardId?: ShardId;
    message: string;
} | {
    kind: "disposed";
};
//# sourceMappingURL=protocol.d.ts.map