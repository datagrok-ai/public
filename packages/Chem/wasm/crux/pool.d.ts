import type { SimilarityResult } from "./collection.js";
export type { SimilarityResult } from "./collection.js";
/** Options for {@link CruxPool.create}. */
export interface PoolOptions {
    /** Worker count. Default `navigator.hardwareConcurrency` (clamped ≥ 1). */
    workers?: number;
    /**
     * How the `.wasm` is obtained once on the main thread and shared to workers
     * (mirrors {@link initCrux}). Omit in a bundler/browser to resolve the sibling
     * `.wasm`. In Node, pass the bytes (no fetch).
     */
    wasm?: BufferSource | Response | WebAssembly.Module | URL | string;
    /** Override the worker entry URL (advanced). */
    workerUrl?: URL | string;
}
/** Progress snapshot during {@link CruxPool.load}. */
export interface LoadProgress {
    /** Molecules parsed across all shards. */
    done: number;
    /** Total molecules (== input length). */
    total: number;
    /** Shards finished building. */
    shardsDone: number;
    shardCount: number;
}
/** Options for {@link CruxPool.load}. */
export interface LoadOptions {
    /** Build a fingerprint index per shard (default `false`; see {@link BuildOptions}). */
    buildIndex?: boolean;
    /** Exact shard count (overrides the default policy and `shardSize`). */
    shards?: number;
    /** Target molecules per shard (overrides the default policy). */
    shardSize?: number;
    /** Molecules per WASM `addMany` call inside a worker. Default 5000. */
    chunkSize?: number;
    /** Aggregated progress across all shards. */
    onProgress?: (p: LoadProgress) => void;
    /** Abort the load (frees any shards that already built). */
    signal?: AbortSignal;
}
/** Per-shard progress of a running search. */
export interface SearchProgress {
    shardsDone: number;
    shardCount: number;
    /** Hits accumulated so far. */
    hits: number;
}
/** Options for a substructure search. */
export interface SubstructureSearchOptions {
    /** Global max hits (0 = unlimited). Default 0. */
    limit?: number;
    signal?: AbortSignal;
    onProgress?: (p: SearchProgress) => void;
}
/** Options for a similarity search. */
export interface SimilaritySearchOptions {
    /** Minimum Tanimoto similarity in [0, 1]. Default 0. */
    threshold?: number;
    /** Global top-N by score (0 = all ≥ threshold). Default 0. */
    limit?: number;
    signal?: AbortSignal;
    onProgress?: (p: SearchProgress) => void;
}
/** One streamed batch of substructure hits, from a single shard. */
export interface SubstructureBatch {
    /** Matching molecules' GLOBAL input positions (sorted ascending). */
    indices: Uint32Array;
    shardsDone: number;
    shardCount: number;
}
/** One streamed batch of similarity hits, from a single shard. */
export interface SimilarityBatch {
    /** Hits' GLOBAL input positions, parallel to `scores`, sorted desc. */
    indices: Uint32Array;
    scores: Float32Array;
    shardsDone: number;
    shardCount: number;
}
/** Options for {@link CruxDataset.session}. */
export interface SearchSessionOptions {
    /** Debounce (ms) before a queued query dispatches, to coalesce draw events. Default 0. */
    debounceMs?: number;
}
/** A query for a {@link SearchSession}. */
export type SessionQuery = {
    kind: "substructure";
    smarts: string;
    options?: Omit<SubstructureSearchOptions, "signal">;
} | {
    kind: "similarity";
    smiles: string;
    options?: Omit<SimilaritySearchOptions, "signal">;
};
interface DatasetMeta {
    size: number;
    failed: number;
    indexed: boolean;
    shardCount: number;
}
export declare class CruxPool {
    #private;
    private constructor();
    /** Number of live worker lanes. */
    get workerCount(): number;
    /**
     * Spawn the workers, compile the WASM once on the main thread and share the
     * compiled module to every worker. Resolves when all workers are ready. In a
     * Node / no-Worker environment, falls back to a single-thread inline backend.
     */
    static create(opts?: PoolOptions): Promise<CruxPool>;
    /** Shard `smiles` across the workers and parse in parallel. */
    load(smiles: string[], opts?: LoadOptions): Promise<CruxDataset>;
    /** Terminate all workers and free their WASM memory. Idempotent. */
    dispose(): Promise<void>;
    [Symbol.asyncDispose](): Promise<void>;
    /** @internal */
    _startSubstructure(smarts: string, opts: SubstructureSearchOptions): AsyncIterable<SubstructureBatch>;
    /** @internal */
    _startSimilarity(smiles: string, opts: SimilaritySearchOptions): AsyncIterable<SimilarityBatch>;
    /** @internal Free the current dataset's shards and cancel its searches. */
    _dropDataset(): Promise<void>;
}
export declare class CruxDataset {
    #private;
    /** Molecules successfully loaded (parse failures excluded). */
    readonly size: number;
    /** Input SMILES that failed to parse. */
    readonly failed: number;
    readonly indexed: boolean;
    readonly shardCount: number;
    /** @internal Created by {@link CruxPool.load}. */
    constructor(pool: CruxPool, meta: DatasetMeta);
    /** Streaming substructure (SMARTS) search — one batch per shard as it finishes. */
    substructureSearchStream(smarts: string, opts?: SubstructureSearchOptions): AsyncIterable<SubstructureBatch>;
    /** Collected substructure search — merged, sorted by ascending index, capped at `limit`. */
    substructureSearch(smarts: string, opts?: SubstructureSearchOptions): Promise<Uint32Array>;
    /** Streaming similarity search — one batch per shard as it finishes. */
    similaritySearchStream(querySmiles: string, opts?: SimilaritySearchOptions): AsyncIterable<SimilarityBatch>;
    /** Collected similarity search — globally merged + re-sorted desc, capped at `limit`. */
    similaritySearch(querySmiles: string, opts?: SimilaritySearchOptions): Promise<SimilarityResult>;
    /** Create a latest-wins session for search-as-you-draw (auto-cancels the prior query). */
    session(opts?: SearchSessionOptions): SearchSession;
    /** Free this dataset's shard collections (the pool stays alive). Idempotent. */
    dispose(): Promise<void>;
    [Symbol.asyncDispose](): Promise<void>;
}
export declare class SearchSession {
    #private;
    /** @internal Created by {@link CruxDataset.session}. */
    constructor(dataset: CruxDataset, opts?: SearchSessionOptions);
    /** Run a query, auto-cancelling the previous one. Resolves to merged results. */
    search(query: SessionQuery): Promise<Uint32Array | SimilarityResult>;
    /** Stream a query, auto-cancelling the previous one. */
    searchStream(query: SessionQuery): AsyncIterable<SubstructureBatch | SimilarityBatch>;
    /** Cancel the current query without starting a new one. */
    cancel(): void;
    dispose(): void;
    [Symbol.dispose](): void;
}
//# sourceMappingURL=pool.d.ts.map