/* tslint:disable */
/* eslint-disable */

export class WasmBitBirch {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Fit fingerprints, resetting any existing tree.
     *
     * `fingerprints` is a packed byte array: n_fps consecutive fingerprints,
     * each ceil(n_features/8) bytes long.
     */
    fit(fingerprints: Uint8Array, n_fps: number): void;
    /**
     * Get cluster assignment for each fitted molecule.
     * Returns an array where result[i] = cluster_id for molecule i.
     */
    get_assignments(): Int32Array;
    /**
     * Get all leaf centroids as a flat packed byte array.
     */
    get_centroids(): Uint8Array;
    /**
     * Get molecule indices grouped by cluster, as JSON string.
     * Format: [[0,5,12],[3,7],...]
     */
    get_cluster_mol_ids_json(): string;
    /**
     * Number of clusters (leaf subclusters).
     */
    n_clusters(): number;
    /**
     * Number of fingerprints fitted so far.
     */
    n_fitted(): number;
    /**
     * Create a new BitBIRCH clusterer.
     *
     * - `threshold`: Tanimoto similarity threshold for merging (0.0–1.0)
     * - `branching_factor`: max subclusters per node before splitting
     * - `n_features`: number of bits in each fingerprint (e.g., 2048 for ECFP4)
     */
    constructor(threshold: number, branching_factor: number, n_features: number);
    /**
     * Incrementally add fingerprints to the existing tree.
     */
    partial_fit(fingerprints: Uint8Array, n_fps: number): void;
    /**
     * Predict cluster assignments for new fingerprints (read-only).
     */
    predict(fingerprints: Uint8Array, n_fps: number): Int32Array;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_wasmbitbirch_free: (a: number, b: number) => void;
    readonly wasmbitbirch_fit: (a: number, b: number, c: number, d: number) => [number, number];
    readonly wasmbitbirch_get_assignments: (a: number) => [number, number];
    readonly wasmbitbirch_get_centroids: (a: number) => [number, number];
    readonly wasmbitbirch_get_cluster_mol_ids_json: (a: number) => [number, number];
    readonly wasmbitbirch_n_clusters: (a: number) => number;
    readonly wasmbitbirch_n_fitted: (a: number) => number;
    readonly wasmbitbirch_new: (a: number, b: number, c: number) => [number, number, number];
    readonly wasmbitbirch_partial_fit: (a: number, b: number, c: number, d: number) => [number, number];
    readonly wasmbitbirch_predict: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
