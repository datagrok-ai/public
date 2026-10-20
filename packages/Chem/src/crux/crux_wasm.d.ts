/* tslint:disable */
/* eslint-disable */

/**
 * A finalised, queryable collection of molecules.
 *
 * `indexed` collections carry fingerprint sidecars and answer queries via the
 * screened matcher / sidecar Tanimoto kernels. Non-indexed ("direct")
 * collections store molecules only and answer queries by brute-force graph
 * matching / on-the-fly fingerprinting over every molecule.
 */
export class Collection {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Build an indexed collection from a single batch of SMILES (convenience
     * for callers that don't need progress reporting).
     */
    static fromSmiles(smiles: string[]): Collection;
    /**
     * Whether this collection has a fingerprint index (vs. direct matching).
     */
    isIndexed(): boolean;
    /**
     * Similarity (ECFP4 Tanimoto) search. If `top_k > 0`, returns the top
     * `top_k` hits with Tanimoto >= `threshold` (a floor); otherwise returns
     * every hit with Tanimoto >= `threshold`. Results are sorted by score
     * descending — a {@link SimilarityResult} of parallel `indices` / `scores`.
     */
    similaritySearch(query_smiles: string, threshold: number, top_k: number): SimilarityResult;
    /**
     * Number of molecules in the collection.
     */
    size(): number;
    /**
     * Substructure (SMARTS) search. Returns a `Uint32Array` of the matching
     * molecules' original input positions, capped at `limit` (0 = unlimited).
     */
    substructureSearch(smarts: string, limit: number): Uint32Array;
}

/**
 * Incrementally builds a [`Collection`] from SMILES strings.
 *
 * When `build_index` is true the builder also computes the screening + ECFP4
 * fingerprints and assembles the sidecars (fast queries, larger memory). When
 * false it stores molecules only — searches then run by direct graph matching
 * / on-the-fly fingerprinting over every molecule (no build cost, slower
 * queries).
 */
export class CollectionBuilder {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Parse + fingerprint a chunk of SMILES, appending each to the index.
     * Returns the number successfully added; unparseable / unsupported SMILES
     * are skipped and counted in [`CollectionBuilder::failed`].
     */
    addMany(smiles: string[]): number;
    /**
     * Number of molecules accepted so far.
     */
    added(): number;
    /**
     * Number of SMILES skipped because they failed to parse / encode.
     */
    failed(): number;
    /**
     * Finalise into a queryable [`Collection`]. Builds the sidecars + rarity
     * table when indexing; otherwise produces a molecules-only collection.
     */
    finish(): Collection;
    constructor(build_index: boolean);
    /**
     * Read a SMILES RDKit's default read rejects the way RDKit reads it with Kekulize left out
     * of its sanitization (MinimalLib's `get_mol(smiles, {"kekulize": false})`, the retry
     * Datagrok's Chem makes): aromatic rings written without their `[nH]`, aromatic
     * phosphazenes. Off by default.
     */
    setLenient(lenient: boolean): void;
}

/**
 * A similarity-search result as two parallel typed arrays: `indices[i]` is the
 * caller's original input position of the i-th hit, `scores[i]` its Tanimoto
 * similarity in `[0, 1]`. Hits are ordered by descending score.
 */
export class SimilarityResult {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Original input positions of the hits (a `Uint32Array` on the JS side).
     */
    readonly indices: Uint32Array;
    /**
     * Tanimoto scores, parallel to `indices` (a `Float32Array` on the JS side).
     */
    readonly scores: Float32Array;
}

/**
 * An in-browser synthon-space collection built from a reaction/synthon CSV
 * (RDKit text format; connectors `[U]`/`[Np]` or `[n*]`). Searches enumerate +
 * verify the combinatorial product space without materialising it, returning
 * product SMILES + provenance. The synthon corpora are tiny (~256 KB) relative
 * to their product space, so the whole space lives in memory; the per-query
 * `SearchEngine` / `SimilarityEngine` (which build a CXMOL / FP index over the
 * synthon pool) are constructed per call — cheap for these pool sizes.
 */
export class SynthonCollection {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Build a collection from the text of a synthon CSV. Builds the substructure
     * index once here so each search is a cheap view, not a rebuild.
     */
    static fromCsv(csv_text: string): SynthonCollection;
    /**
     * Nominal product-space size (sum of per-reaction product upper bounds).
     */
    numProducts(): bigint;
    /**
     * Number of reactions.
     */
    numReactions(): number;
    /**
     * Number of distinct synthons in the pool.
     */
    numSynthons(): number;
    /**
     * Similarity search (Crux ECFP4 over assembled products). Returns the highest
     * `top_k` products with Tanimoto ≥ `cutoff` (`top_k = 0` = all ≥ cutoff),
     * sorted descending.
     */
    similaritySearch(query_smiles: string, cutoff: number, top_k: number): SynthonSimHits;
    /**
     * Substructure search: the query is parsed as SMILES (element + aromaticity +
     * bond order). `limit` caps the listing (0 = unlimited).
     */
    substructureSearch(query_smiles: string, limit: number): SynthonHits;
}

/**
 * A synthon-space substructure result: parallel `names` / `smiles` arrays, one
 * entry per hit. `names[i]` is the RDKit-compatible product name
 * `synthonId0;…;reactionId`; `smiles[i]` the assembled product SMILES.
 */
export class SynthonHits {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    len(): number;
    readonly names: string[];
    readonly smiles: string[];
}

/**
 * A synthon-space similarity result: parallel `names` / `smiles` / `scores`
 * arrays, ordered by descending Tanimoto (`scores[i]` in `[0, 1]`).
 */
export class SynthonSimHits {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    len(): number;
    readonly names: string[];
    readonly scores: Float32Array;
    readonly smiles: string[];
}

/**
 * Install a panic hook that forwards Rust panics to `console.error`. Called
 * once at module init on wasm; a no-op on native targets.
 */
export function start(): void;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_collection_free: (a: number, b: number) => void;
    readonly __wbg_collectionbuilder_free: (a: number, b: number) => void;
    readonly __wbg_similarityresult_free: (a: number, b: number) => void;
    readonly __wbg_synthoncollection_free: (a: number, b: number) => void;
    readonly __wbg_synthonhits_free: (a: number, b: number) => void;
    readonly __wbg_synthonsimhits_free: (a: number, b: number) => void;
    readonly collection_fromSmiles: (a: number, b: number) => number;
    readonly collection_isIndexed: (a: number) => number;
    readonly collection_similaritySearch: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly collection_size: (a: number) => number;
    readonly collection_substructureSearch: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly collectionbuilder_addMany: (a: number, b: number, c: number) => number;
    readonly collectionbuilder_added: (a: number) => number;
    readonly collectionbuilder_failed: (a: number) => number;
    readonly collectionbuilder_finish: (a: number) => number;
    readonly collectionbuilder_new: (a: number) => number;
    readonly collectionbuilder_setLenient: (a: number, b: number) => void;
    readonly similarityresult_indices: (a: number) => [number, number];
    readonly similarityresult_scores: (a: number) => [number, number];
    readonly synthoncollection_fromCsv: (a: number, b: number) => [number, number, number];
    readonly synthoncollection_numProducts: (a: number) => bigint;
    readonly synthoncollection_numReactions: (a: number) => number;
    readonly synthoncollection_numSynthons: (a: number) => number;
    readonly synthoncollection_similaritySearch: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly synthoncollection_substructureSearch: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly synthonhits_len: (a: number) => number;
    readonly synthonhits_names: (a: number) => [number, number];
    readonly synthonhits_smiles: (a: number) => [number, number];
    readonly synthonsimhits_names: (a: number) => [number, number];
    readonly synthonsimhits_scores: (a: number) => [number, number];
    readonly synthonsimhits_smiles: (a: number) => [number, number];
    readonly start: () => void;
    readonly synthonsimhits_len: (a: number) => number;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __externref_table_alloc: () => number;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __externref_drop_slice: (a: number, b: number) => void;
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
