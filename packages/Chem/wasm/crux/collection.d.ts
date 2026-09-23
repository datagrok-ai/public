/**
 * A similarity-search result as two parallel arrays, ordered by descending
 * score: `indices[i]` is a hit's original input position, `scores[i]` its
 * Tanimoto similarity in `[0, 1]`.
 */
export interface SimilarityResult {
    indices: Uint32Array;
    scores: Float32Array;
}
/** Options for {@link CruxCollection.fromSmiles}. */
export interface BuildOptions {
    /** Molecules processed per WASM call between progress ticks. Default 5000. */
    chunkSize?: number;
    /** Called after each chunk with the running counts. */
    onProgress?: (done: number, total: number) => void;
    /**
     * Build a fingerprint index (default `false`). When `true`: fast queries,
     * more memory and build time. When `false` (the default): store molecules
     * only and search by direct graph matching / on-the-fly fingerprinting — no
     * build cost, slower queries. Repeated querying of one dataset favours `true`.
     */
    buildIndex?: boolean;
}
/** Options for {@link CruxCollection.substructureSearch}. */
export interface SubstructureOptions {
    /** Max hits to return (0 = unlimited). Default 0. */
    limit?: number;
}
/** Options for {@link CruxCollection.similaritySearch}. */
export interface SimilarityOptions {
    /** Minimum Tanimoto similarity in [0, 1]. Default 0. */
    threshold?: number;
    /** If > 0, return only the top-N hits by score (with `threshold` as a floor). Default 0. */
    limit?: number;
}
/**
 * An in-memory, searchable collection of molecules.
 *
 * Build it with {@link CruxCollection.fromSmiles}; the heavy work (parsing and
 * fingerprinting every molecule) happens there. Each search is then a fast scan.
 */
export declare class CruxCollection {
    #private;
    private constructor();
    /**
     * Build a collection from a list of SMILES. Unparseable / unsupported SMILES
     * are skipped (the resulting `size` may be smaller than `smiles.length`), but
     * hit indices still refer to positions in the original `smiles` array.
     * Work is chunked so a host (e.g. a Web Worker) can report progress.
     */
    static fromSmiles(smiles: string[], opts?: BuildOptions): Promise<CruxCollection>;
    /** Number of molecules in the collection. */
    get size(): number;
    /** Whether this collection has a fingerprint index (vs. direct matching). */
    get indexed(): boolean;
    /**
     * Substructure (SMARTS) search. Returns the matching molecules' original
     * input positions, capped at `opts.limit` (0 = unlimited).
     */
    substructureSearch(smarts: string, opts?: SubstructureOptions): Uint32Array;
    /**
     * Similarity (ECFP4 Tanimoto) search against a query SMILES. Returns parallel
     * `indices` / `scores` arrays sorted by descending score.
     */
    similaritySearch(querySmiles: string, opts?: SimilarityOptions): SimilarityResult;
    /** Release the underlying WASM memory. The instance is unusable afterwards. */
    free(): void;
}
//# sourceMappingURL=collection.d.ts.map