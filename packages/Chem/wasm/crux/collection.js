/**
 * Ergonomic TypeScript wrapper around a single WASM `Collection` — load a set
 * of molecules once, then run substructure / similarity search against it.
 *
 * This is the low-level **single-index engine**: it parses + searches on the
 * calling thread. Most callers want {@link CruxPool} (in `./pool.ts`), which
 * shards a dataset across Web Workers and reuses this engine inside each one.
 *
 * Results are returned as typed arrays of the molecules' **original input
 * positions** (the index into the `smiles` array you passed to
 * {@link CruxCollection.fromSmiles}); map those back to your own data. SMILES
 * are not echoed back — you already have them.
 */
import { initCrux, CollectionBuilder, Collection as WasmCollection } from "./wasm.js";
/**
 * An in-memory, searchable collection of molecules.
 *
 * Build it with {@link CruxCollection.fromSmiles}; the heavy work (parsing and
 * fingerprinting every molecule) happens there. Each search is then a fast scan.
 */
export class CruxCollection {
    #inner;
    constructor(inner) {
        this.#inner = inner;
    }
    /**
     * Build a collection from a list of SMILES. Unparseable / unsupported SMILES
     * are skipped (the resulting `size` may be smaller than `smiles.length`), but
     * hit indices still refer to positions in the original `smiles` array.
     * Work is chunked so a host (e.g. a Web Worker) can report progress.
     */
    static async fromSmiles(smiles, opts = {}) {
        await initCrux();
        const chunkSize = opts.chunkSize ?? 5000;
        const total = smiles.length;
        const builder = new CollectionBuilder(opts.buildIndex ?? false);
        try {
            for (let i = 0; i < total; i += chunkSize) {
                builder.addMany(smiles.slice(i, i + chunkSize));
                opts.onProgress?.(Math.min(i + chunkSize, total), total);
                // Yield to the event loop so progress messages flush and the host
                // stays responsive between chunks.
                await new Promise((resolve) => setTimeout(resolve, 0));
            }
            return new CruxCollection(builder.finish());
        }
        catch (err) {
            builder.free();
            throw err;
        }
    }
    /** Number of molecules in the collection. */
    get size() {
        return this.#inner.size();
    }
    /** Whether this collection has a fingerprint index (vs. direct matching). */
    get indexed() {
        return this.#inner.isIndexed();
    }
    /**
     * Substructure (SMARTS) search. Returns the matching molecules' original
     * input positions, capped at `opts.limit` (0 = unlimited).
     */
    substructureSearch(smarts, opts = {}) {
        return this.#inner.substructureSearch(smarts, opts.limit ?? 0);
    }
    /**
     * Similarity (ECFP4 Tanimoto) search against a query SMILES. Returns parallel
     * `indices` / `scores` arrays sorted by descending score.
     */
    similaritySearch(querySmiles, opts = {}) {
        const res = this.#inner.similaritySearch(querySmiles, opts.threshold ?? 0, opts.limit ?? 0);
        // The getters copy out of WASM memory, so the arrays outlive `res`.
        const out = { indices: res.indices, scores: res.scores };
        res.free();
        return out;
    }
    /** Release the underlying WASM memory. The instance is unusable afterwards. */
    free() {
        this.#inner.free();
    }
}
//# sourceMappingURL=collection.js.map