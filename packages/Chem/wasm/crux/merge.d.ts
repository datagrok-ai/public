/**
 * Pure merge helpers that combine the per-shard partial results (already
 * carrying GLOBAL input positions) into a single globally-correct result.
 * Side-effect-free and typed-array-based, so they're cheap and unit-testable.
 */
import type { SimilarityResult } from "./collection.js";
/** Concatenate typed-array segments into one `Uint32Array`. */
export declare function concatUint32(parts: Uint32Array[]): Uint32Array;
/**
 * Merge per-shard substructure hits. Shards are disjoint, so there are no
 * duplicates: concatenate, sort by ascending index (stable display order), and
 * cap at `limit` (0 = unlimited).
 */
export declare function mergeSubstructure(parts: Uint32Array[], limit: number): Uint32Array;
/**
 * Merge per-shard similarity hits into the global top-`limit`. A global top-N
 * member is necessarily in its own shard's top-N (scores are query-relative and
 * shard-independent), so concatenating the shard tops, sorting by descending
 * score, and slicing to `limit` is exactly the single-collection result.
 */
export declare function mergeSimilarity(parts: SimilarityResult[], limit: number): SimilarityResult;
//# sourceMappingURL=merge.d.ts.map