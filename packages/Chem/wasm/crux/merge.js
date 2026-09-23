/** Concatenate typed-array segments into one `Uint32Array`. */
export function concatUint32(parts) {
    let n = 0;
    for (const p of parts)
        n += p.length;
    const out = new Uint32Array(n);
    let off = 0;
    for (const p of parts) {
        out.set(p, off);
        off += p.length;
    }
    return out;
}
/**
 * Merge per-shard substructure hits. Shards are disjoint, so there are no
 * duplicates: concatenate, sort by ascending index (stable display order), and
 * cap at `limit` (0 = unlimited).
 */
export function mergeSubstructure(parts, limit) {
    const all = concatUint32(parts);
    all.sort(); // typed arrays sort numerically
    return limit > 0 && all.length > limit ? all.slice(0, limit) : all;
}
/**
 * Merge per-shard similarity hits into the global top-`limit`. A global top-N
 * member is necessarily in its own shard's top-N (scores are query-relative and
 * shard-independent), so concatenating the shard tops, sorting by descending
 * score, and slicing to `limit` is exactly the single-collection result.
 */
export function mergeSimilarity(parts, limit) {
    const pairs = [];
    for (const p of parts) {
        for (let k = 0; k < p.indices.length; k++)
            pairs.push({ i: p.indices[k], s: p.scores[k] });
    }
    pairs.sort((a, b) => b.s - a.s);
    const n = limit > 0 ? Math.min(limit, pairs.length) : pairs.length;
    const indices = new Uint32Array(n);
    const scores = new Float32Array(n);
    for (let k = 0; k < n; k++) {
        indices[k] = pairs[k].i;
        scores[k] = pairs[k].s;
    }
    return { indices, scores };
}
//# sourceMappingURL=merge.js.map