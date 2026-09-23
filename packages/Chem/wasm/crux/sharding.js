/**
 * Shard-planning policy: how to split a dataset of `total` molecules across
 * `workers` lanes. Pure and side-effect-free, so it can be unit-tested without
 * spinning up any workers.
 *
 * Shards are JS-managed (one WASM `Collection` each), NOT crux-core's internal
 * multipart shards. We aim for a few shards per worker so the worker has a queue
 * it can stream from and cancel between, while keeping each shard large enough
 * that an indexed shard's rarity table stays statistically useful.
 */
/** Target molecules per shard when nothing is overridden. */
export const TARGET_SHARD_MOLECULES = 250_000;
/** Cap on shards per worker (bounds per-shard fixed overhead). */
export const MAX_SHARDS_PER_WORKER = 4;
/** Don't slice below this many molecules per shard (keeps shards worthwhile). */
export const MIN_SHARD_MOLECULES = 10_000;
/** Default shard count given `total` molecules and `workers` lanes. */
export function defaultShardCount(total, workers) {
    const w = Math.max(1, Math.floor(workers));
    if (total <= MIN_SHARD_MOLECULES)
        return 1;
    // Never make a shard smaller than the floor...
    const maxBySize = Math.max(1, Math.floor(total / MIN_SHARD_MOLECULES));
    // ...aim for ~one shard per TARGET molecules, but use at least every worker
    // and at most MAX_SHARDS_PER_WORKER per worker.
    const byTarget = Math.max(1, Math.ceil(total / TARGET_SHARD_MOLECULES));
    const desired = Math.min(Math.max(byTarget, w), w * MAX_SHARDS_PER_WORKER);
    return Math.max(1, Math.min(desired, maxBySize));
}
/**
 * Plan the shards for `total` molecules across `workers` lanes. `opts.shards`
 * (exact count) wins over `opts.shardSize` (target molecules per shard), which
 * wins over the default policy. Shards are contiguous, near-equal input ranges.
 */
export function planShards(total, workers, opts = {}) {
    let count;
    if (opts.shards && opts.shards > 0)
        count = Math.floor(opts.shards);
    else if (opts.shardSize && opts.shardSize > 0)
        count = Math.max(1, Math.ceil(total / opts.shardSize));
    else
        count = defaultShardCount(total, workers);
    // Can't have more shards than molecules (and always at least one shard).
    count = Math.max(1, Math.min(count, Math.max(1, total)));
    const size = Math.ceil(total / count);
    const shards = [];
    for (let start = 0; start < total; start += size) {
        shards.push({ start, end: Math.min(start + size, total) });
    }
    // An empty dataset still yields a single empty shard so progress/searches
    // have something to complete against.
    if (shards.length === 0)
        shards.push({ start: 0, end: 0 });
    return { shardCount: shards.length, shards };
}
//# sourceMappingURL=sharding.js.map