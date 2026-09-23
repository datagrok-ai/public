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
export declare const TARGET_SHARD_MOLECULES = 250000;
/** Cap on shards per worker (bounds per-shard fixed overhead). */
export declare const MAX_SHARDS_PER_WORKER = 4;
/** Don't slice below this many molecules per shard (keeps shards worthwhile). */
export declare const MIN_SHARD_MOLECULES = 10000;
/** A contiguous half-open input range `[start, end)` assigned to one shard. */
export interface ShardRange {
    start: number;
    end: number;
}
export interface ShardPlan {
    shardCount: number;
    shards: ShardRange[];
}
/** Options that override the default shard count / size. */
export interface ShardOptions {
    shards?: number;
    shardSize?: number;
}
/** Default shard count given `total` molecules and `workers` lanes. */
export declare function defaultShardCount(total: number, workers: number): number;
/**
 * Plan the shards for `total` molecules across `workers` lanes. `opts.shards`
 * (exact count) wins over `opts.shardSize` (target molecules per shard), which
 * wins over the default policy. Shards are contiguous, near-equal input ranges.
 */
export declare function planShards(total: number, workers: number, opts?: ShardOptions): ShardPlan;
//# sourceMappingURL=sharding.d.ts.map