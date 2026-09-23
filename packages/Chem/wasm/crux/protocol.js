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
export {};
//# sourceMappingURL=protocol.js.map