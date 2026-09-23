/**
 * crux-js — browser/WASM consumer of the Crux chemical-search engine.
 *
 * Wraps the `crux-wasm` artifact built from `../crux-core` and exposes
 * substructure / similarity search over an in-memory collection of molecules.
 *
 * For most apps, use {@link CruxPool}: it shards a dataset across all CPU cores
 * (Web Workers), parses in parallel, and broadcasts each query to every shard —
 * with progress, cancellation, streaming and a search-as-you-draw session.
 *
 * ```ts
 * import { CruxPool } from "@datagrok/crux-js";
 *
 * await using pool = await CruxPool.create();
 * await using ds = await pool.load(smiles, { onProgress: (p) => bar(p.done / p.total) });
 * const live = ds.session({ debounceMs: 50 });
 * for await (const b of live.searchStream({ kind: "substructure", smarts }))
 *   render(b.indices); // GLOBAL input positions; map to your data
 * ```
 *
 * {@link CruxCollection} is the low-level single-thread engine (one in-memory
 * index, no workers); the pool reuses it inside each worker.
 */
export const VERSION = "0.2.0";
export { initCrux } from "./wasm.js";
export { CruxCollection, } from "./collection.js";
export { CruxPool, CruxDataset, SearchSession, } from "./pool.js";
//# sourceMappingURL=index.js.map