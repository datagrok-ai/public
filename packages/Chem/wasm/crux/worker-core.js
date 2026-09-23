/**
 * The per-lane engine, independent of how messages get to it. The real Web
 * Worker (`pool.worker.ts`) wires this to `self`; the single-thread
 * {@link InlineBackend} wires it to an in-process channel. Either way it owns a
 * lane's shards (one WASM `Collection` each) and runs searches across them,
 * streaming a partial per shard and checking for cancellation between shards.
 */
import { initCrux } from "./wasm.js";
import { CruxCollection } from "./collection.js";
const errMessage = (e) => (e instanceof Error ? e.message : String(e));
/**
 * A macrotask yield. Crucially NOT a microtask: incoming `postMessage`s (e.g. a
 * `cancel`) are delivered as macrotasks, so only a macrotask checkpoint between
 * shards lets a queued cancellation be observed mid-search. `setTimeout(0)` is
 * used for portability (Node / worker / main thread); its small clamp only adds
 * a few ms between shards, after the first shard's results have already streamed.
 */
function yieldToEvents() {
    return new Promise((resolve) => setTimeout(resolve, 0));
}
export class WorkerCore {
    #post;
    #shards = new Map();
    /** Build requests run one-at-a-time per lane (FIFO), not interleaved. */
    #loadChain = Promise.resolve();
    /** Query ids explicitly cancelled/superseded; checked between shards. */
    #cancelled = new Set();
    constructor(post) {
        this.#post = post;
    }
    handle(msg) {
        switch (msg.kind) {
            case "init":
                this.#init(msg.module, msg.laneId);
                break;
            case "loadShard":
                this.#loadChain = this.#loadChain.then(() => this.#buildShard(msg));
                break;
            case "searchSub":
                void this.#runSub(msg.queryId, msg.smarts, msg.limit);
                break;
            case "searchSim":
                void this.#runSim(msg.queryId, msg.smiles, msg.threshold, msg.limit);
                break;
            case "cancel":
                this.#cancelled.add(msg.queryId);
                break;
            case "dropAll":
                this.#dropAll();
                break;
            case "dispose":
                this.#dropAll();
                this.#post({ kind: "disposed" });
                break;
        }
    }
    #init(module, laneId) {
        initCrux(module).then(() => this.#post({ kind: "ready", laneId }), (e) => this.#post({ kind: "error", message: errMessage(e) }));
    }
    async #buildShard(msg) {
        try {
            const text = new TextDecoder().decode(new Uint8Array(msg.payload));
            const smiles = text.length > 0 ? text.split("\n") : [];
            const col = await CruxCollection.fromSmiles(smiles, {
                buildIndex: msg.buildIndex,
                chunkSize: msg.chunkSize,
                onProgress: (done, total) => this.#post({ kind: "shardProgress", shardId: msg.shardId, done, total }),
            });
            this.#shards.set(msg.shardId, { offset: msg.shardOffset, col });
            this.#post({
                kind: "shardLoaded",
                shardId: msg.shardId,
                size: col.size,
                total: smiles.length,
                indexed: col.indexed,
            });
        }
        catch (e) {
            this.#post({ kind: "error", shardId: msg.shardId, message: errMessage(e) });
        }
    }
    async #runSub(queryId, smarts, limit) {
        try {
            for (const [shardId, shard] of this.#shards) {
                if (this.#cancelled.has(queryId)) {
                    this.#post({ kind: "searchDone", queryId, cancelled: true });
                    return;
                }
                const idx = shard.col.substructureSearch(smarts, { limit });
                for (let i = 0; i < idx.length; i++)
                    idx[i] += shard.offset;
                this.#post({ kind: "subPartial", queryId, shardId, indices: idx }, [idx.buffer]);
                await yieldToEvents();
            }
            this.#post({ kind: "searchDone", queryId, cancelled: this.#cancelled.has(queryId) });
        }
        catch (e) {
            this.#post({ kind: "error", queryId, message: errMessage(e) });
            this.#post({ kind: "searchDone", queryId, cancelled: true });
        }
        finally {
            this.#cancelled.delete(queryId);
        }
    }
    async #runSim(queryId, smiles, threshold, limit) {
        try {
            for (const [shardId, shard] of this.#shards) {
                if (this.#cancelled.has(queryId)) {
                    this.#post({ kind: "searchDone", queryId, cancelled: true });
                    return;
                }
                const { indices, scores } = shard.col.similaritySearch(smiles, { threshold, limit });
                for (let i = 0; i < indices.length; i++)
                    indices[i] += shard.offset;
                this.#post({ kind: "simPartial", queryId, shardId, indices, scores }, [
                    indices.buffer,
                    scores.buffer,
                ]);
                await yieldToEvents();
            }
            this.#post({ kind: "searchDone", queryId, cancelled: this.#cancelled.has(queryId) });
        }
        catch (e) {
            this.#post({ kind: "error", queryId, message: errMessage(e) });
            this.#post({ kind: "searchDone", queryId, cancelled: true });
        }
        finally {
            this.#cancelled.delete(queryId);
        }
    }
    #dropAll() {
        for (const shard of this.#shards.values())
            shard.col.free();
        this.#shards.clear();
    }
}
//# sourceMappingURL=worker-core.js.map