/* @ts-self-types="./crux_wasm.d.ts" */

/**
 * A finalised, queryable collection of molecules.
 *
 * `indexed` collections carry fingerprint sidecars and answer queries via the
 * screened matcher / sidecar Tanimoto kernels. Non-indexed ("direct")
 * collections store molecules only and answer queries by brute-force graph
 * matching / on-the-fly fingerprinting over every molecule.
 */
export class Collection {
    static __wrap(ptr) {
        const obj = Object.create(Collection.prototype);
        obj.__wbg_ptr = ptr;
        CollectionFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CollectionFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_collection_free(ptr, 0);
    }
    /**
     * Build an indexed collection from a single batch of SMILES (convenience
     * for callers that don't need progress reporting).
     * @param {string[]} smiles
     * @returns {Collection}
     */
    static fromSmiles(smiles) {
        const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collection_fromSmiles(ptr0, len0);
        return Collection.__wrap(ret);
    }
    /**
     * Whether this collection has a fingerprint index (vs. direct matching).
     * @returns {boolean}
     */
    isIndexed() {
        const ret = wasm.collection_isIndexed(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Similarity (ECFP4 Tanimoto) search. If `top_k > 0`, returns the top
     * `top_k` hits with Tanimoto >= `threshold` (a floor); otherwise returns
     * every hit with Tanimoto >= `threshold`. Results are sorted by score
     * descending — a {@link SimilarityResult} of parallel `indices` / `scores`.
     * @param {string} query_smiles
     * @param {number} threshold
     * @param {number} top_k
     * @returns {SimilarityResult}
     */
    similaritySearch(query_smiles, threshold, top_k) {
        const ptr0 = passStringToWasm0(query_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collection_similaritySearch(this.__wbg_ptr, ptr0, len0, threshold, top_k);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return SimilarityResult.__wrap(ret[0]);
    }
    /**
     * Number of molecules in the collection.
     * @returns {number}
     */
    size() {
        const ret = wasm.collection_size(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Substructure (SMARTS) search. Returns a `Uint32Array` of the matching
     * molecules' original input positions, capped at `limit` (0 = unlimited).
     * @param {string} smarts
     * @param {number} limit
     * @returns {Uint32Array}
     */
    substructureSearch(smarts, limit) {
        const ptr0 = passStringToWasm0(smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collection_substructureSearch(this.__wbg_ptr, ptr0, len0, limit);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v2;
    }
}
if (Symbol.dispose) Collection.prototype[Symbol.dispose] = Collection.prototype.free;

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
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CollectionBuilderFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_collectionbuilder_free(ptr, 0);
    }
    /**
     * Parse + fingerprint a chunk of SMILES, appending each to the index.
     * Returns the number successfully added; unparseable / unsupported SMILES
     * are skipped and counted in [`CollectionBuilder::failed`].
     * @param {string[]} smiles
     * @returns {number}
     */
    addMany(smiles) {
        const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collectionbuilder_addMany(this.__wbg_ptr, ptr0, len0);
        return ret >>> 0;
    }
    /**
     * Number of molecules accepted so far.
     * @returns {number}
     */
    added() {
        const ret = wasm.collectionbuilder_added(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of SMILES skipped because they failed to parse / encode.
     * @returns {number}
     */
    failed() {
        const ret = wasm.collectionbuilder_failed(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Finalise into a queryable [`Collection`]. Builds the sidecars + rarity
     * table when indexing; otherwise produces a molecules-only collection.
     * @returns {Collection}
     */
    finish() {
        const ptr = this.__destroy_into_raw();
        const ret = wasm.collectionbuilder_finish(ptr);
        return Collection.__wrap(ret);
    }
    /**
     * @param {boolean} build_index
     */
    constructor(build_index) {
        const ret = wasm.collectionbuilder_new(build_index);
        this.__wbg_ptr = ret;
        CollectionBuilderFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
}
if (Symbol.dispose) CollectionBuilder.prototype[Symbol.dispose] = CollectionBuilder.prototype.free;

/**
 * A similarity-search result as two parallel typed arrays: `indices[i]` is the
 * caller's original input position of the i-th hit, `scores[i]` its Tanimoto
 * similarity in `[0, 1]`. Hits are ordered by descending score.
 */
export class SimilarityResult {
    static __wrap(ptr) {
        const obj = Object.create(SimilarityResult.prototype);
        obj.__wbg_ptr = ptr;
        SimilarityResultFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SimilarityResultFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_similarityresult_free(ptr, 0);
    }
    /**
     * Original input positions of the hits (a `Uint32Array` on the JS side).
     * @returns {Uint32Array}
     */
    get indices() {
        const ret = wasm.similarityresult_indices(this.__wbg_ptr);
        var v1 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Tanimoto scores, parallel to `indices` (a `Float32Array` on the JS side).
     * @returns {Float32Array}
     */
    get scores() {
        const ret = wasm.similarityresult_scores(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
}
if (Symbol.dispose) SimilarityResult.prototype[Symbol.dispose] = SimilarityResult.prototype.free;

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
    static __wrap(ptr) {
        const obj = Object.create(SynthonCollection.prototype);
        obj.__wbg_ptr = ptr;
        SynthonCollectionFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SynthonCollectionFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_synthoncollection_free(ptr, 0);
    }
    /**
     * Build a collection from the text of a synthon CSV. Builds the substructure
     * index once here so each search is a cheap view, not a rebuild.
     * @param {string} csv_text
     * @returns {SynthonCollection}
     */
    static fromCsv(csv_text) {
        const ptr0 = passStringToWasm0(csv_text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.synthoncollection_fromCsv(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return SynthonCollection.__wrap(ret[0]);
    }
    /**
     * Nominal product-space size (sum of per-reaction product upper bounds).
     * @returns {bigint}
     */
    numProducts() {
        const ret = wasm.synthoncollection_numProducts(this.__wbg_ptr);
        return BigInt.asUintN(64, ret);
    }
    /**
     * Number of reactions.
     * @returns {number}
     */
    numReactions() {
        const ret = wasm.synthoncollection_numReactions(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of distinct synthons in the pool.
     * @returns {number}
     */
    numSynthons() {
        const ret = wasm.synthoncollection_numSynthons(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Similarity search (Crux ECFP4 over assembled products). Returns the highest
     * `top_k` products with Tanimoto ≥ `cutoff` (`top_k = 0` = all ≥ cutoff),
     * sorted descending.
     * @param {string} query_smiles
     * @param {number} cutoff
     * @param {number} top_k
     * @returns {SynthonSimHits}
     */
    similaritySearch(query_smiles, cutoff, top_k) {
        const ptr0 = passStringToWasm0(query_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.synthoncollection_similaritySearch(this.__wbg_ptr, ptr0, len0, cutoff, top_k);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return SynthonSimHits.__wrap(ret[0]);
    }
    /**
     * Substructure search: the query is parsed as SMILES (element + aromaticity +
     * bond order). `limit` caps the listing (0 = unlimited).
     * @param {string} query_smiles
     * @param {number} limit
     * @returns {SynthonHits}
     */
    substructureSearch(query_smiles, limit) {
        const ptr0 = passStringToWasm0(query_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.synthoncollection_substructureSearch(this.__wbg_ptr, ptr0, len0, limit);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return SynthonHits.__wrap(ret[0]);
    }
}
if (Symbol.dispose) SynthonCollection.prototype[Symbol.dispose] = SynthonCollection.prototype.free;

/**
 * A synthon-space substructure result: parallel `names` / `smiles` arrays, one
 * entry per hit. `names[i]` is the RDKit-compatible product name
 * `synthonId0;…;reactionId`; `smiles[i]` the assembled product SMILES.
 */
export class SynthonHits {
    static __wrap(ptr) {
        const obj = Object.create(SynthonHits.prototype);
        obj.__wbg_ptr = ptr;
        SynthonHitsFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SynthonHitsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_synthonhits_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    len() {
        const ret = wasm.synthonhits_len(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {string[]}
     */
    get names() {
        const ret = wasm.synthonhits_names(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @returns {string[]}
     */
    get smiles() {
        const ret = wasm.synthonhits_smiles(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
}
if (Symbol.dispose) SynthonHits.prototype[Symbol.dispose] = SynthonHits.prototype.free;

/**
 * A synthon-space similarity result: parallel `names` / `smiles` / `scores`
 * arrays, ordered by descending Tanimoto (`scores[i]` in `[0, 1]`).
 */
export class SynthonSimHits {
    static __wrap(ptr) {
        const obj = Object.create(SynthonSimHits.prototype);
        obj.__wbg_ptr = ptr;
        SynthonSimHitsFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SynthonSimHitsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_synthonsimhits_free(ptr, 0);
    }
    /**
     * @returns {number}
     */
    len() {
        const ret = wasm.synthonsimhits_len(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @returns {string[]}
     */
    get names() {
        const ret = wasm.synthonsimhits_names(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @returns {Float32Array}
     */
    get scores() {
        const ret = wasm.synthonsimhits_scores(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @returns {string[]}
     */
    get smiles() {
        const ret = wasm.synthonsimhits_smiles(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
}
if (Symbol.dispose) SynthonSimHits.prototype[Symbol.dispose] = SynthonSimHits.prototype.free;

/**
 * Install a panic hook that forwards Rust panics to `console.error`. Called
 * once at module init on wasm; a no-op on native targets.
 */
export function start() {
    wasm.start();
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_string_get_72bdf95d3ae505b1: function(arg0, arg1) {
            const obj = arg1;
            const ret = typeof(obj) === 'string' ? obj : undefined;
            var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_throw_1506f2235d1bdba0: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_error_a6fa202b58aa1cd3: function(arg0, arg1) {
            let deferred0_0;
            let deferred0_1;
            try {
                deferred0_0 = arg0;
                deferred0_1 = arg1;
                console.error(getStringFromWasm0(arg0, arg1));
            } finally {
                wasm.__wbindgen_free(deferred0_0, deferred0_1, 1);
            }
        },
        __wbg_new_227d7c05414eb861: function() {
            const ret = new Error();
            return ret;
        },
        __wbg_stack_3b0d974bbf31e44f: function(arg0, arg1) {
            const ret = arg1.stack;
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbindgen_cast_0000000000000001: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./crux_wasm_bg.js": import0,
    };
}

const CollectionFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_collection_free(ptr, 1));
const CollectionBuilderFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_collectionbuilder_free(ptr, 1));
const SimilarityResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_similarityresult_free(ptr, 1));
const SynthonCollectionFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_synthoncollection_free(ptr, 1));
const SynthonHitsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_synthonhits_free(ptr, 1));
const SynthonSimHitsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_synthonsimhits_free(ptr, 1));

function addToExternrefTable0(obj) {
    const idx = wasm.__externref_table_alloc();
    wasm.__wbindgen_externrefs.set(idx, obj);
    return idx;
}

function getArrayF32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayJsValueFromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    const mem = getDataViewMemory0();
    const result = [];
    for (let i = ptr; i < ptr + 4 * len; i += 4) {
        result.push(wasm.__wbindgen_externrefs.get(mem.getUint32(i, true)));
    }
    wasm.__externref_drop_slice(ptr, len);
    return result;
}

function getArrayU32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

let cachedFloat32ArrayMemory0 = null;
function getFloat32ArrayMemory0() {
    if (cachedFloat32ArrayMemory0 === null || cachedFloat32ArrayMemory0.byteLength === 0) {
        cachedFloat32ArrayMemory0 = new Float32Array(wasm.memory.buffer);
    }
    return cachedFloat32ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint32ArrayMemory0 = null;
function getUint32ArrayMemory0() {
    if (cachedUint32ArrayMemory0 === null || cachedUint32ArrayMemory0.byteLength === 0) {
        cachedUint32ArrayMemory0 = new Uint32Array(wasm.memory.buffer);
    }
    return cachedUint32ArrayMemory0;
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

function passArrayJsValueToWasm0(array, malloc) {
    const ptr = malloc(array.length * 4, 4) >>> 0;
    for (let i = 0; i < array.length; i++) {
        const add = addToExternrefTable0(array[i]);
        getDataViewMemory0().setUint32(ptr + 4 * i, add, true);
    }
    WASM_VECTOR_LEN = array.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    };
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedDataViewMemory0 = null;
    cachedFloat32ArrayMemory0 = null;
    cachedUint32ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('crux_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
