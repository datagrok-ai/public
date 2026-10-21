/* @ts-self-types="./chem_wasm.d.ts" */

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
     * The engine revision this module was built from (see [`ENGINE_REV`]).
     * @returns {string}
     */
    static engineRev() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.collection_engineRev();
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Restore a collection written by [`Collection::to_bytes`]. Rejects (with
     * an error, never a panic) blobs from another engine revision or format
     * version, and any truncated / corrupted blob (checksum + structural
     * checks of every section).
     * @param {Uint8Array} bytes
     * @returns {Collection}
     */
    static fromBytes(bytes) {
        const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collection_fromBytes(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return Collection.__wrap(ret[0]);
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
     * The persisted-collection container version (see [`PERSIST_FORMAT_VERSION`]).
     * @returns {number}
     */
    static persistFormatVersion() {
        const ret = wasm.collection_persistFormatVersion();
        return ret >>> 0;
    }
    /**
     * Number of scan positions for [`Collection::substructure_search_range`]
     * (== `size`): callers walk `[0, scanLength)` in consecutive ranges.
     * @returns {number}
     */
    scanLength() {
        const ret = wasm.collection_scanLength(this.__wbg_ptr);
        return ret >>> 0;
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
     * Substructure search with a query molecule (RDKit-js's `get_substruct_match(query)`, as
     * Datagrok's Chem searches with the query its `getQueryMolSafe` builds: a `get_qmol`, or a
     * `get_mol` with `mergeQueryHs`): every molecule `get_substruct_match` would match, as
     * [`Collection::substructure_search`] reports them. The query's atoms and bonds match as
     * RDKit matches them (query features, plain atoms by element, charge, isotope and
     * radicals, H atoms merged into counts as the query holds them); stereo is ignored, as
     * `get_substruct_match` ignores it (spike fast-search-rdkit).
     * @param {JSMol} query
     * @param {number} limit
     * @returns {Uint32Array}
     */
    substructureSearchMol(query, limit) {
        _assertClass(query, JSMol);
        const ret = wasm.collection_substructureSearchMol(this.__wbg_ptr, query.__wbg_ptr, limit);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Resumable substructure search over scan positions `[start, end)`
     * (screen-sidecar sorted positions when indexed, input order in direct
     * mode), capped at `limit` hits (0 = unlimited). Walking consecutive
     * ranges from 0 to [`Collection::scan_length`] and concatenating the
     * results (passing the remaining `limit` each time, stopping when it hits
     * 0) yields exactly [`Collection::substructure_search`]'s output — same
     * hits, same order — while letting the host yield between ranges.
     * @param {string} smarts
     * @param {number} start
     * @param {number} end
     * @param {number} limit
     * @returns {Uint32Array}
     */
    substructureSearchRange(smarts, start, end, limit) {
        const ptr0 = passStringToWasm0(smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collection_substructureSearchRange(this.__wbg_ptr, ptr0, len0, start, end, limit);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v2;
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
    /**
     * Serialise an **indexed** collection (CXMOL store, both sidecars, rarity
     * table, original input ids) into one self-validating blob that
     * [`Collection::from_bytes`] restores without re-parsing any SMILES.
     * Layout: `magic[8] | format u32 | flags u32 | rev_len u32 | rev bytes |
     * n u64 | N_SECTIONS × (len u64, bytes) | checksum u64` (all LE; the
     * checksum covers everything before it).
     * @returns {Uint8Array}
     */
    toBytes() {
        const ret = wasm.collection_toBytes(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
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
     * are skipped and counted in [`CollectionBuilder::failed`]. In a build with
     * RDKit-js's `JSMol` API (chem-wasm, crux-wasm's chem module) an item may be a
     * MOL block: read as Datagrok Chem's `getMolSafe` reads it, its H atoms kept.
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
     * Stream SD-file bytes in (any chunking): each record read so far is added like a
     * SMILES (RDKit's MOL-file reading: sanitize, removeHs). Returns the number added.
     * @param {Uint8Array} chunk
     * @returns {number}
     */
    addSdf(chunk) {
        const ptr0 = passArray8ToWasm0(chunk, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.collectionbuilder_addSdf(this.__wbg_ptr, ptr0, len0);
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
     * End of the SD input: reads the last record(s). Returns the number added.
     * @returns {number}
     */
    finishSdf() {
        const ret = wasm.collectionbuilder_finishSdf(this.__wbg_ptr);
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
    /**
     * Read a SMILES RDKit's default read rejects the way RDKit reads it with Kekulize left out
     * of its sanitization (MinimalLib's `get_mol(smiles, {"kekulize": false})`, the retry
     * Datagrok's Chem makes): aromatic rings written without their `[nH]`, aromatic
     * phosphazenes. Off by default.
     * @param {boolean} lenient
     */
    setLenient(lenient) {
        wasm.collectionbuilder_setLenient(this.__wbg_ptr, lenient);
    }
}
if (Symbol.dispose) CollectionBuilder.prototype[Symbol.dispose] = CollectionBuilder.prototype.free;

/**
 * RDKit-js `JSMol`.
 */
export class JSMol {
    static __wrap(ptr) {
        const obj = Object.create(JSMol.prototype);
        obj.__wbg_ptr = ptr;
        JSMolFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        JSMolFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_jsmol_free(ptr, 0);
    }
    /**
     * PORT(minilib.cpp `add_hs`): the MOL block of a copy with its Hs as atoms
     * (`MolOps::addHs(mol, false, addCoords=true)`).
     * @returns {string}
     */
    add_hs() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_add_hs(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp `add_hs_in_place`): `addHs` in place, the stereo perceived again
     * ([`perceived_again`]); `true`.
     * @returns {boolean}
     */
    add_hs_in_place() {
        const ret = wasm.jsmol_add_hs_in_place(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(common.h:894-904): `canonical_tautomer(details)` —
     * `TautomerEnumerator(params).canonicalize`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    canonical_tautomer(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_canonical_tautomer(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(common.h:918-934): `charge_parent(details)` — `chargeParent(mol, params,
     * skipStandardize)`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    charge_parent(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_charge_parent(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(common.h:873-880): `cleanup(details)` — `MolStandardize::cleanup`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    cleanup(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_cleanup(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * `clear_prop(key)`: `true` when it was set.
     * @param {string} key
     * @returns {boolean}
     */
    clear_prop(key) {
        const ptr0 = passStringToWasm0(key, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_clear_prop(this.__wbg_ptr, ptr0, len0);
        return ret !== 0;
    }
    /**
     * PORT(minilib.cpp `convert_to_aromatic_form`): `updatePropertyCache`, `setAromaticity`.
     */
    convert_to_aromatic_form() {
        const ret = wasm.jsmol_convert_to_aromatic_form(this.__wbg_ptr);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * PORT(minilib.cpp `convert_to_kekule_form`): `MolOps::Kekulize` in place (its canonical
     * traversal, `rk::mol_ops::kekulize_canonical`).
     */
    convert_to_kekule_form() {
        const ret = wasm.jsmol_convert_to_kekule_form(this.__wbg_ptr);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * `copy()`: a copy of the molecule.
     * @returns {JSMol}
     */
    copy() {
        const ret = wasm.jsmol_copy(this.__wbg_ptr);
        return JSMol.__wrap(ret);
    }
    /**
     * Embind's `delete()`: frees the molecule (the same as wasm-bindgen's `free()`).
     */
    delete() {
        const ptr = this.__destroy_into_raw();
        wasm.jsmol_delete(ptr);
    }
    /**
     * Crux extension: `MetalDisconnector(MetalDisconnectorOptions).disconnect`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    disconnect_metals(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_disconnect_metals(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `disconnectOrganometallics(mol, MetalDisconnectorOptions)`
     * (the Syngenta defaults, overridden by the details).
     * @param {string | null} [details]
     * @returns {boolean}
     */
    disconnect_organometallics(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_disconnect_organometallics(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(jswrapper.cpp:83-86): `draw_to_canvas(canvas, width, height)`
     * (-1: the canvas's).
     * @param {any} canvas
     * @param {number} width
     * @param {number} height
     * @returns {string}
     */
    draw_to_canvas(canvas, width, height) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_draw_to_canvas(this.__wbg_ptr, canvas, width, height);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(jswrapper.cpp:88-96): `draw_to_canvas_with_highlights(canvas,
     * details)`.
     * @param {any} canvas
     * @param {string} details
     * @returns {string}
     */
    draw_to_canvas_with_highlights(canvas, details) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_draw_to_canvas_with_highlights(this.__wbg_ptr, canvas, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(jswrapper.cpp:67-81): `draw_to_canvas_with_offset`.
     * @param {any} canvas
     * @param {number} offsetx
     * @param {number} offsety
     * @param {number} width
     * @param {number} height
     * @returns {string}
     */
    draw_to_canvas_with_offset(canvas, offsetx, offsety, width, height) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_draw_to_canvas_with_offset(this.__wbg_ptr, canvas, offsetx, offsety, width, height);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(common.h:936-952): `fragment_parent(details)` — `fragmentParent(mol,
     * params, skipStandardize)`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    fragment_parent(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_fragment_parent(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(common.h `generate_aligned_coords`): the molecule laid out with the part matching
     * `template` at the template's coordinates (`generateDepictionMatching2DStructure`):
     * `allowRGroups`, `acceptFailure`, `alignOnly`, `referenceSmarts`, `useCoordGen`
     * (`RDDepict::preferCoordGen` for the call: CoordGen, decision `coordgen`). The match as
     * `{"atoms":[...],"bonds":[...]}`; `""` when the template has no coordinates or nothing
     * matches (the coordinates kept); `"{}"` when `acceptFailure` laid it out freely. A layout
     * CoordGen cannot make yet (a part not ported) throws "not supported by Crux: CoordGen ...",
     * where RDKit-js would have laid it out.
     * @param {JSMol} template
     * @param {string | null} [details]
     * @returns {string}
     */
    generate_aligned_coords(template, details) {
        let deferred3_0;
        let deferred3_1;
        try {
            _assertClass(template, JSMol);
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_generate_aligned_coords(this.__wbg_ptr, template.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp `get_aromatic_form`): the MOL block of the aromatic form.
     * @returns {string}
     */
    get_aromatic_form() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_aromatic_form(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * @param {string | null} [details]
     * @returns {Uint8Array}
     */
    get_atom_pair_fp_as_uint8array(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_atom_pair_fp_as_uint8array(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * MinimalLib `get_atom_pair_fp(details)`.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_atom_pair_fp(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_atom_pair_fp(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * @param {string | null} [details]
     * @returns {Uint8Array}
     */
    get_avalon_fp_as_uint8array(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_avalon_fp_as_uint8array(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * MinimalLib `get_avalon_fp(details)` (`nBits` 512, `avalonSSSBits`, as there); throws in
     * chem-wasm's one file ([`AVALON_AND_MAPS`]).
     * @param {string | null} [details]
     * @returns {string}
     */
    get_avalon_fp(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_avalon_fp(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(SmartsWrite.cpp `MolToCXSmarts`): `get_cxsmarts(details)`: the SMARTS, then
     * `getCXExtensions(mol, CX_ALL)` over the SMARTS writer's atom and bond output order.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_cxsmarts(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_cxsmarts(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:104-110): `get_cxsmiles(details)`: `MolToCXSmiles(mol, params,
     * flags, restoreBondDirs)`; the coordinates written when the molecule has them, as RDKit.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_cxsmiles(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_cxsmiles(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(common.h `get_descriptors`): RDKit's `Descriptors::Properties` as JSON.
     * @returns {string}
     */
    get_descriptors() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_descriptors(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * Crux: any kind (`maccs|rdkit|tt|pattern|morgan|ap|avalon`) in any form, as JSON.
     * @param {string} kind
     * @param {string | null} [details]
     * @returns {string}
     */
    get_fp_json(kind, details) {
        let deferred4_0;
        let deferred4_1;
        try {
            const ptr0 = passStringToWasm0(kind, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            var ptr1 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len1 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_fp_json(this.__wbg_ptr, ptr0, len0, ptr1, len1);
            var ptr3 = ret[0];
            var len3 = ret[1];
            if (ret[3]) {
                ptr3 = 0; len3 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred4_0 = ptr3;
            deferred4_1 = len3;
            return getStringFromWasm0(ptr3, len3);
        } finally {
            wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
        }
    }
    /**
     * PORT(jswrapper.cpp `get_frags_helper`): `{molList, mappings}`: the fragments
     * (`getMolFrags(asMols, sanitizeFrags)`) and `{"frags":[...],"fragsMolAtomMapping":[[...]]}`.
     * `sanitizeFrags` (default true); `copyConformers` only on.
     * @param {string | null} [details]
     * @returns {any}
     */
    get_frags(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_frags(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * `get_inchi(options)`: InChI lives in inchi-wasm (`molToInchi` of `get_molblock()`), the lazy
     * module Chem already ships (decision `chem-wasm`): throws `not supported by Crux`.
     * @param {string | null} [_options]
     * @returns {string}
     */
    get_inchi(_options) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(_options) ? 0 : passStringToWasm0(_options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_inchi(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:158-160): `get_json()`: `MolInterchange::MolToJSONData(mol)`, RDKit's
     * commonchem writer with its extensions (`molfile::mol_interchange`, decision `chem-wasm`).
     * @returns {string}
     */
    get_json() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_json(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp `get_kekule_form`): the MOL block of the Kekulé form.
     * @returns {string}
     */
    get_kekule_form() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_kekule_form(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * @returns {Uint8Array}
     */
    get_maccs_fp_as_uint8array() {
        const ret = wasm.jsmol_get_maccs_fp_as_uint8array(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
    }
    /**
     * MinimalLib `get_maccs_fp()`: the 167 MACCS keys as a bit string.
     * @returns {string}
     */
    get_maccs_fp() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_maccs_fp(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:592-603, jswrapper.cpp:449-458 `get_mmpa_frags`): `{cores, sidechains}`, two `MolList`s
     * of RDKit's `MMPA::fragmentMol(mol, minCuts, maxCuts, maxCutBonds)` (its default pattern) results, empty when
     * no cut has `minCuts` bonds; both `null` where `fragmentMol` returns false (no cuttable bond, or more than
     * `maxCutBonds`). A single cut's core is empty (`null` from `at` / `next`).
     * @param {number} min_cuts
     * @param {number} max_cuts
     * @param {number} max_cut_bonds
     * @returns {any}
     */
    get_mmpa_frags(min_cuts, max_cuts, max_cut_bonds) {
        const ret = wasm.jsmol_get_mmpa_frags(this.__wbg_ptr, min_cuts, max_cuts, max_cut_bonds);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * `get_molblock(details)`: `MolToMolBlock`; with MinimalLib's details
     * (`includeStereo`, `kekulize`, `useMolBlockWedging`, `addChiralHs`; decision
     * `chem-wasm`: [`crate::minilib::molblock_with_details`]).
     * @param {string | null} [details]
     * @returns {string}
     */
    get_molblock(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_molblock(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * @param {string | null} [details]
     * @returns {Uint8Array}
     */
    get_morgan_fp_as_uint8array(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_morgan_fp_as_uint8array(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * MinimalLib `get_morgan_fp(details)` (radius 2 by default, as there).
     * @param {string | null} [details]
     * @returns {string}
     */
    get_morgan_fp(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_morgan_fp(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:396-409): `get_new_coords(useCoordGen)`: the mol block of a copy
     * with new 2D coordinates, `RDDepict::preferCoordGen` = `useCoordGen` for the call.
     * @param {boolean | null} [use_coord_gen]
     * @returns {string}
     */
    get_new_coords(use_coord_gen) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_new_coords(this.__wbg_ptr, isLikeNone(use_coord_gen) ? 0xFFFFFF : use_coord_gen ? 1 : 0);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * `get_num_atoms(heavyOnly)`.
     * @param {boolean | null} [heavy_only]
     * @returns {number}
     */
    get_num_atoms(heavy_only) {
        const ret = wasm.jsmol_get_num_atoms(this.__wbg_ptr, isLikeNone(heavy_only) ? 0xFFFFFF : heavy_only ? 1 : 0);
        return ret >>> 0;
    }
    /**
     * `get_num_bonds()`.
     * @returns {number}
     */
    get_num_bonds() {
        const ret = wasm.jsmol_get_num_bonds(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {string | null} [details]
     * @returns {Uint8Array}
     */
    get_pattern_fp_as_uint8array(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_pattern_fp_as_uint8array(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * MinimalLib `get_pattern_fp(details)`.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_pattern_fp(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_pattern_fp(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * `get_prop(key)`: the molecule property as text (`""` when unset).
     * @param {string} key
     * @returns {string}
     */
    get_prop(key) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ptr0 = passStringToWasm0(key, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_prop(this.__wbg_ptr, ptr0, len0);
            deferred2_0 = ret[0];
            deferred2_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * `get_prop_list(includePrivate, includeComputed)`: the molecule property names (an
     * array; private names, `_`-prefixed, only with `includePrivate`).
     * @param {boolean | null} [include_private]
     * @param {boolean | null} [_include_computed]
     * @returns {string[]}
     */
    get_prop_list(include_private, _include_computed) {
        const ret = wasm.jsmol_get_prop_list(this.__wbg_ptr, isLikeNone(include_private) ? 0xFFFFFF : include_private ? 1 : 0, isLikeNone(_include_computed) ? 0xFFFFFF : _include_computed ? 1 : 0);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @param {string | null} [details]
     * @returns {Uint8Array}
     */
    get_rdkit_fp_as_uint8array(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_rdkit_fp_as_uint8array(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * MinimalLib `get_rdkit_fp(details)`.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_rdkit_fp(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_rdkit_fp(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:112-118): `get_smarts(details)`: `MolToSmarts(mol, params)`
     * (`doIsomericSmiles`, `rootedAtAtom`).
     * @param {string | null} [details]
     * @returns {string}
     */
    get_smarts(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_smarts(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * `get_smiles(details)`: RDKit's canonical SMILES; with MinimalLib's details
     * (`SmilesWriteParams` JSON, decision `chem-wasm`: [`crate::minilib::smiles_write_params`]).
     * @param {string | null} [details]
     * @returns {string}
     */
    get_smiles(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_smiles(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp `get_stereo_tags`): `{"CIP_atoms":[[atom,"(R)"],...],"CIP_bonds":[[begin,end,"(E)"],...]}`:
     * `assignStereochemistry(cleanIt, force, flagPossible)`, then the new CIP labeler
     * (`CIPLabeler::assignCIPLabels`, crate cip-labeler) on that molecule; a possible centre it
     * leaves unlabelled is `(?)`.
     * @returns {string}
     */
    get_stereo_tags() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_stereo_tags(this.__wbg_ptr);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:189-209): `get_substruct_match(qmol, details)`: the first match as
     * `{"atoms":[...],"bonds":[...]}`, `"{}"` for none.
     * @param {JSMol} q
     * @param {string | null} [details]
     * @returns {string}
     */
    get_substruct_match(q, details) {
        let deferred3_0;
        let deferred3_1;
        try {
            _assertClass(q, JSMol);
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_substruct_match(this.__wbg_ptr, q.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:211-233): `get_substruct_matches(qmol, details)`: every match
     * (uniquified, at most 1000) as a JSON array, `"{}"` for none.
     * @param {JSMol} q
     * @param {string | null} [details]
     * @returns {string}
     */
    get_substruct_matches(q, details) {
        let deferred3_0;
        let deferred3_1;
        try {
            _assertClass(q, JSMol);
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_substruct_matches(this.__wbg_ptr, q.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * RDKit-js `get_substruct_matches(qmol)` with the query given as SMARTS
     * (there is no query `JSMol` yet): JSON
     * `[{"atoms":[...],"bonds":[...]}, ...]`, atoms in query-atom order, bonds
     * in query-bond order (PORT(minilib.cpp:189-209) `get_substruct_matches`).
     * @param {string} smarts
     * @returns {string}
     */
    get_substruct_matches_smarts(smarts) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_substruct_matches_smarts(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:130-132): `get_svg(width, height)` (250 x 200 by
     * default).
     * @param {number | null} [width]
     * @param {number | null} [height]
     * @returns {string}
     */
    get_svg(width, height) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_svg(this.__wbg_ptr, isLikeNone(width) ? Number.MAX_SAFE_INTEGER : (width) >> 0, isLikeNone(height) ? Number.MAX_SAFE_INTEGER : (height) >> 0);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:133-138): `get_svg_with_highlights(details)`.
     * @param {string} details
     * @returns {string}
     */
    get_svg_with_highlights(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_svg_with_highlights(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * @param {string | null} [details]
     * @returns {Uint8Array}
     */
    get_topological_torsion_fp_as_uint8array(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_get_topological_torsion_fp_as_uint8array(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v2;
    }
    /**
     * MinimalLib `get_topological_torsion_fp(details)`.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_topological_torsion_fp(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_topological_torsion_fp(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * `get_v2Kmolblock(details)`: `MolToV2KMolBlock` (an error when V2000 cannot hold it).
     * @param {string | null} [details]
     * @returns {string}
     */
    get_v2Kmolblock(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_v2Kmolblock(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:150-153): `get_v3Kmolblock(details)`: `MolToV3KMolBlock`.
     * @param {string | null} [details]
     * @returns {string}
     */
    get_v3Kmolblock(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_v3Kmolblock(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp:518-523): `has_coords()`: 0, 2 or 3.
     * @returns {number}
     */
    has_coords() {
        const ret = wasm.jsmol_has_coords(this.__wbg_ptr);
        return ret;
    }
    /**
     * `has_prop(key)`: the molecule property is set.
     * @param {string} key
     * @returns {boolean}
     */
    has_prop(key) {
        const ptr0 = passStringToWasm0(key, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_has_prop(this.__wbg_ptr, ptr0, len0);
        return ret !== 0;
    }
    /**
     * `is_valid()`.
     * @returns {boolean}
     */
    is_valid() {
        const ret = wasm.jsmol_is_valid(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Crux extension: `isotopeParent(mol, params, skipStandardize)`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    isotope_parent(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_isotope_parent(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `LargestFragmentChooser(preferOrganic, useAtomCount,
     * countHeavyAtomsOnly).choose` (defaults from the params:
     * `preferOrganic`, `largestFragmentChooserUseAtomCount`,
     * `largestFragmentChooserCountHeavyAtomsOnly`).
     * @param {string | null} [details]
     * @returns {boolean}
     */
    largest_fragment(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_largest_fragment(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(common.h:906-916): `neutralize(details)` —
     * `Uncharger(params.doCanonical).uncharge`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    neutralize(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_neutralize(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(minilib.cpp `normalize_depiction`): RDKit's `normalizeDepiction(mol, -1,
     * canonicalize, scaleFactor)` with MinimalLib's canonical transform (its build without
     * Eigen: `depictor::CanonicalTransformBuild::PowerSolver`); the scale factor, or -1
     * (nothing done). Throws where MinimalLib's does (no principal axis).
     * @param {number | null} [canonicalize]
     * @param {number | null} [scale_factor]
     * @returns {number}
     */
    normalize_depiction(canonicalize, scale_factor) {
        const ret = wasm.jsmol_normalize_depiction(this.__wbg_ptr, isLikeNone(canonicalize) ? Number.MAX_SAFE_INTEGER : (canonicalize) >> 0, !isLikeNone(scale_factor), isLikeNone(scale_factor) ? 0 : scale_factor);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0];
    }
    /**
     * PORT(common.h:882-886): `normalize(details)` — `MolStandardize::normalize`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    normalize(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_normalize(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `normalizeInPlace(mol, params)` — RDKit's in-place Normalizer
     * (`runReactant` in place), which can give another result than `normalize`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    normalize_in_place(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_normalize_in_place(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(common.h:888-892): `reionize(details)` — `MolStandardize::reionize`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    reionize(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_reionize(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `removeFragments(mol, params)`; with `leave_last` /
     * `skip_if_all_match` in the details, `FragmentRemover(params' catalog,
     * leave_last, skip_if_all_match).remove` (defaults true / false).
     * @param {string | null} [details]
     * @returns {boolean}
     */
    remove_fragments(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_remove_fragments(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(minilib.cpp `remove_hs`): the MOL block of a copy without its Hs: RDKit's
     * `removeAllHs` (isotopic Hs too, as RDKit-js), or with details `removeHs` with RDKit's
     * defaults (details only `{}`). A copy without coordinates is laid out with the atom ranks
     * the removal kept (RDKit's depictor reads them); the writer then reads the molecule without
     * its computed stereo ([`without_computed_stereo`]), so a potential stereo double bond is
     * written "either" (3), as RDKit-js's `remove_hs()` writes it.
     * @param {string | null} [details]
     * @returns {string}
     */
    remove_hs(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_remove_hs(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * PORT(minilib.cpp `remove_hs_in_place`): [`Self::remove_hs`]'s removal in place, the stereo
     * perceived again ([`perceived_again`]); `true`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    remove_hs_in_place(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_remove_hs_in_place(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(minilib.cpp:383-394): `set_new_coords(useCoordGen)`: `compute2DCoords` (C++
     * defaults: no canonical orientation) with `RDDepict::preferCoordGen` = `useCoordGen` for
     * the call ([`with_prefer_coord_gen`]: CoordGen, crate `coordgen`, decision `coordgen`),
     * under the stereo perception RDKit-js lays out with ([`compute_new_coords`]).
     * @param {boolean | null} [use_coord_gen]
     * @returns {boolean}
     */
    set_new_coords(use_coord_gen) {
        const ret = wasm.jsmol_set_new_coords(this.__wbg_ptr, isLikeNone(use_coord_gen) ? 0xFFFFFF : use_coord_gen ? 1 : 0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * `set_prop(key, value, computed)`: sets a molecule property (text; `_MolFileChiralFlag`
     * as the number the MOL writer reads); `true`.
     * @param {string} key
     * @param {string} value
     * @param {boolean | null} [_computed]
     * @returns {boolean}
     */
    set_prop(key, value, _computed) {
        const ptr0 = passStringToWasm0(key, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(value, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_set_prop(this.__wbg_ptr, ptr0, len0, ptr1, len1, isLikeNone(_computed) ? 0xFFFFFF : _computed ? 1 : 0);
        return ret !== 0;
    }
    /**
     * Crux extension: `stereoParent(mol, params, skipStandardize)`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    stereo_parent(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_stereo_parent(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * PORT(minilib.cpp `straighten_depiction`): RDKit's `straightenDepiction(mol, -1,
     * minimizeRotation)`.
     * @param {boolean | null} [minimize_rotation]
     */
    straighten_depiction(minimize_rotation) {
        wasm.jsmol_straighten_depiction(this.__wbg_ptr, isLikeNone(minimize_rotation) ? 0xFFFFFF : minimize_rotation ? 1 : 0);
    }
    /**
     * Crux extension: `superParent(mol, params, skipStandardize)`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    super_parent(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_super_parent(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `tautomerParent(mol, params, skipStandardize)`.
     * @param {string | null} [details]
     * @returns {boolean}
     */
    tautomer_parent(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_tautomer_parent(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `Uncharger(canonicalOrder, force, protonationOnly).uncharge`
     * (`canonicalOrder` defaults to the params' `doCanonical`).
     * @param {string | null} [details]
     * @returns {boolean}
     */
    uncharge(details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.jsmol_uncharge(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] !== 0;
    }
    /**
     * Crux extension: `ValidationMethod::validate(mol, reportAllFailures)` — a JSON
     * array of the messages. Details: `validations` (RDKit class names, run as a
     * `CompositeValidation`; default `["RDKitValidation"]`), `reportAllFailures`
     * (default false), and the method options (`allowEmptyMolecules`, `strict`,
     * `allowedAtoms` / `disallowedAtoms` as atomic numbers, the
     * `FeaturesValidation` `allow*` flags). RDKitValidation reports valence errors
     * only on an unsanitized molecule: `get_mol(smiles, '{"sanitize":false}')`.
     * @param {string | null} [details]
     * @returns {string}
     */
    validate(details) {
        let deferred3_0;
        let deferred3_1;
        try {
            var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_validate(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
}
if (Symbol.dispose) JSMol.prototype[Symbol.dispose] = JSMol.prototype.free;

/**
 * Embind's `std::vector<JSMolList *>` as `run_reactants` returns it.
 */
export class JSMolListList {
    static __wrap(ptr) {
        const obj = Object.create(JSMolListList.prototype);
        obj.__wbg_ptr = ptr;
        JSMolListListFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        JSMolListListFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_jsmollistlist_free(ptr, 0);
    }
    /**
     * Embind's `delete()`.
     */
    delete() {
        const ptr = this.__destroy_into_raw();
        wasm.jsmollistlist_delete(ptr);
    }
    /**
     * `get(i)`; `undefined` past the end.
     * @param {number} i
     * @returns {MolList | undefined}
     */
    get(i) {
        const ret = wasm.jsmollistlist_get(this.__wbg_ptr, i);
        return ret === 0 ? undefined : MolList.__wrap(ret);
    }
    /**
     * `size()`.
     * @returns {number}
     */
    size() {
        const ret = wasm.jsmollistlist_size(this.__wbg_ptr);
        return ret >>> 0;
    }
}
if (Symbol.dispose) JSMolListList.prototype[Symbol.dispose] = JSMolListList.prototype.free;

/**
 * RDKit-js `MolList` (MinimalLib `JSMolList`). An entry can be empty (`null`, as RDKit-js):
 * `get_mmpa_frags`'s core of a single cut.
 */
export class MolList {
    static __wrap(ptr) {
        const obj = Object.create(MolList.prototype);
        obj.__wbg_ptr = ptr;
        MolListFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    static __unwrap(jsValue) {
        if (!(jsValue instanceof MolList)) {
            return 0;
        }
        return jsValue.__destroy_into_raw();
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MolListFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_mollist_free(ptr, 0);
    }
    /**
     * PORT(minilib.cpp:697-700): `append(mol)` — a copy of `mol`; the new size.
     * @param {JSMol} mol
     * @returns {number}
     */
    append(mol) {
        _assertClass(mol, JSMol);
        const ret = wasm.mollist_append(this.__wbg_ptr, mol.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * PORT(minilib.cpp:674-683): `at(idx)`; `null` past the end or for an empty entry.
     * @param {number} idx
     * @returns {JSMol | null}
     */
    at(idx) {
        const ret = wasm.mollist_at(this.__wbg_ptr, idx);
        return ret;
    }
    /**
     * `at_end()`.
     * @returns {boolean}
     */
    at_end() {
        const ret = wasm.mollist_at_end(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Embind's `delete()` (wasm-bindgen's `free()`).
     */
    delete() {
        const ptr = this.__destroy_into_raw();
        wasm.mollist_delete(ptr);
    }
    /**
     * PORT(minilib.cpp:702-706): `insert(idx, mol)` — at `min(idx, size)`; the new
     * size.
     * @param {number} idx
     * @param {JSMol} mol
     * @returns {number}
     */
    insert(idx, mol) {
        _assertClass(mol, JSMol);
        const ret = wasm.mollist_insert(this.__wbg_ptr, idx, mol.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * `new MolList()`.
     */
    constructor() {
        const ret = wasm.mollist_new();
        this.__wbg_ptr = ret;
        MolListFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * PORT(minilib.cpp:666-672): `next()` — the molecule at the cursor, then the
     * cursor moves on; `null` at the end (or for an empty entry), as RDKit-js.
     * @returns {JSMol | null}
     */
    next() {
        const ret = wasm.mollist_next(this.__wbg_ptr);
        return ret;
    }
    /**
     * PORT(minilib.cpp:685-695): `pop(idx)` — removes and returns the molecule (`null` past the
     * end); the cursor stays on the same next molecule.
     * @param {number} idx
     * @returns {JSMol | null}
     */
    pop(idx) {
        const ret = wasm.mollist_pop(this.__wbg_ptr, idx);
        return ret;
    }
    /**
     * `reset()` — the cursor back to the start.
     */
    reset() {
        wasm.mollist_reset(this.__wbg_ptr);
    }
    /**
     * `size()`.
     * @returns {number}
     */
    size() {
        const ret = wasm.mollist_size(this.__wbg_ptr);
        return ret >>> 0;
    }
}
if (Symbol.dispose) MolList.prototype[Symbol.dispose] = MolList.prototype.free;

/**
 * MinimalLib's `RGroupDecomposition` object.
 */
export class RGroupDecomposition {
    static __wrap(ptr) {
        const obj = Object.create(RGroupDecomposition.prototype);
        obj.__wbg_ptr = ptr;
        RGroupDecompositionFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        RGroupDecompositionFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_rgroupdecomposition_free(ptr, 0);
    }
    /**
     * PORT(minilib.cpp:997-999): RDKit-js `add(mol)` — the molecule's row, -1 when no core
     * matches, -2 when it matches with no side chains.
     * @param {JSMol} mol
     * @returns {number}
     */
    add(mol) {
        _assertClass(mol, JSMol);
        const ret = wasm.rgroupdecomposition_add(this.__wbg_ptr, mol.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0];
    }
    /**
     * `add` of a molecule given as SMILES or a MOL block.
     * @param {string} mol
     * @returns {number}
     */
    add_text(mol) {
        const ptr0 = passStringToWasm0(mol, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.rgroupdecomposition_add_text(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0];
    }
    /**
     * Embind's `delete()` (wasm-bindgen's `free()`).
     */
    delete() {
        const ptr = this.__destroy_into_raw();
        wasm.rgroupdecomposition_delete(ptr);
    }
    /**
     * PORT(jswrapper.cpp:489-500): RDKit-js `get_rgroups_as_columns()` — `{Core: MolList, R1:
     * MolList, ...}`, each piece a `JSMol` with RDKit's props (`_rgroupTargetAtoms`, ...).
     * @returns {any}
     */
    get_rgroups_as_columns() {
        const ret = wasm.rgroupdecomposition_get_rgroups_as_columns(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * `get_rgroups_as_columns()` as canonical SMILES: `{Core: [smiles, ...], R1: [...], ...}`
     * (`""` where a molecule has no group at a label).
     * @returns {any}
     */
    get_rgroups_as_columns_smiles() {
        const ret = wasm.rgroupdecomposition_get_rgroups_as_columns_smiles(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * PORT(jswrapper.cpp:502-518): RDKit-js `get_rgroups_as_rows()` — an array of
     * `{Core: JSMol, R1: JSMol, ...}`.
     * @returns {any}
     */
    get_rgroups_as_rows() {
        const ret = wasm.rgroupdecomposition_get_rgroups_as_rows(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * `get_rgroups_as_rows()` as canonical SMILES: an array of `{Core: smiles, R1: smiles, ...}`.
     * @returns {any}
     */
    get_rgroups_as_rows_smiles() {
        const ret = wasm.rgroupdecomposition_get_rgroups_as_rows_smiles(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * PORT(minilib.cpp:1001): `process()`.
     * @returns {boolean}
     */
    process() {
        const ret = wasm.rgroupdecomposition_process(this.__wbg_ptr);
        return ret !== 0;
    }
}
if (Symbol.dispose) RGroupDecomposition.prototype[Symbol.dispose] = RGroupDecomposition.prototype.free;

/**
 * RDKit-js `Reaction` (MinimalLib `JSReaction`).
 */
export class Reaction {
    static __wrap(ptr) {
        const obj = Object.create(Reaction.prototype);
        obj.__wbg_ptr = ptr;
        ReactionFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ReactionFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_reaction_free(ptr, 0);
    }
    /**
     * Embind's `delete()`.
     */
    delete() {
        const ptr = this.__destroy_into_raw();
        wasm.reaction_delete(ptr);
    }
    /**
     * See [`Self::get_svg`].
     * @param {any} _canvas
     * @param {number} _width
     * @param {number} _height
     * @returns {string}
     */
    draw_to_canvas(_canvas, _width, _height) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.reaction_draw_to_canvas(this.__wbg_ptr, _canvas, _width, _height);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * See [`Self::get_svg`].
     * @param {any} _canvas
     * @param {string} _details
     * @returns {string}
     */
    draw_to_canvas_with_highlights(_canvas, _details) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(_details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.reaction_draw_to_canvas_with_highlights(this.__wbg_ptr, _canvas, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * RDKit `EnumerateLibrary(rxn, reagents, strategy, params)` (not in
     * MinimalLib): an enumerator over `reagentLists` — one array of SMILES per
     * reactant template (unparsable SMILES are skipped). `details` (JSON, all
     * optional): `strategy` (`"cartesian"` default, `"random"`, `"randomAllBBs"`,
     * `"evenPairs"`), `seed` (the random strategies' `minstd_rand`; RDKit's default
     * is 1), `reagentMaxMatchCount`, `sanePartialProducts`, `state` (a
     * `get_state()` text to resume from).
     * @param {any} reagent_lists
     * @param {string | null} [details]
     * @returns {ReactionEnumerator}
     */
    enumerate(reagent_lists, details) {
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.reaction_enumerate(this.__wbg_ptr, reagent_lists, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ReactionEnumerator.__wrap(ret[0]);
    }
    /**
     * `ChemicalReactionToRxnSmarts` (not in MinimalLib).
     * @returns {string}
     */
    get_smarts() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.reaction_get_smarts(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * RDKit-js's reaction drawing (`get_svg`, `get_svg_with_highlights`, `draw_to_canvas*`) is not
     * ported (moldraw2d draws molecules only): these throw `not supported by Crux` rather than
     * being undefined (decision `chem-wasm`).
     * @param {number | null} [_width]
     * @param {number | null} [_height]
     * @returns {string}
     */
    get_svg(_width, _height) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.reaction_get_svg(this.__wbg_ptr, isLikeNone(_width) ? Number.MAX_SAFE_INTEGER : (_width) >> 0, isLikeNone(_height) ? Number.MAX_SAFE_INTEGER : (_height) >> 0);
            var ptr1 = ret[0];
            var len1 = ret[1];
            if (ret[3]) {
                ptr1 = 0; len1 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred2_0 = ptr1;
            deferred2_1 = len1;
            return getStringFromWasm0(ptr1, len1);
        } finally {
            wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
        }
    }
    /**
     * See [`Self::get_svg`].
     * @param {string} _details
     * @returns {string}
     */
    get_svg_with_highlights(_details) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(_details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.reaction_get_svg_with_highlights(this.__wbg_ptr, ptr0, len0);
            var ptr2 = ret[0];
            var len2 = ret[1];
            if (ret[3]) {
                ptr2 = 0; len2 = 0;
                throw takeFromExternrefTable0(ret[2]);
            }
            deferred3_0 = ptr2;
            deferred3_1 = len2;
            return getStringFromWasm0(ptr2, len2);
        } finally {
            wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
        }
    }
    /**
     * `is_valid()`.
     * @returns {boolean}
     */
    is_valid() {
        const ret = wasm.reaction_is_valid(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * PORT(minilib.cpp:644-663): `run_reactants(reactants, maxProducts)` —
     * `initReactantMatchers`, then `runReactants`: one `MolList` of (unsanitized)
     * products per combination of reactant matches. `maxProducts` defaults to
     * 1000 (`JSReaction::maxProducts`).
     * @param {MolList} reactants
     * @param {number | null} [max_products]
     * @returns {JSMolListList}
     */
    run_reactants(reactants, max_products) {
        _assertClass(reactants, MolList);
        const ret = wasm.reaction_run_reactants(this.__wbg_ptr, reactants.__wbg_ptr, isLikeNone(max_products) ? Number.MAX_SAFE_INTEGER : (max_products) >>> 0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return JSMolListList.__wrap(ret[0]);
    }
}
if (Symbol.dispose) Reaction.prototype[Symbol.dispose] = Reaction.prototype.free;

/**
 * RDKit `EnumerateLibrary` ([`JSReaction::enumerate`]; `rxn::enumerate`).
 */
export class ReactionEnumerator {
    static __wrap(ptr) {
        const obj = Object.create(ReactionEnumerator.prototype);
        obj.__wbg_ptr = ptr;
        ReactionEnumeratorFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ReactionEnumeratorFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_reactionenumerator_free(ptr, 0);
    }
    /**
     * Embind's `delete()`.
     */
    delete() {
        const ptr = this.__destroy_into_raw();
        wasm.reactionenumerator_delete(ptr);
    }
    /**
     * `getEnumerator().getNumPermutations()` (`2^64 - 1` when it overflows).
     * @returns {number}
     */
    get_num_permutations() {
        const ret = wasm.reactionenumerator_get_num_permutations(this.__wbg_ptr);
        return ret;
    }
    /**
     * `getEnumerator().getPermutationIdx()`: combinations handed out so far.
     * @returns {number}
     */
    get_permutation_idx() {
        const ret = wasm.reactionenumerator_get_permutation_idx(this.__wbg_ptr);
        return ret;
    }
    /**
     * `getPosition()`: the indices of the last combination into each reagent
     * list the reaction accepts (`get_reagent_counts`).
     * @returns {Float64Array}
     */
    get_position() {
        const ret = wasm.reactionenumerator_get_position(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * `getReagents()` sizes: the reagents kept per reactant template.
     * @returns {Float64Array}
     */
    get_reagent_counts() {
        const ret = wasm.reactionenumerator_get_reagent_counts(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * `getState()`: RDKit's state text (resume with `set_state` or the
     * `state` detail).
     * @returns {string}
     */
    get_state() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.reactionenumerator_get_state(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * `bool(library)`: are there combinations left (the random strategies never
     * halt).
     * @returns {boolean}
     */
    has_next() {
        const ret = wasm.reactionenumerator_has_next(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * `next()`: the product sets of the next reagent combination (one `MolList`
     * per combination of reactant matches, one product per product template);
     * `undefined` when the enumeration is exhausted.
     * @returns {JSMolListList | undefined}
     */
    next() {
        const ret = wasm.reactionenumerator_next(this.__wbg_ptr);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] === 0 ? undefined : JSMolListList.__wrap(ret[0]);
    }
    /**
     * `nextSmiles()`: [`ReactionEnumerator::next`] as a JSON array of arrays of
     * SMILES; `undefined` when exhausted.
     * @returns {string | undefined}
     */
    next_smiles() {
        const ret = wasm.reactionenumerator_next_smiles(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        let v1;
        if (ret[0] !== 0) {
            v1 = getStringFromWasm0(ret[0], ret[1]);
            wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        }
        return v1;
    }
    /**
     * `resetState()`: back to the start.
     */
    reset_state() {
        wasm.reactionenumerator_reset_state(this.__wbg_ptr);
    }
    /**
     * `setState(state)`.
     * @param {string} state
     */
    set_state(state) {
        const ptr0 = passStringToWasm0(state, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.reactionenumerator_set_state(this.__wbg_ptr, ptr0, len0);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
}
if (Symbol.dispose) ReactionEnumerator.prototype[Symbol.dispose] = ReactionEnumerator.prototype.free;

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
 * `calcDescriptors(molecules, names)`: the JSON of the container's `POST /chem/descriptors` (module docs). A
 * `null` (or any non-string) molecule is a row of nulls, as chem.py's `"M  END" in None` raises.
 * @param {(string | null)[]} molecules
 * @param {string[]} names
 * @returns {string}
 */
export function calcDescriptors(molecules, names) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passArrayJsValueToWasm0(names, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.calcDescriptors(molecules, ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * `chemProperties(molecules, names)`: the JSON of Chem's property columns (module docs); throws on a name not in
 * its table.
 * @param {(string | null)[]} molecules
 * @param {string[]} names
 * @returns {string}
 */
export function chemProperties(molecules, names) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passArrayJsValueToWasm0(names, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.chemProperties(molecules, ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Which parts this build holds: `"core"`, `"heavy"`, `"core,heavy"`, ...; the one file Chem loads (feature `chem`)
 * is `"core,heavy,search,descriptors,scaffold-tree"`.
 * @returns {string}
 */
export function chemWasmModules() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.chemWasmModules();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * `descriptorTree()`: the JSON of the container's `GET /chem/descriptors/tree`.
 * @returns {string}
 */
export function descriptorTree() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.descriptorTree();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Crux extension: `enumerateTautomerSmiles(smiles, params)` — a JSON array of the
 * tautomer SMILES.
 * @param {string} smiles
 * @param {string | null} [details]
 * @returns {string}
 */
export function enumerate_tautomer_smiles(smiles, details) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        var ptr1 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len1 = WASM_VECTOR_LEN;
        const ret = wasm.enumerate_tautomer_smiles(ptr0, len0, ptr1, len1);
        var ptr3 = ret[0];
        var len3 = ret[1];
        if (ret[3]) {
            ptr3 = 0; len3 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
    }
}

/**
 * One-shot `RGroupDecompose(cores, mols, asSmiles=True)` as JSON: `{"rows": [...],
 * "unmatched": [mol index, ...], "row_of_mol": [row or null, ...]}`. Cores are
 * SMILES / MOL blocks (or SMARTS with `cores_are_smarts`), molecules SMILES / MOL
 * blocks; `details` is RDKit's RGD parameter JSON.
 * @param {string[]} cores
 * @param {string[]} mols
 * @param {string | null} [details]
 * @param {boolean | null} [cores_are_smarts]
 * @returns {string}
 */
export function getRGroupDecomposition(cores, mols, details, cores_are_smarts) {
    let deferred5_0;
    let deferred5_1;
    try {
        const ptr0 = passArrayJsValueToWasm0(cores, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayJsValueToWasm0(mols, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        var ptr2 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len2 = WASM_VECTOR_LEN;
        const ret = wasm.getRGroupDecomposition(ptr0, len0, ptr1, len1, ptr2, len2, isLikeNone(cores_are_smarts) ? 0xFFFFFF : cores_are_smarts ? 1 : 0);
        var ptr4 = ret[0];
        var len4 = ret[1];
        if (ret[3]) {
            ptr4 = 0; len4 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred5_0 = ptr4;
        deferred5_1 = len4;
        return getStringFromWasm0(ptr4, len4);
    } finally {
        wasm.__wbindgen_free(deferred5_0, deferred5_1, 1);
    }
}

/**
 * RDKit-js `get_inchikey_for_inchi(inchi)`: in inchi-wasm (`inchiToInchiKey`): throws.
 * @param {string} _inchi
 * @returns {string}
 */
export function get_inchikey_for_inchi(_inchi) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(_inchi, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.get_inchikey_for_inchi(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * `get_mces(mol1, mol2, details)`: RDKit's `rascalMCES(mol1, mol2, opts)`: the Maximum
 * Common Edge Substructure of two molecules (the largest set of bonds they share, possibly
 * in several fragments), when their similarity estimate reaches `similarityThreshold`.
 *
 * `details`: a JSON object of `RascalOptions` fields by RDKit's names
 * (`similarityThreshold`, `completeAromaticRings`, `ringMatchesRingOnly`,
 * `completeSmallestRings`, `exactConnectionsMatch`, `singleLargestFrag`, `minFragSize`,
 * `maxFragSeparation`, `allBestMCESs`, `maxBestMCESs`, `timeout`, `doEquivBondPruning`,
 * `returnEmptyMCES`, `maxBondMatchPairs`, `equivalentAtoms`, `ignoreBondOrders`,
 * `ignoreAtomAromaticity`, `minCliqueSize`); each defaults to RDKit's except
 * `similarityThreshold`, 0.5 here (RDKit: 0.7). Other keys are ignored. `timeout` (seconds,
 * 60 by default, -1 for none) runs on the page's clock.
 *
 * Returns a JSON array of results, as RDKit's `FindMCES` returns a list: empty when the
 * pair is below the threshold (unless `returnEmptyMCES`), one result by default, every best
 * MCES with `allBestMCESs`. Each result, with the Python wrapper's names:
 * `{"smartsString", "atomMatches": [[atom of mol1, atom of mol2], ...], "bondMatches": [[bond
 * of mol1, bond of mol2], ...], "similarity", "tier1Sim", "tier2Sim", "timedOut",
 * "numFragments", "largestFragmentSize"}` (the indices each molecule's own: what a
 * highlight takes). Throws on details that are not JSON or a bad `equivalentAtoms`.
 * @param {JSMol} mol1
 * @param {JSMol} mol2
 * @param {string | null} [details]
 * @returns {string}
 */
export function get_mces(mol1, mol2, details) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(mol1, JSMol);
        _assertClass(mol2, JSMol);
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.get_mces(mol1.__wbg_ptr, mol2.__wbg_ptr, ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * RDKit-js `get_mcs_as_json(mols, details)`: `{"smarts","canceled","numAtoms","numBonds"}`,
 * in RDKit-js's key order.
 * @param {MolList} mols
 * @param {string | null} [details]
 * @returns {string}
 */
export function get_mcs_as_json(mols, details) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(mols, MolList);
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.get_mcs_as_json(mols.__wbg_ptr, ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * RDKit-js `get_mcs_as_mol(mols, details)`: the MCS as a query molecule (`get_qmol` of its
 * SMARTS); `null` when there is none.
 * @param {MolList} mols
 * @param {string | null} [details]
 * @returns {JSMol | null}
 */
export function get_mcs_as_mol(mols, details) {
    _assertClass(mols, MolList);
    var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.get_mcs_as_mol(mols.__wbg_ptr, ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * RDKit-js `get_mcs_as_smarts(mols, details)`: the MCS's SMARTS.
 * @param {MolList} mols
 * @param {string | null} [details]
 * @returns {string}
 */
export function get_mcs_as_smarts(mols, details) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(mols, MolList);
        var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.get_mcs_as_smarts(mols.__wbg_ptr, ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * RDKit-js `get_mol(input, details)` (decision `chem-wasm`): `null` when the input does
 * not parse, as RDKit-js 2026.03.6; throws for details that are not JSON or that ask
 * what Crux cannot do (`not supported by Crux: <key>`). Details keys: [`MolInputDetails`].
 * @param {string} input
 * @param {string | null} [details]
 * @returns {JSMol | null}
 */
export function get_mol(input, details) {
    const ptr0 = passStringToWasm0(input, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.get_mol(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * PORT(minilib.cpp:827-829): `get_mol_copy(other)`.
 * @param {JSMol} other
 * @returns {JSMol}
 */
export function get_mol_copy(other) {
    _assertClass(other, JSMol);
    const ret = wasm.get_mol_copy(other.__wbg_ptr);
    return JSMol.__wrap(ret);
}

/**
 * PORT(common.h `qmol_from_input`): RDKit-js `get_qmol(input)`: a SMARTS (`SmartsToMol`) or
 * a MOL block read as a query (`MolBlockToMol(sanitize=false, removeHs=false)`, its query
 * atoms and bonds); `null` when it does not parse.
 * @param {string} input
 * @returns {JSMol | null}
 */
export function get_qmol(input) {
    const ptr0 = passStringToWasm0(input, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.get_qmol(ptr0, len0);
    return ret;
}

/**
 * PORT(jswrapper.cpp:464-487): RDKit-js `get_rgd(core, details)`: the decomposition around
 * the core molecule (a `get_mol` / `get_qmol` `JSMol`).
 * @param {JSMol} core
 * @param {string | null} [details]
 * @returns {RGroupDecomposition}
 */
export function get_rgd(core, details) {
    _assertClass(core, JSMol);
    var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.get_rgd(core.__wbg_ptr, ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return RGroupDecomposition.__wrap(ret[0]);
}

/**
 * `get_rgd` over SMARTS cores (MinimalLib: `get_rgd(get_qmol(smarts), details)`).
 * @param {string[]} cores
 * @param {string | null} [details]
 * @returns {RGroupDecomposition}
 */
export function get_rgd_smarts(cores, details) {
    const ptr0 = passArrayJsValueToWasm0(cores, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.get_rgd_smarts(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return RGroupDecomposition.__wrap(ret[0]);
}

/**
 * `get_rgd` of one or more cores given as text (SMILES or MOL block each), `details` =
 * RDKit's RGD parameter JSON (`updateRGroupDecompositionParametersFromJSON`).
 * @param {string[]} cores
 * @param {string | null} [details]
 * @returns {RGroupDecomposition}
 */
export function get_rgd_text(cores, details) {
    const ptr0 = passArrayJsValueToWasm0(cores, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.get_rgd_text(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return RGroupDecomposition.__wrap(ret[0]);
}

/**
 * PORT(minilib.cpp:856-859): `get_rxn(input, details)` — `undefined` when the
 * input does not read.
 * @param {string} input
 * @param {string | null} [details]
 * @returns {Reaction | undefined}
 */
export function get_rxn(input, details) {
    const ptr0 = passStringToWasm0(input, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.get_rxn(ptr0, len0, ptr1, len1);
    return ret === 0 ? undefined : Reaction.__wrap(ret);
}

/**
 * PORT(jswrapper.cpp:156-212, 921-928): RDKit-js `molzip`, its three forms and Crux's one-molecule form (module
 * docs): `molzip(a, b, details?)`, `molzip(row, details?)`, `molzip(mol, details?)`. A `JSMol`, or `null`.
 * @param {JSMol | Record<string, JSMol>} first
 * @param {JSMol | string} [second]
 * @param {string} [details]
 * @returns {JSMol | null}
 */
export function molzip(first, second, details) {
    var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.molzip(first, isLikeNone(second) ? 0 : addToExternrefTable0(second), ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Python's `Chem.molzipFragments(mols, params)` (RDKit's `molzip(decomposition, params)`, which RDKit-js does not
 * bind): the list's molecules, the core first. A `JSMol`, or `null` for an empty list.
 * @param {MolList} mols
 * @param {string | null} [details]
 * @returns {JSMol | null}
 */
export function molzipFragments(mols, details) {
    _assertClass(mols, MolList);
    var ptr0 = isLikeNone(details) ? 0 : passStringToWasm0(details, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.molzipFragments(mols.__wbg_ptr, ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * PORT(minilib.cpp:864-868): RDKit-js `prefer_coordgen(prefer)`: `RDDepict::preferCoordGen`,
 * process-wide (CoordGen, crate `coordgen`, decision `coordgen`); `true` throws in a build without
 * CoordGen.
 * @param {boolean} prefer
 */
export function prefer_coordgen(prefer) {
    const ret = wasm.prefer_coordgen(prefer);
    if (ret[1]) {
        throw takeFromExternrefTable0(ret[0]);
    }
}

/**
 * `scaffoldTree(smiles, ringCutoff, discharge)`: the JSON string Datagrok Chem's `scaffold_tree.py` returns for
 * these SMILES (`GenerateScaffoldTree(data, smilesColumn, ringCutoff, dischargeAndDeradicalize)`); throws where
 * the script fails.
 * @param {string[]} smiles
 * @param {number} ring_cutoff
 * @param {boolean} discharge
 * @returns {string}
 */
export function scaffoldTree(smiles, ring_cutoff, discharge) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.scaffoldTree(ptr0, len0, ring_cutoff, discharge);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * RDKit-js `set_log_capture(logName)` / `set_log_tee(logName)`: Crux keeps no RDKit log (its
 * messages are not RDKit's), so these throw `not supported by Crux` rather than being undefined.
 * @param {string} log_name
 * @returns {any}
 */
export function set_log_capture(log_name) {
    const ptr0 = passStringToWasm0(log_name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.set_log_capture(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * See [`set_log_capture`].
 * @param {string} log_name
 * @returns {any}
 */
export function set_log_tee(log_name) {
    const ptr0 = passStringToWasm0(log_name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.set_log_tee(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * RDKit `SimilarityMaps.GetSimilarityMapForFingerprint(ref, probe, fp, draw2d)`
 * on a `MolDraw2DSVG(width, height)`: the SVG text.
 * @param {JSMol} reference
 * @param {JSMol} probe
 * @param {string | null} [options]
 * @returns {string}
 */
export function similarityMapSvg(reference, probe, options) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(reference, JSMol);
        _assertClass(probe, JSMol);
        var ptr0 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.similarityMapSvg(reference.__wbg_ptr, probe.__wbg_ptr, ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * RDKit `SimilarityMaps.GetAtomicWeightsForFingerprint(ref, probe, fp,
 * metric)` as JSON `{"weights":[...],"max_weight":x}` (raw weights; the
 * map draws them divided by `max_weight`).
 * @param {JSMol} reference
 * @param {JSMol} probe
 * @param {string | null} [options]
 * @returns {string}
 */
export function similarityMapWeights(reference, probe, options) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(reference, JSMol);
        _assertClass(probe, JSMol);
        var ptr0 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.similarityMapWeights(reference.__wbg_ptr, probe.__wbg_ptr, ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Crux extension: `standardizeSmiles(smiles)` — the SMILES read unsanitized,
 * `cleanup` with the default parameters, `MolToSmiles`.
 * @param {string} smiles
 * @returns {string}
 */
export function standardize_smiles(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.standardize_smiles(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Panics as `Error`s with their message in the console (not a bare `unreachable`).
 */
export function start() {
    wasm.start();
}

/**
 * RDKit-js `use_legacy_stereo_perception(value)`: RDKit's legacy stereo perception on or
 * off (process-wide, legacy by default as in RDKit); the value it had.
 * @param {boolean} value
 * @returns {boolean}
 */
export function use_legacy_stereo_perception(value) {
    const ret = wasm.use_legacy_stereo_perception(value);
    return ret !== 0;
}

/**
 * Crux extension: `validateSmiles(smiles)` — a JSON array of the messages.
 * @param {string} smiles
 * @returns {string}
 */
export function validate_smiles(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.validate_smiles(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * RDKit-js `version()`: the RDKit release whose behaviour Crux reproduces (`RDKIT_VERSION`).
 * @returns {string}
 */
export function version() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.version();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg_Error_30c8987f7c2ed4e2: function(arg0, arg1) {
            const ret = Error(getStringFromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_Function_a437dbcf1b80ea5e: function() { return handleError(function (arg0, arg1, arg2, arg3) {
            const ret = Function(getStringFromWasm0(arg0, arg1), getStringFromWasm0(arg2, arg3));
            return ret;
        }, arguments); },
        __wbg_String_8564e559799eccda: function(arg0, arg1) {
            const ret = String(arg1);
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_boolean_get_5b446f51afd21013: function(arg0) {
            const v = arg0;
            const ret = typeof(v) === 'boolean' ? v : undefined;
            return isLikeNone(ret) ? 0xFFFFFF : ret ? 1 : 0;
        },
        __wbg___wbindgen_debug_string_4687d8d8c2017d52: function(arg0, arg1) {
            const ret = debugString(arg1);
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_is_function_1f9d30630b8b1d3d: function(arg0) {
            const ret = typeof(arg0) === 'function';
            return ret;
        },
        __wbg___wbindgen_is_object_3c45d4f2dde4e749: function(arg0) {
            const val = arg0;
            const ret = typeof(val) === 'object' && val !== null;
            return ret;
        },
        __wbg___wbindgen_is_string_90b56bc79aad6f6c: function(arg0) {
            const ret = typeof(arg0) === 'string';
            return ret;
        },
        __wbg___wbindgen_is_undefined_8865fb403f8fe9d8: function(arg0) {
            const ret = arg0 === undefined;
            return ret;
        },
        __wbg___wbindgen_jsval_loose_eq_677f21e468d6b461: function(arg0, arg1) {
            const ret = arg0 == arg1;
            return ret;
        },
        __wbg___wbindgen_number_get_2e0e7dee9f701a71: function(arg0, arg1) {
            const obj = arg1;
            const ret = typeof(obj) === 'number' ? obj : undefined;
            getDataViewMemory0().setFloat64(arg0 + 8 * 1, isLikeNone(ret) ? 0 : ret, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, !isLikeNone(ret), true);
        },
        __wbg___wbindgen_string_get_0380ccaa2f57f0d9: function(arg0, arg1) {
            const obj = arg1;
            const ret = typeof(obj) === 'string' ? obj : undefined;
            var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_throw_41e9ee4f547fc59a: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_apply_a910804df6e1e433: function() { return handleError(function (arg0, arg1, arg2) {
            const ret = arg0.apply(arg1, arg2);
            return ret;
        }, arguments); },
        __wbg_call_187d372bd5fdd4aa: function() { return handleError(function (arg0, arg1, arg2) {
            const ret = arg0.call(arg1, arg2);
            return ret;
        }, arguments); },
        __wbg_call_6137034ef55c9d0f: function() { return handleError(function (arg0, arg1) {
            const ret = arg0.call(arg1);
            return ret;
        }, arguments); },
        __wbg_call_939a2607c4484b0b: function() { return handleError(function (arg0, arg1, arg2, arg3, arg4) {
            const ret = arg0.call(arg1, arg2, arg3, arg4);
            return ret;
        }, arguments); },
        __wbg_done_b41a1d26cdb37fb6: function(arg0) {
            const ret = arg0.done;
            return ret;
        },
        __wbg_error_757e9472f8410341: function(arg0, arg1) {
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
        __wbg_get_31af05bd4842a84f: function() { return handleError(function (arg0, arg1) {
            const ret = Reflect.get(arg0, arg1);
            return ret;
        }, arguments); },
        __wbg_get_658f6698067d9515: function() { return handleError(function (arg0, arg1) {
            const ret = Reflect.get(arg0, arg1);
            return ret;
        }, arguments); },
        __wbg_get_unchecked_288889d017702237: function(arg0, arg1) {
            const ret = arg0[arg1 >>> 0];
            return ret;
        },
        __wbg_instanceof_ArrayBuffer_a99f175873e5d9b8: function(arg0) {
            let result;
            try {
                result = arg0 instanceof ArrayBuffer;
            } catch (_) {
                result = false;
            }
            const ret = result;
            return ret;
        },
        __wbg_instanceof_Object_67a83cdc00c5d141: function(arg0) {
            let result;
            try {
                result = arg0 instanceof Object;
            } catch (_) {
                result = false;
            }
            const ret = result;
            return ret;
        },
        __wbg_instanceof_Uint8Array_828cef2aaacafc31: function(arg0) {
            let result;
            try {
                result = arg0 instanceof Uint8Array;
            } catch (_) {
                result = false;
            }
            const ret = result;
            return ret;
        },
        __wbg_isArray_e15a2ff68ffdbef2: function(arg0) {
            const ret = Array.isArray(arg0);
            return ret;
        },
        __wbg_iterator_e3c31c892080e444: function() {
            const ret = Symbol.iterator;
            return ret;
        },
        __wbg_jsmol_new: function(arg0) {
            const ret = JSMol.__wrap(arg0);
            return ret;
        },
        __wbg_keys_440be172f17c0265: function(arg0) {
            const ret = Object.keys(arg0);
            return ret;
        },
        __wbg_length_7f3c00c40364105e: function(arg0) {
            const ret = arg0.length;
            return ret;
        },
        __wbg_length_d4bdea10311bd9cf: function(arg0) {
            const ret = arg0.length;
            return ret;
        },
        __wbg_mollist_new: function(arg0) {
            const ret = MolList.__wrap(arg0);
            return ret;
        },
        __wbg_mollist_unwrap: function(arg0) {
            const ret = MolList.__unwrap(arg0);
            return ret;
        },
        __wbg_new_1dbf7428bba60a42: function(arg0) {
            const ret = new Uint8Array(arg0);
            return ret;
        },
        __wbg_new_227d7c05414eb861: function() {
            const ret = new Error();
            return ret;
        },
        __wbg_new_28744009d011f847: function() {
            const ret = new Map();
            return ret;
        },
        __wbg_new_343a093a3c2ffb4e: function(arg0, arg1) {
            const ret = new Error(getStringFromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_new_617a8cdb8bb1130e: function() {
            const ret = new Object();
            return ret;
        },
        __wbg_new_ee2291f50781bf1d: function() {
            const ret = new Array();
            return ret;
        },
        __wbg_new_from_slice_37a3c3a194266cb8: function(arg0, arg1) {
            const ret = new Uint32Array(getArrayU32FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_new_from_slice_f545fd22ddc142b8: function(arg0, arg1) {
            const ret = new Float64Array(getArrayF64FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_next_33784799010f1bbe: function(arg0) {
            const ret = arg0.next;
            return ret;
        },
        __wbg_next_f4aac29c42af995c: function() { return handleError(function (arg0) {
            const ret = arg0.next();
            return ret;
        }, arguments); },
        __wbg_now_aa4ccb83129e9e55: function() {
            const ret = Date.now();
            return ret;
        },
        __wbg_prototypesetcall_bc27214492979395: function(arg0, arg1, arg2) {
            Uint8Array.prototype.set.call(getArrayU8FromWasm0(arg0, arg1), arg2);
        },
        __wbg_push_2baf45db356cf468: function(arg0, arg1) {
            const ret = arg0.push(arg1);
            return ret;
        },
        __wbg_set_145a351398b48c65: function() { return handleError(function (arg0, arg1, arg2) {
            const ret = Reflect.set(arg0, arg1, arg2);
            return ret;
        }, arguments); },
        __wbg_set_6ae97e73113c4f0b: function(arg0, arg1, arg2) {
            const ret = arg0.set(arg1, arg2);
            return ret;
        },
        __wbg_set_6be42768c690e380: function(arg0, arg1, arg2) {
            arg0[arg1] = arg2;
        },
        __wbg_set_bea140a88be9b277: function(arg0, arg1, arg2) {
            arg0[arg1 >>> 0] = arg2;
        },
        __wbg_stack_3b0d974bbf31e44f: function(arg0, arg1) {
            const ret = arg1.stack;
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg_value_f3c585ee8f5ba40c: function(arg0) {
            const ret = arg0.value;
            return ret;
        },
        __wbindgen_generic_0000000000000001: function(arg0) {
            // Cast intrinsic for `F64 -> Externref`.
            const ret = arg0;
            return ret;
        },
        __wbindgen_generic_0000000000000002: function(arg0) {
            // Cast intrinsic for `I64 -> Externref`.
            const ret = arg0;
            return ret;
        },
        __wbindgen_generic_0000000000000003: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_generic_0000000000000004: function(arg0) {
            // Cast intrinsic for `U64 -> Externref`.
            const ret = BigInt.asUintN(64, arg0);
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
        "./chem_wasm_bg.js": import0,
    };
}

const CollectionFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_collection_free(ptr, 1));
const CollectionBuilderFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_collectionbuilder_free(ptr, 1));
const JSMolFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_jsmol_free(ptr, 1));
const JSMolListListFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_jsmollistlist_free(ptr, 1));
const MolListFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mollist_free(ptr, 1));
const RGroupDecompositionFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_rgroupdecomposition_free(ptr, 1));
const ReactionFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_reaction_free(ptr, 1));
const ReactionEnumeratorFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_reactionenumerator_free(ptr, 1));
const SimilarityResultFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_similarityresult_free(ptr, 1));

function addToExternrefTable0(obj) {
    const idx = wasm.__externref_table_alloc();
    wasm.__wbindgen_externrefs.set(idx, obj);
    return idx;
}

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
}

function debugString(val) {
    // primitive types
    const type = typeof val;
    if (type == 'number' || type == 'boolean' || val == null) {
        return  `${val}`;
    }
    if (type == 'string') {
        return `"${val}"`;
    }
    if (type == 'symbol') {
        const description = val.description;
        if (description == null) {
            return 'Symbol';
        } else {
            return `Symbol(${description})`;
        }
    }
    if (type == 'function') {
        const name = val.name;
        if (typeof name == 'string' && name.length > 0) {
            return `Function(${name})`;
        } else {
            return 'Function';
        }
    }
    // objects
    if (Array.isArray(val)) {
        const length = val.length;
        let debug = '[';
        if (length > 0) {
            debug += debugString(val[0]);
        }
        for(let i = 1; i < length; i++) {
            debug += ', ' + debugString(val[i]);
        }
        debug += ']';
        return debug;
    }
    // Test for built-in
    const builtInMatches = /\[object ([^\]]+)\]/.exec(toString.call(val));
    let className;
    if (builtInMatches && builtInMatches.length > 1) {
        className = builtInMatches[1];
    } else {
        // Failed to match the standard '[object ClassName]'
        return toString.call(val);
    }
    if (className == 'Object') {
        // we're a user defined class or Object
        // JSON.stringify avoids problems with cycles, and is generally much
        // easier than looping through ownProperties of `val`.
        try {
            return 'Object(' + JSON.stringify(val) + ')';
        } catch (_) {
            return 'Object';
        }
    }
    // errors
    if (val instanceof Error) {
        return `${val.name}: ${val.message}\n${val.stack}`;
    }
    // TODO we could test for more things here, like `Set`s and `Map`s.
    return className;
}

function getArrayF32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
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

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
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

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
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

function handleError(f, args) {
    try {
        return f.apply(this, args);
    } catch (e) {
        const idx = addToExternrefTable0(e);
        wasm.__wbindgen_exn_store(idx);
    }
}

function isLikeNone(x) {
    return x === undefined || x === null;
}

function passArray8ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 1, 1) >>> 0;
    getUint8ArrayMemory0().set(arg, ptr / 1);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
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
    cachedFloat64ArrayMemory0 = null;
    cachedUint32ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (!module.ok) {
            throw new Error(`failed to fetch Wasm: ${module.status} ${module.statusText} fetching '${module.url}'`);
        }

        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = expectedResponseType(module.type);

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
        module_or_path = new URL('chem_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
