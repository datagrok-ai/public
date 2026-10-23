/* @ts-self-types="./crux_wasm.d.ts" */

/**
 * A CSV held in this wasm instance's memory (copied in once; a leading BOM dropped).
 */
export class CsvInput {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        CsvInputFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_csvinput_free(ptr, 0);
    }
    /**
     * Pass 2 of verified chunk `k` (on a worker) → an Arrow IPC stream of its batches.
     * @param {any} options
     * @param {any} verified
     * @param {number} k
     * @returns {Uint8Array}
     */
    buildChunk(options, verified, k) {
        const ret = wasm.csvinput_buildChunk(this.__wbg_ptr, options, verified, k);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
    }
    /**
     * Pass 1 of chunk `k` (on a worker).
     * @param {any} options
     * @param {any} plan
     * @param {number} k
     * @returns {any}
     */
    inferChunk(options, plan, k) {
        const ret = wasm.csvinput_inferChunk(this.__wbg_ptr, options, plan, k);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * Bytes held (after the BOM).
     * @returns {number}
     */
    get length() {
        const ret = wasm.csvinput_length(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {Uint8Array} bytes
     */
    constructor(bytes) {
        const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.csvinput_new(ptr0, len0);
        this.__wbg_ptr = ret;
        CsvInputFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Header + chunk guesses (`null` for an input without rows). `engine`: "Fast" (default) or
     * "Compat" (after `verify` returned `null`).
     * @param {any} options
     * @param {number} n_chunks
     * @param {string | null} [engine]
     * @returns {any}
     */
    plan(options, n_chunks, engine) {
        var ptr0 = isLikeNone(engine) ? 0 : passStringToWasm0(engine, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.csvinput_plan(this.__wbg_ptr, options, n_chunks, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
    /**
     * Checks the chunks (re-running misplaced ones here) and merges the types: the verified plan
     * for `buildChunk`, or `null` (non-standard quoting: plan again with "Compat").
     * @param {any} options
     * @param {any} plan
     * @param {any} infos
     * @returns {any}
     */
    verify(options, plan, infos) {
        const ret = wasm.csvinput_verify(this.__wbg_ptr, options, plan, infos);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return takeFromExternrefTable0(ret[0]);
    }
}
if (Symbol.dispose) CsvInput.prototype[Symbol.dispose] = CsvInput.prototype.free;

/**
 * The result of [`descriptors_table`]: the requested columns (in request order)
 * and the per-row `error`, as a `RecordBatch`.
 */
export class DescriptorTable {
    static __wrap(ptr) {
        const obj = Object.create(DescriptorTable.prototype);
        obj.__wbg_ptr = ptr;
        DescriptorTableFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        DescriptorTableFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_descriptortable_free(ptr, 0);
    }
    /**
     * Column `i` as numbers (null → NaN; see [`DescriptorTable::nulls`]).
     * @param {number} i
     * @returns {Float64Array}
     */
    column(i) {
        const ret = wasm.descriptortable_column(this.__wbg_ptr, i);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Per row: why it has no values (`""` when it has).
     * @returns {string[]}
     */
    errors() {
        const ret = wasm.descriptortable_errors(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * The descriptor column names, in request order (without `error`).
     * @returns {string[]}
     */
    names() {
        const ret = wasm.descriptortable_names(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Column `i`'s null mask (1 = null: RDKit has no value, or the row failed).
     * @param {number} i
     * @returns {Uint8Array}
     */
    nulls(i) {
        const ret = wasm.descriptortable_nulls(this.__wbg_ptr, i);
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
    }
    /**
     * @returns {number}
     */
    get numRows() {
        const ret = wasm.descriptortable_numRows(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * The table as Arrow IPC stream bytes.
     * @returns {Uint8Array}
     */
    toIpc() {
        const ret = wasm.descriptortable_toIpc(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
    }
    /**
     * `"double"` / `"bigint"` per column.
     * @returns {string[]}
     */
    types() {
        const ret = wasm.descriptortable_types(this.__wbg_ptr);
        var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
}
if (Symbol.dispose) DescriptorTable.prototype[Symbol.dispose] = DescriptorTable.prototype.free;

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
     * `get_autocorr2d()`: `rdMolDescriptors.CalcAUTOCORR2D(mol)` (192 values,
     * `AUTOCORR2D_1..192`).
     * @returns {Float64Array}
     */
    get_autocorr2d() {
        const ret = wasm.jsmol_get_autocorr2d(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
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
     * `get_descriptors_3d()`: `Descriptors3D.CalcMolDescriptors3D(mol)` on the
     * molecule's conformer, in [`descriptors_3d_names`] order. Throws without one.
     * @returns {Float64Array}
     */
    get_descriptors_3d() {
        const ret = wasm.jsmol_get_descriptors_3d(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * `get_descriptors_3d_vectors()`: `CalcAUTOCORR3D` (80), `CalcRDF` (210),
     * `CalcMORSE` (224) and `CalcWHIM` (114, `thresh = 0.001`) on the molecule's
     * conformer, concatenated (628 values). Throws without one.
     * @returns {Float64Array}
     */
    get_descriptors_3d_vectors() {
        const ret = wasm.jsmol_get_descriptors_3d_vectors(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * `get_efgs_json()`: the molecule's functional groups (module docs).
     * @returns {string}
     */
    get_efgs_json() {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_efgs_json(this.__wbg_ptr);
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
     * `get_efgs_svg(width, height)`: the molecule drawn with its functional
     * groups coloured as in Ertl's paper (EFGs `col_mol`; 300 x 250 by default).
     * @param {number | null} [width]
     * @param {number | null} [height]
     * @returns {string}
     */
    get_efgs_svg(width, height) {
        let deferred2_0;
        let deferred2_1;
        try {
            const ret = wasm.jsmol_get_efgs_svg(this.__wbg_ptr, isLikeNone(width) ? Number.MAX_SAFE_INTEGER : (width) >> 0, isLikeNone(height) ? Number.MAX_SAFE_INTEGER : (height) >> 0);
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
     * `get_gasteiger_charges()`: `ComputeGasteigerCharges(mol)` (12 iterations),
     * each atom's `_GasteigerCharge`.
     * @returns {Float64Array}
     */
    get_gasteiger_charges() {
        const ret = wasm.jsmol_get_gasteiger_charges(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * `get_gasteiger_h_charges()`: each atom's `_GasteigerHCharge` (the summed
     * charge of its implicit hydrogens).
     * @returns {Float64Array}
     */
    get_gasteiger_h_charges() {
        const ret = wasm.jsmol_get_gasteiger_h_charges(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
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
     * `rdMolHash.MolHash(mol, func, useCXSmiles, cxFlagsToSkip)`.
     * @param {string} func
     * @param {boolean} use_cx_smiles
     * @param {number} cx_flags_to_skip
     * @returns {string}
     */
    get_molhash(func, use_cx_smiles, cx_flags_to_skip) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(func, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.jsmol_get_molhash(this.__wbg_ptr, ptr0, len0, use_cx_smiles, cx_flags_to_skip);
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
     * `get_rdkit_descriptors()`: `Descriptors.CalcMolDescriptors(mol)` values in
     * [`rdkit_descriptor_names`] order (NaN where RDKit has none).
     * @returns {Float64Array}
     */
    get_rdkit_descriptors() {
        const ret = wasm.jsmol_get_rdkit_descriptors(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
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
 * `MAP4Calculator` (module docs).
 */
export class Map4Calculator {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        Map4CalculatorFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_map4calculator_free(ptr, 0);
    }
    /**
     * `calculate(mol)`.
     * @param {JSMol} mol
     * @returns {Uint32Array}
     */
    calculateMol(mol) {
        _assertClass(mol, JSMol);
        const ret = wasm.map4calculator_calculateMol(this.__wbg_ptr, mol.__wbg_ptr);
        var v1 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * `calculate(Chem.MolFromSmiles(smiles))`.
     * @param {string} smiles
     * @returns {Uint32Array}
     */
    calculate(smiles) {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.map4calculator_calculate(this.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v2;
    }
    /**
     * @param {string | null} [options]
     */
    constructor(options) {
        var ptr0 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.map4calculator_new(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        Map4CalculatorFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
}
if (Symbol.dispose) Map4Calculator.prototype[Symbol.dispose] = Map4Calculator.prototype.free;

/**
 * `MHFPEncoder` (module docs).
 */
export class MhfpEncoder {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MhfpEncoderFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_mhfpencoder_free(ptr, 0);
    }
    /**
     * `encode_mol(in_mol, …)`.
     * @param {JSMol} mol
     * @param {string | null} [options]
     * @returns {Uint32Array}
     */
    encodeMol(mol, options) {
        _assertClass(mol, JSMol);
        var ptr0 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.mhfpencoder_encodeMol(this.__wbg_ptr, mol.__wbg_ptr, ptr0, len0);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v2;
    }
    /**
     * `encode(in_smiles, radius, rings, kekulize, min_radius, sanitize=True)`.
     * @param {string} smiles
     * @param {string | null} [options]
     * @returns {Uint32Array}
     */
    encode(smiles, options) {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        var ptr1 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len1 = WASM_VECTOR_LEN;
        const ret = wasm.mhfpencoder_encode(this.__wbg_ptr, ptr0, len0, ptr1, len1);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v3 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v3;
    }
    /**
     * `{"nPermutations": 2048, "seed": 42}`.
     * @param {string | null} [options]
     */
    constructor(options) {
        var ptr0 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        var len0 = WASM_VECTOR_LEN;
        const ret = wasm.mhfpencoder_new(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        MhfpEncoderFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
}
if (Symbol.dispose) MhfpEncoder.prototype[Symbol.dispose] = MhfpEncoder.prototype.free;

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
 * [`info`] as JSON text: `tools/parity/wasm-modules/modules.py build` writes the loader's
 * manifest.json (function → module) from the modules it built.
 * @returns {string}
 */
export function cruxModuleInfo() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.cruxModuleInfo();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * The modules this wasm was built with, comma-separated (`"ml"`).
 * @returns {string}
 */
export function cruxModules() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.cruxModules();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * The whole file on this thread → an Arrow IPC stream.
 * @param {Uint8Array} bytes
 * @param {any} options
 * @returns {Uint8Array}
 */
export function csvToArrow(bytes, options) {
    const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.csvToArrow(ptr0, len0, options);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v2 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v2;
}

/**
 * `descriptorNames()`: every name `descriptorsTable` accepts, in order:
 * `[{name, family, type: "double" | "bigint", needs3d, doc}]` — for pickers.
 * @returns {any}
 */
export function descriptorNames() {
    const ret = wasm.descriptorNames();
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * `descriptors3dNames()`: `Descriptors3D.descList` names, in order.
 * @returns {string[]}
 */
export function descriptors3dNames() {
    const ret = wasm.descriptors3dNames();
    var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
    wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
    return v1;
}

/**
 * `descriptorsArrow(smiles, names, seed?)`: [`descriptors_table`] as Arrow IPC
 * stream bytes (`application/vnd.apache.arrow.stream`): the requested columns
 * (Float64; Int64 for counts; null where RDKit has none), then `error`.
 * @param {string[]} smiles
 * @param {string[]} names
 * @param {number | null} [seed]
 * @returns {Uint8Array}
 */
export function descriptorsArrow(smiles, names, seed) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayJsValueToWasm0(names, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.descriptorsArrow(ptr0, len0, ptr1, len1, isLikeNone(seed) ? Number.MAX_SAFE_INTEGER : (seed) >> 0);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v3 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v3;
}

/**
 * `descriptorsTable(smiles, names, seed?)`: RDKit descriptors by name
 * ([`descriptor_names_info`]; `all`, `@<family>`, `@all`) over a SMILES list —
 * only the requested ones are computed. Throws on an unknown name (with the near
 * matches). 3D names embed `AddHs(mol)` with ETKDGv3 (`seed`, default 42).
 * @param {string[]} smiles
 * @param {string[]} names
 * @param {number | null} [seed]
 * @returns {DescriptorTable}
 */
export function descriptorsTable(smiles, names, seed) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArrayJsValueToWasm0(names, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.descriptorsTable(ptr0, len0, ptr1, len1, isLikeNone(seed) ? Number.MAX_SAFE_INTEGER : (seed) >> 0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return DescriptorTable.__wrap(ret[0]);
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
 * The entries of a catalog selection, in match order (`entry` is the index
 * `structuralAlerts` reports).
 * @param {string | null} [catalogs]
 * @returns {any}
 */
export function filterCatalogEntries(catalogs) {
    var ptr0 = isLikeNone(catalogs) ? 0 : passStringToWasm0(catalogs, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.filterCatalogEntries(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Register resource `name` as a dataset read in place: `manifestJson` (its
 * `crux-dataset.json`) and `base`: a base URL (`<base>/<path>` per file, HTTP range), or an
 * object mapping each file's relative path to a source (an OPFS sync access handle, a
 * `File` / `Blob`, bytes, a URL, or a synchronous `{size, read(offset, length)}` reader, e.g.
 * node's `fs.readSync`). Reads block: call functions that read it in a Web Worker.
 * Returns the dataset id (`sha256:` of the canonical manifest), its cache identity.
 * @param {string} name
 * @param {string} manifest_json
 * @param {any} base
 * @returns {string}
 */
export function fnDataset(name, manifest_json, base) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(manifest_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.fnDataset(ptr0, len0, ptr1, len1, base);
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
 * One function's spec (JSON).
 * @param {string} name
 * @returns {string}
 */
export function fnDescribe(name) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.fnDescribe(ptr0, len0);
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
 * Every registered function's spec: `{"crux_fn_version", "functions": [...]}` (JSON).
 * @returns {string}
 */
export function fnList() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.fnList();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Register resource `name` (e.g. `model`) as in-memory bytes for later
 * `callFn`s; `id` identifies them in the loaded-object cache (a URL or hash).
 * @param {string} name
 * @param {string} id
 * @param {Uint8Array} bytes
 */
export function fnResource(name, id, bytes) {
    const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    wasm.fnResource(ptr0, len0, ptr1, len1, ptr2, len2);
}

/**
 * Register resource `name` as one object read in place (`src`: a URL; in a worker also an
 * OPFS sync access handle, a `File` / `Blob` or bytes); `id` identifies it in caches.
 * @param {string} name
 * @param {string} id
 * @param {any} src
 */
export function fnSource(name, id, src) {
    const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.fnSource(ptr0, len0, ptr1, len1, src);
    if (ret[1]) {
        throw takeFromExternrefTable0(ret[0]);
    }
}

/**
 * Generic Murcko scaffold of one molecule (all atoms C, all bonds single) —
 * `MolToSmiles(MakeScaffoldGeneric(GetScaffoldForMol(MolFromSmiles(s))))`.
 * @param {string} smiles
 * @returns {string}
 */
export function genericScaffold(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.genericScaffold(ptr0, len0);
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
 * The functional groups of a column of SMILES in one call: an array of the
 * same length holding each molecule's JSON (as [`JSMol::get_efgs_json`]), or
 * `null` where it does not parse.
 * @param {string[]} smiles
 * @returns {any}
 */
export function getEfgsJson(smiles) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.getEfgsJson(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
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
 * The identity key of `smiles` under the named policy (`molstandardize::identity`).
 * @param {string} smiles
 * @param {string} policy
 * @returns {string}
 */
export function identityKey(smiles, policy) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(policy, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.identityKey(ptr0, len0, ptr1, len1);
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
 * [`identity_key`] of a column: `null` where a molecule fails.
 * @param {string[]} smiles
 * @param {string} policy
 * @returns {any}
 */
export function identityKeys(smiles, policy) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(policy, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.identityKeys(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * The identity policies: `[{name, description}]`.
 * @returns {any}
 */
export function identityPolicies() {
    const ret = wasm.identityPolicies();
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Load the NP_Score model (`npscore_publicnp.cxft` bytes). Once per page;
 * later calls are no-ops.
 * @param {Uint8Array} bytes
 */
export function initNpScore(bytes) {
    const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.initNpScore(ptr0, len0);
    if (ret[1]) {
        throw takeFromExternrefTable0(ret[0]);
    }
}

/**
 * Load the SA_Score fragment table (`sascore_fpscores.cxft` bytes). Once per
 * page; later calls are no-ops.
 * @param {Uint8Array} bytes
 */
export function initSaScore(bytes) {
    const ptr0 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.initSaScore(ptr0, len0);
    if (ret[1]) {
        throw takeFromExternrefTable0(ret[0]);
    }
}

/**
 * @returns {number}
 */
export function leafAbi() {
    const ret = wasm.leafAbi();
    return ret >>> 0;
}

/**
 * @param {string} name
 * @param {any} table
 * @param {any} params
 * @param {string | null | undefined} opts_json
 * @param {any} keep
 * @param {number | null} [call_id]
 * @param {Function | null} [on_progress]
 * @param {Function | null} [on_log]
 * @returns {any}
 */
export function leafCall(name, table, params, opts_json, keep, call_id, on_progress, on_log) {
    const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(opts_json) ? 0 : passStringToWasm0(opts_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.leafCall(ptr0, len0, table, params, ptr1, len1, keep, isLikeNone(call_id) ? Number.MAX_SAFE_INTEGER : (call_id) >>> 0, isLikeNone(on_progress) ? 0 : addToExternrefTable0(on_progress), isLikeNone(on_log) ? 0 : addToExternrefTable0(on_log));
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * @param {string} name
 * @param {any} table
 * @param {any} params
 * @param {string | null | undefined} opts_json
 * @param {any} keep
 * @param {number | null} [call_id]
 * @param {Function | null} [on_progress]
 * @param {Function | null} [on_log]
 * @returns {Promise<any>}
 */
export function leafCallAsync(name, table, params, opts_json, keep, call_id, on_progress, on_log) {
    const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(opts_json) ? 0 : passStringToWasm0(opts_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.leafCallAsync(ptr0, len0, table, params, ptr1, len1, keep, isLikeNone(call_id) ? Number.MAX_SAFE_INTEGER : (call_id) >>> 0, isLikeNone(on_progress) ? 0 : addToExternrefTable0(on_progress), isLikeNone(on_log) ? 0 : addToExternrefTable0(on_log));
    return ret;
}

/**
 * @param {number} call_id
 * @returns {boolean}
 */
export function leafCancel(call_id) {
    const ret = wasm.leafCancel(call_id);
    return ret !== 0;
}

/**
 * `{instance, cap}` of the handle table, set from `cfg` first where given (`instance`: the
 * host's number for this leaf instance, before any handle; `cap`: kept bytes).
 * @param {any} cfg
 * @returns {any}
 */
export function leafConfig(cfg) {
    const ret = wasm.leafConfig(cfg);
    return ret;
}

/**
 * @param {number} handle
 * @returns {any}
 */
export function leafDescribe(handle) {
    const ret = wasm.leafDescribe(handle);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * @param {number} handle
 * @returns {boolean}
 */
export function leafFree(handle) {
    const ret = wasm.leafFree(handle);
    return ret !== 0;
}

/**
 * @param {any} spec
 * @returns {any}
 */
export function leafGather(spec) {
    const ret = wasm.leafGather(spec);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * `{abi, leaf, functions, probes}`: the module's functions (specs: its fns.json) and, per
 * variant with a device probe, whether it runs here.
 * @returns {any}
 */
export function leafInfo() {
    const ret = wasm.leafInfo();
    return ret;
}

/**
 * Pins (`on`) or unpins a kept result: a pinned one is never evicted (a run's live data).
 * @param {number} handle
 * @param {boolean} on
 * @returns {boolean}
 */
export function leafPin(handle, on) {
    const ret = wasm.leafPin(handle, on);
    return ret !== 0;
}

/**
 * A call checked without running it: `{fn, want, shape}` (crux_leaf::abi::plan).
 * @param {string} name
 * @param {any} table
 * @param {any} params
 * @param {string | null | undefined} opts_json
 * @param {any} resources
 * @param {boolean | null} [relax_variant]
 * @returns {any}
 */
export function leafPlan(name, table, params, opts_json, resources, relax_variant) {
    const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(opts_json) ? 0 : passStringToWasm0(opts_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.leafPlan(ptr0, len0, table, params, ptr1, len1, resources, isLikeNone(relax_variant) ? 0xFFFFFF : relax_variant ? 1 : 0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * @param {any} table
 * @param {any} keep
 * @returns {any}
 */
export function leafPut(table, keep) {
    const ret = wasm.leafPut(table, keep);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * @param {number} handle
 * @param {boolean} view
 * @param {string | null} [which]
 * @returns {any}
 */
export function leafRead(handle, view, which) {
    var ptr0 = isLikeNone(which) ? 0 : passStringToWasm0(which, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.leafRead(handle, view, ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Ends lease `lease` (its stored resources evictable again); a call naming it ends it too.
 * @param {number} lease
 * @returns {boolean}
 */
export function leafRelease(lease) {
    const ret = wasm.leafRelease(lease);
    return ret !== 0;
}

/**
 * `fnResource` under the ABI's name.
 * @param {string} name
 * @param {string} id
 * @param {Uint8Array} bytes
 */
export function leafResource(name, id, bytes) {
    const ptr0 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
    const len2 = WASM_VECTOR_LEN;
    wasm.leafResource(ptr0, len0, ptr1, len1, ptr2, len2);
}

/**
 * @returns {any}
 */
export function leafStats() {
    const ret = wasm.leafStats();
    return ret;
}

/**
 * Stores resource bytes by identity: calls bind them by `resources: {name: {id}}` only.
 * @param {string} id
 * @param {Uint8Array} bytes
 * @param {number | null} [lease]
 */
export function leafStore(id, bytes, lease) {
    const ptr0 = passStringToWasm0(id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray8ToWasm0(bytes, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.leafStore(ptr0, len0, ptr1, len1, !isLikeNone(lease), isLikeNone(lease) ? 0 : lease);
    if (ret[1]) {
        throw takeFromExternrefTable0(ret[0]);
    }
}

/**
 * Which of the resource identities `ids` are not stored here (evicted, unstored, never stored);
 * the others marked used.
 * @param {any} ids
 * @param {number | null} [lease]
 * @returns {any}
 */
export function leafStored(ids, lease) {
    const ret = wasm.leafStored(ids, !isLikeNone(lease), isLikeNone(lease) ? 0 : lease);
    return ret;
}

/**
 * Drops stored resource bytes `id`.
 * @param {string} id
 * @returns {boolean}
 */
export function leafUnstore(id) {
    const ptr0 = passStringToWasm0(id, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.leafUnstore(ptr0, len0);
    return ret !== 0;
}

/**
 * Load RDKit Contrib's `SubstructureFilter_HitTriaging_wPubChemExamples.csv`
 * (its text, fetched by the caller; nothing is built in). Returns the number
 * of filters; replaces a previously loaded table.
 * @param {string} csv_text
 * @returns {number}
 */
export function loadNibrCatalog(csv_text) {
    const ptr0 = passStringToWasm0(csv_text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.loadNibrCatalog(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return ret[0] >>> 0;
}

/**
 * tmap `Minhash.get_distance(a, b)` (MAP4).
 * @param {Uint32Array} a
 * @param {Uint32Array} b
 * @returns {number}
 */
export function map4Distance(a, b) {
    const ptr0 = passArray32ToWasm0(a, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray32ToWasm0(b, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.map4Distance(ptr0, len0, ptr1, len1);
    return ret;
}

/**
 * `MHFPEncoder.distance(a, b)`.
 * @param {Uint32Array} a
 * @param {Uint32Array} b
 * @returns {number}
 */
export function mhfpDistance(a, b) {
    const ptr0 = passArray32ToWasm0(a, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passArray32ToWasm0(b, wasm.__wbindgen_malloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.mhfpDistance(ptr0, len0, ptr1, len1);
    return ret;
}

/**
 * `rdMolHash.MolHash(Chem.MolFromSmiles(smiles), func, useCXSmiles, cxFlagsToSkip)`;
 * throws on a SMILES that does not read.
 * @param {string} smiles
 * @param {string} func
 * @param {boolean} use_cx_smiles
 * @param {number} cx_flags_to_skip
 * @returns {string}
 */
export function molHash(smiles, func, use_cx_smiles, cx_flags_to_skip) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(func, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.molHash(ptr0, len0, ptr1, len1, use_cx_smiles, cx_flags_to_skip);
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
 * The MolHash function names, RDKit's order.
 * @returns {string[]}
 */
export function molHashNames() {
    const ret = wasm.molHashNames();
    var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
    wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
    return v1;
}

/**
 * [`mol_hash`] of a column: an array of the same length, `null` where a molecule fails.
 * @param {string[]} smiles
 * @param {string} func
 * @param {boolean} use_cx_smiles
 * @param {number} cx_flags_to_skip
 * @returns {any}
 */
export function molHashes(smiles, func, use_cx_smiles, cx_flags_to_skip) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(func, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.molHashes(ptr0, len0, ptr1, len1, use_cx_smiles, cx_flags_to_skip);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Bemis–Murcko scaffold of one molecule — RDKit
 * `MurckoScaffold.MurckoScaffoldSmiles(smiles, includeChirality)`; `""` for an
 * acyclic molecule. `isomeric = false` is RDKit's default (no stereo, no
 * isotopes). Throws on a SMILES that does not parse.
 * @param {string} smiles
 * @param {boolean} isomeric
 * @returns {string}
 */
export function murckoScaffold(smiles, isomeric) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.murckoScaffold(ptr0, len0, isomeric);
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
 * Murcko (or, with `generic`, generic) scaffolds of a column of SMILES in one
 * call: an array of the same length, `null` where a molecule fails.
 * @param {string[]} smiles
 * @param {boolean} generic
 * @param {boolean} isomeric
 * @returns {any}
 */
export function murckoScaffolds(smiles, generic, isomeric) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.murckoScaffolds(ptr0, len0, generic, isomeric);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * RDKit Contrib `assignFilters` for one molecule, plus each matching NIBR
 * filter with its annotations and match atoms. Needs `loadNibrCatalog` first.
 * @param {string} smiles
 * @returns {any}
 */
export function nibrFilters(smiles) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.nibrFilters(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * NP-likeness with confidence — `npscorer.scoreMolWConfidence`:
 * `{score, confidence, rawScore, fragmentsFound, fragmentsTotal}`, or `null`
 * for a molecule without atoms. Throws on a SMILES that does not parse, or
 * before [`init_np_score`].
 * @param {string} smiles
 * @returns {any}
 */
export function npScore(smiles) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.npScore(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * NP scores of a column of SMILES: an array of `{score, confidence, …}` or
 * `null`. Throws before [`init_np_score`].
 * @param {string[]} smiles
 * @returns {any}
 */
export function npScores(smiles) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.npScores(ptr0, len0);
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
 * QED of one molecule: the three weightings, the eight properties and each
 * property's desirability.
 * @param {string} smiles
 * @returns {any}
 */
export function qed(smiles) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.qed(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * QED (mean weights) per SMILES; `NaN` where a SMILES does not parse.
 * @param {string[]} smiles
 * @returns {Float64Array}
 */
export function qedMany(smiles) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.qedMany(ptr0, len0);
    var v2 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
    return v2;
}

/**
 * `rdkitDescriptorNames()`: `Descriptors._descList` names, in order.
 * @returns {string[]}
 */
export function rdkitDescriptorNames() {
    const ret = wasm.rdkitDescriptorNames();
    var v1 = getArrayJsValueFromWasm0(ret[0], ret[1]);
    wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
    return v1;
}

/**
 * `RegistrationHash.GetMolLayers(mol, enable_tautomer_hash_v2)` and `GetMolHash(layers,
 * scheme)` (`ALL_LAYERS`, `STEREO_INSENSITIVE_LAYERS`, `TAUTOMER_INSENSITIVE_LAYERS`) as
 * `{"layers": {"CANONICAL_SMILES": …, …}, "hash": "<sha1>"}`.
 * @param {string} smiles
 * @param {string} scheme
 * @param {boolean} tautomer_v2
 * @returns {any}
 */
export function registrationHash(smiles, scheme, tautomer_v2) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(scheme, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.registrationHash(ptr0, len0, ptr1, len1, tautomer_v2);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * SA score (1 easy … 10 hard) — `sascorer.calculateScore(MolFromSmiles(s))`;
 * `undefined` for a molecule without atoms. Throws on a SMILES that does not
 * parse, or before [`init_sa_score`].
 * @param {string} smiles
 * @returns {number | undefined}
 */
export function saScore(smiles) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.saScore(ptr0, len0);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    return ret[0] === 0 ? undefined : ret[1];
}

/**
 * SA score with its terms: `{score, fragmentScore, complexityScore,
 * densityCorrection, numAtoms, numChiralCenters, numSpiroAtoms,
 * numBridgeheadAtoms, numMacrocycles}`, or `null`.
 * @param {string} smiles
 * @returns {any}
 */
export function saScoreDetails(smiles) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.saScoreDetails(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * SA scores of a column of SMILES: an array of the same length, `null` where a
 * molecule fails or has no atoms. Throws before [`init_sa_score`].
 * @param {string[]} smiles
 * @returns {any}
 */
export function saScores(smiles) {
    const ptr0 = passArrayJsValueToWasm0(smiles, wasm.__wbindgen_malloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.saScores(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * `MHFPEncoder.secfp_from_smiles(in_smiles, length, radius, rings, kekulize, sanitize=False)`:
 * the on bits.
 * @param {string} smiles
 * @param {string | null} [options]
 * @returns {Uint32Array}
 */
export function secfp(smiles, options) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.secfp(ptr0, len0, ptr1, len1);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v3 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
    return v3;
}

/**
 * `MHFPEncoder.secfp_from_mol(in_mol, …)`: the on bits.
 * @param {JSMol} mol
 * @param {string | null} [options]
 * @returns {Uint32Array}
 */
export function secfpMol(mol, options) {
    _assertClass(mol, JSMol);
    var ptr0 = isLikeNone(options) ? 0 : passStringToWasm0(options, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len0 = WASM_VECTOR_LEN;
    const ret = wasm.secfpMol(mol.__wbg_ptr, ptr0, len0);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v2 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
    return v2;
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
 * Install a panic hook that forwards Rust panics to `console.error`. Called
 * once at module init on wasm; a no-op on native targets.
 */
export function start() {
    wasm.start();
}

/**
 * Every matching FilterCatalog entry for one molecule, with its matches.
 * @param {string} smiles
 * @param {string | null} [catalogs]
 * @returns {any}
 */
export function structuralAlerts(smiles, catalogs) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    var ptr1 = isLikeNone(catalogs) ? 0 : passStringToWasm0(catalogs, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    var len1 = WASM_VECTOR_LEN;
    const ret = wasm.structuralAlerts(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
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
        __wbg_Function_2c2311dcf2911062: function() { return handleError(function (arg0, arg1, arg2, arg3) {
            const ret = Function(getStringFromWasm0(arg0, arg1), getStringFromWasm0(arg2, arg3));
            return ret;
        }, arguments); },
        __wbg_Number_14af1003b8dd5ead: function(arg0) {
            const ret = Number(arg0);
            return ret;
        },
        __wbg_String_8564e559799eccda: function(arg0, arg1) {
            const ret = String(arg1);
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_bigint_get_as_i64_a2383202b9353e4c: function(arg0, arg1) {
            const v = arg1;
            const ret = typeof(v) === 'bigint' ? v : undefined;
            getDataViewMemory0().setBigInt64(arg0 + 8 * 1, isLikeNone(ret) ? BigInt(0) : ret, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, !isLikeNone(ret), true);
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
        __wbg___wbindgen_in_92f62ee1427d9e49: function(arg0, arg1) {
            const ret = arg0 in arg1;
            return ret;
        },
        __wbg___wbindgen_is_bigint_b123553bed3bb382: function(arg0) {
            const ret = typeof(arg0) === 'bigint';
            return ret;
        },
        __wbg___wbindgen_is_falsy_16bd49b68658263e: function(arg0) {
            const ret = !arg0;
            return ret;
        },
        __wbg___wbindgen_is_function_1f9d30630b8b1d3d: function(arg0) {
            const ret = typeof(arg0) === 'function';
            return ret;
        },
        __wbg___wbindgen_is_null_e343b7d08827ba72: function(arg0) {
            const ret = arg0 === null;
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
        __wbg___wbindgen_jsval_eq_02babf21faa37971: function(arg0, arg1) {
            const ret = arg0 === arg1;
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
        __wbg__wbg_cb_unref_dcc1a90847f04c41: function(arg0) {
            arg0._wbg_cb_unref();
        },
        __wbg_apply_a910804df6e1e433: function() { return handleError(function (arg0, arg1, arg2) {
            const ret = arg0.apply(arg1, arg2);
            return ret;
        }, arguments); },
        __wbg_call_1875a20c43a36133: function() { return handleError(function (arg0, arg1, arg2, arg3) {
            const ret = arg0.call(arg1, arg2, arg3);
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
        __wbg_entries_fb6397112b1de25f: function(arg0) {
            const ret = Object.entries(arg0);
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
        __wbg_getResponseHeader_9147b2579988a88a: function() { return handleError(function (arg0, arg1, arg2, arg3) {
            const ret = arg1.getResponseHeader(getStringFromWasm0(arg2, arg3));
            var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        }, arguments); },
        __wbg_getSize_b3c9adeccff7e656: function() { return handleError(function (arg0) {
            const ret = arg0.getSize();
            return ret;
        }, arguments); },
        __wbg_get_31af05bd4842a84f: function() { return handleError(function (arg0, arg1) {
            const ret = Reflect.get(arg0, arg1);
            return ret;
        }, arguments); },
        __wbg_get_658f6698067d9515: function() { return handleError(function (arg0, arg1) {
            const ret = Reflect.get(arg0, arg1);
            return ret;
        }, arguments); },
        __wbg_get_6c896e0571ddae51: function(arg0, arg1) {
            const ret = arg0[arg1 >>> 0];
            return ret;
        },
        __wbg_get_unchecked_288889d017702237: function(arg0, arg1) {
            const ret = arg0[arg1 >>> 0];
            return ret;
        },
        __wbg_get_with_ref_key_6412cf3094599694: function(arg0, arg1) {
            const ret = arg0[arg1];
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
        __wbg_instanceof_Blob_7c6191a4cd5834d3: function(arg0) {
            let result;
            try {
                result = arg0 instanceof Blob;
            } catch (_) {
                result = false;
            }
            const ret = result;
            return ret;
        },
        __wbg_instanceof_FileSystemSyncAccessHandle_b51afaa1de170464: function(arg0) {
            let result;
            try {
                result = arg0 instanceof FileSystemSyncAccessHandle;
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
        __wbg_instanceof_WorkerGlobalScope_8966e756dd00e8b7: function(arg0) {
            let result;
            try {
                result = arg0 instanceof WorkerGlobalScope;
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
        __wbg_isSafeInteger_717808ad6a54bd9e: function(arg0) {
            const ret = Number.isSafeInteger(arg0);
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
        __wbg_new_9ffe4dece67ce9fc: function() { return handleError(function () {
            const ret = new FileReaderSync();
            return ret;
        }, arguments); },
        __wbg_new_d1d1f48568532d74: function() { return handleError(function () {
            const ret = new XMLHttpRequest();
            return ret;
        }, arguments); },
        __wbg_new_ee2291f50781bf1d: function() {
            const ret = new Array();
            return ret;
        },
        __wbg_new_from_slice_37a3c3a194266cb8: function(arg0, arg1) {
            const ret = new Uint32Array(getArrayU32FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_new_from_slice_9a868026ffa4208a: function(arg0, arg1) {
            const ret = new Uint8Array(getArrayU8FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_new_from_slice_f545fd22ddc142b8: function(arg0, arg1) {
            const ret = new Float64Array(getArrayF64FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_new_typed_b01cb72a8af741a3: function(arg0, arg1) {
            try {
                var state0 = {a: arg0, b: arg1};
                var cb0 = (arg0, arg1) => {
                    const a = state0.a;
                    state0.a = 0;
                    try {
                        return wasm_bindgen__convert__closures_____invoke__h274c5cb09364b90a(a, state0.b, arg0, arg1);
                    } finally {
                        state0.a = a;
                    }
                };
                const ret = new Promise(cb0);
                return ret;
            } finally {
                state0.a = 0;
            }
        },
        __wbg_new_with_length_3da0ad195f6f63ba: function(arg0) {
            const ret = new Uint8Array(arg0 >>> 0);
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
        __wbg_of_d5dbecdc59e02e34: function(arg0) {
            const ret = Array.of(arg0);
            return ret;
        },
        __wbg_open_f0c8628338345808: function() { return handleError(function (arg0, arg1, arg2, arg3, arg4, arg5) {
            arg0.open(getStringFromWasm0(arg1, arg2), getStringFromWasm0(arg3, arg4), arg5 !== 0);
        }, arguments); },
        __wbg_overrideMimeType_3ac5d01ed44cd7fd: function() { return handleError(function (arg0, arg1, arg2) {
            arg0.overrideMimeType(getStringFromWasm0(arg1, arg2));
        }, arguments); },
        __wbg_parse_0fc53dead14b3b42: function() { return handleError(function (arg0, arg1) {
            const ret = JSON.parse(getStringFromWasm0(arg0, arg1));
            return ret;
        }, arguments); },
        __wbg_prototypesetcall_bc27214492979395: function(arg0, arg1, arg2) {
            Uint8Array.prototype.set.call(getArrayU8FromWasm0(arg0, arg1), arg2);
        },
        __wbg_push_2baf45db356cf468: function(arg0, arg1) {
            const ret = arg0.push(arg1);
            return ret;
        },
        __wbg_queueMicrotask_9833f9a49df95a49: function(arg0) {
            const ret = arg0.queueMicrotask;
            return ret;
        },
        __wbg_queueMicrotask_a72f977e97f23c5f: function(arg0) {
            queueMicrotask(arg0);
        },
        __wbg_random_5a4cafd2f02395ff: function() {
            const ret = Math.random();
            return ret;
        },
        __wbg_readAsArrayBuffer_14bc98cc9cdf7c0b: function() { return handleError(function (arg0, arg1) {
            const ret = arg0.readAsArrayBuffer(arg1);
            return ret;
        }, arguments); },
        __wbg_read_fcd52999af006e79: function() { return handleError(function (arg0, arg1, arg2, arg3) {
            const ret = arg0.read(getArrayU8FromWasm0(arg1, arg2), arg3);
            return ret;
        }, arguments); },
        __wbg_reject_1548d2ad77237dc9: function(arg0) {
            const ret = Promise.reject(arg0);
            return ret;
        },
        __wbg_resolve_0076e10020304ede: function(arg0) {
            const ret = Promise.resolve(arg0);
            return ret;
        },
        __wbg_responseText_cc18b37f8501c24d: function() { return handleError(function (arg0, arg1) {
            const ret = arg1.responseText;
            var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        }, arguments); },
        __wbg_response_ce3feafb4716c448: function() { return handleError(function (arg0) {
            const ret = arg0.response;
            return ret;
        }, arguments); },
        __wbg_send_9d00bd9d7a2c90b5: function() { return handleError(function (arg0) {
            arg0.send();
        }, arguments); },
        __wbg_setRequestHeader_733b74e99bf9b791: function() { return handleError(function (arg0, arg1, arg2, arg3, arg4) {
            arg0.setRequestHeader(getStringFromWasm0(arg1, arg2), getStringFromWasm0(arg3, arg4));
        }, arguments); },
        __wbg_set_145a351398b48c65: function() { return handleError(function (arg0, arg1, arg2) {
            const ret = Reflect.set(arg0, arg1, arg2);
            return ret;
        }, arguments); },
        __wbg_set_34d08fd992d43d61: function(arg0, arg1, arg2) {
            arg0.set(arg1, arg2 >>> 0);
        },
        __wbg_set_6ae97e73113c4f0b: function(arg0, arg1, arg2) {
            const ret = arg0.set(arg1, arg2);
            return ret;
        },
        __wbg_set_6be42768c690e380: function(arg0, arg1, arg2) {
            arg0[arg1] = arg2;
        },
        __wbg_set_at_a00f11a442de8f74: function(arg0, arg1) {
            arg0.at = arg1;
        },
        __wbg_set_bea140a88be9b277: function(arg0, arg1, arg2) {
            arg0[arg1 >>> 0] = arg2;
        },
        __wbg_set_responseType_1f631c1316f4c16e: function(arg0, arg1) {
            arg0.responseType = __wbindgen_enum_XmlHttpRequestResponseType[arg1];
        },
        __wbg_size_8f1c0b1a1fbb3810: function(arg0) {
            const ret = arg0.size;
            return ret;
        },
        __wbg_slice_afd8274aed4149c0: function() { return handleError(function (arg0, arg1, arg2) {
            const ret = arg0.slice(arg1, arg2);
            return ret;
        }, arguments); },
        __wbg_stack_3b0d974bbf31e44f: function(arg0, arg1) {
            const ret = arg1.stack;
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg_static_accessor_GLOBAL_266715b9d96ba635: function() {
            const ret = typeof global === 'undefined' ? null : global;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_static_accessor_GLOBAL_THIS_10fb7dc1ae063179: function() {
            const ret = typeof globalThis === 'undefined' ? null : globalThis;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_static_accessor_SELF_0b583911f537483a: function() {
            const ret = typeof self === 'undefined' ? null : self;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_static_accessor_WINDOW_d7f903d1508cbdc4: function() {
            const ret = typeof window === 'undefined' ? null : window;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_status_734f92f48f0f8b78: function() { return handleError(function (arg0) {
            const ret = arg0.status;
            return ret;
        }, arguments); },
        __wbg_stringify_52ff602c1cc4fbb6: function() { return handleError(function (arg0) {
            const ret = JSON.stringify(arg0);
            return ret;
        }, arguments); },
        __wbg_then_e71170d78fcf8954: function(arg0, arg1) {
            const ret = arg0.then(arg1);
            return ret;
        },
        __wbg_value_f3c585ee8f5ba40c: function(arg0) {
            const ret = arg0.value;
            return ret;
        },
        __wbindgen_generic_0000000000000001: function(arg0, arg1) {
            // Cast intrinsic for `Closure(Closure { owned: true, function: Function { arguments: [Externref], shim_idx: 1550, ret: Result(Unit), inner_ret: Some(Result(Unit)) }, mutable: true }) -> Externref`.
            const ret = makeMutClosure(arg0, arg1, wasm_bindgen__convert__closures_____invoke__h86605f80b5250a21);
            return ret;
        },
        __wbindgen_generic_0000000000000002: function(arg0) {
            // Cast intrinsic for `F64 -> Externref`.
            const ret = arg0;
            return ret;
        },
        __wbindgen_generic_0000000000000003: function(arg0) {
            // Cast intrinsic for `I64 -> Externref`.
            const ret = arg0;
            return ret;
        },
        __wbindgen_generic_0000000000000004: function(arg0, arg1) {
            // Cast intrinsic for `Ref(Slice(U8)) -> NamedExternref("Uint8Array")`.
            const ret = getArrayU8FromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_generic_0000000000000005: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_generic_0000000000000006: function(arg0) {
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
        "./crux_wasm_bg.js": import0,
    };
}

function wasm_bindgen__convert__closures_____invoke__h86605f80b5250a21(arg0, arg1, arg2) {
    const ret = wasm.wasm_bindgen__convert__closures_____invoke__h86605f80b5250a21(arg0, arg1, arg2);
    if (ret[1]) {
        throw takeFromExternrefTable0(ret[0]);
    }
}

function wasm_bindgen__convert__closures_____invoke__h274c5cb09364b90a(arg0, arg1, arg2, arg3) {
    wasm.wasm_bindgen__convert__closures_____invoke__h274c5cb09364b90a(arg0, arg1, arg2, arg3);
}


const __wbindgen_enum_XmlHttpRequestResponseType = ["", "arraybuffer", "blob", "document", "json", "text"];
const CsvInputFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_csvinput_free(ptr, 1));
const DescriptorTableFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_descriptortable_free(ptr, 1));
const JSMolFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_jsmol_free(ptr, 1));
const Map4CalculatorFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_map4calculator_free(ptr, 1));
const MhfpEncoderFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mhfpencoder_free(ptr, 1));
const MolListFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mollist_free(ptr, 1));

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

const CLOSURE_DTORS = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(state => wasm.__wbindgen_destroy_closure(state.a, state.b));

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

function makeMutClosure(arg0, arg1, f) {
    const state = { a: arg0, b: arg1, cnt: 1 };
    const real = (...args) => {

        // First up with a closure we increment the internal reference
        // count. This ensures that the Rust closure environment won't
        // be deallocated while we're invoking it.
        state.cnt++;
        const a = state.a;
        state.a = 0;
        try {
            return f(a, state.b, ...args);
        } finally {
            state.a = a;
            real._wbg_cb_unref();
        }
    };
    real._wbg_cb_unref = () => {
        if (--state.cnt === 0) {
            wasm.__wbindgen_destroy_closure(state.a, state.b);
            state.a = 0;
            CLOSURE_DTORS.unregister(state);
        }
    };
    CLOSURE_DTORS.register(real, state, state);
    return real;
}

function passArray32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getUint32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
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
