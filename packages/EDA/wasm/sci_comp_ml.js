/* @ts-self-types="./sci_comp_ml.d.ts" */

/**
 * Stateful Elastic Net wrapper for the JS boundary.
 *
 * The EDA package uses this with `l1 = l2 = 0` (plain OLS, matching the
 * retired C++ `fitLinearRegressionParamsWithDataNormalizing`), but the
 * regularised path is available too. Columns cross the boundary as one
 * flat `Float64Array` (`[col0, col1, …]`) plus `n_rows`.
 */
export class WasmElasticNet {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmElasticNetFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmelasticnet_free(ptr, 0);
    }
    /**
     * Coefficients followed by the bias as the last element — a single
     * `Float64Array` of length `n_features + 1`. With μ/σ set, the
     * values are in the original feature space. This is exactly the
     * `[w₀…w_{m-1}, b]` layout the EDA `_fitLinearRegressionParams…`
     * adapter expects.
     * @returns {Float64Array}
     */
    coefficientsWithBias() {
        const ret = wasm.wasmelasticnet_coefficientsWithBias(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Train on standardised columns (`flat_columns` = concatenation,
     * `n_rows` rows each) against response `y`.
     * @param {Float64Array} flat_columns
     * @param {number} n_rows
     * @param {Float64Array} y
     */
    fit(flat_columns, n_rows, y) {
        const ptr0 = passArrayF64ToWasm0(flat_columns, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmelasticnet_fit(this.__wbg_ptr, ptr0, len0, n_rows, ptr1, len1);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * Construct with the most-used hyperparameters; the rest take their
     * spec defaults (full-batch, `seed = 42`, `fit_intercept = true`).
     *
     * `tol` is the early-stopping threshold on |Δloss|. The EDA OLS path
     * passes a very small `tol` (with a large `epochs`) so full-batch GD
     * converges to the closed-form least-squares solution the retired
     * C++ produced; the interactive paths can pass a looser `tol`.
     * @param {number} lr
     * @param {number} epochs
     * @param {number} l1
     * @param {number} l2
     * @param {number} tol
     */
    constructor(lr, epochs, l1, l2, tol) {
        const ret = wasm.wasmelasticnet_new(lr, epochs, l1, l2, tol);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        WasmElasticNetFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Predict on standardised columns; returns a `Float64Array`.
     * @param {Float64Array} flat_columns
     * @param {number} n_rows
     * @returns {Float64Array}
     */
    predict(flat_columns, n_rows) {
        const ptr0 = passArrayF64ToWasm0(flat_columns, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmelasticnet_predict(this.__wbg_ptr, ptr0, len0, n_rows);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v2;
    }
    /**
     * Optional μ/σ (one per feature) for unwinding `coefficients()` into
     * the original space. Does not affect `fit`/`predict`.
     * @param {Float64Array} means
     * @param {Float64Array} stds
     */
    setFeatureStats(means, stds) {
        const ptr0 = passArrayF64ToWasm0(means, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(stds, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmelasticnet_setFeatureStats(this.__wbg_ptr, ptr0, len0, ptr1, len1);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
}
if (Symbol.dispose) WasmElasticNet.prototype[Symbol.dispose] = WasmElasticNet.prototype.free;

/**
 * Stateful PCA wrapper. The EDA `_principalComponentAnalysisNipals…`
 * adapter calls [`fit`](WasmPca::fit) then reads
 * [`scores`](WasmPca::scores) to build the `k` score columns.
 */
export class WasmPca {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmPcaFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmpca_free(ptr, 0);
    }
    /**
     * Fraction of variance explained per component.
     * @returns {Float64Array}
     */
    explainedVarianceRatio() {
        const ret = wasm.wasmpca_explainedVarianceRatio(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Fit on already-centred columns (`flat_columns` concatenation,
     * `n_rows` rows each).
     * @param {Float64Array} flat_columns
     * @param {number} n_rows
     */
    fit(flat_columns, n_rows) {
        const ptr0 = passArrayF64ToWasm0(flat_columns, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmpca_fit(this.__wbg_ptr, ptr0, len0, n_rows);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * Loadings `P`, flat `k × m` row-major.
     * @returns {Float64Array}
     */
    loadings() {
        const ret = wasm.wasmpca_loadings(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Number of components actually extracted (≤ requested).
     * @returns {number}
     */
    nComponents() {
        const ret = wasm.wasmpca_nComponents(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * `n_components = 0` means `min(n_rows−1, m)`. `tol`/`max_iter` take
     * the NIPALS convergence controls; `reorthogonalize`/`strict` use
     * their spec defaults (`false`).
     * @param {number} n_components
     * @param {number} tol
     * @param {number} max_iter
     */
    constructor(n_components, tol, max_iter) {
        const ret = wasm.wasmpca_new(n_components, tol, max_iter);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        WasmPcaFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Cached training scores, flat `k × n_rows` row-major (component by
     * component). The EDA adapter slices this into `k` score columns of
     * length `n_rows`.
     * @returns {Float64Array}
     */
    scores() {
        const ret = wasm.wasmpca_scores(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Optional μ/σ for unwinding loadings (only σ used).
     * @param {Float64Array} means
     * @param {Float64Array} stds
     */
    setFeatureStats(means, stds) {
        const ptr0 = passArrayF64ToWasm0(means, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(stds, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmpca_setFeatureStats(this.__wbg_ptr, ptr0, len0, ptr1, len1);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * Project new centred columns; flat `k × n_rows` row-major.
     * @param {Float64Array} flat_columns
     * @param {number} n_rows
     * @returns {Float64Array}
     */
    transform(flat_columns, n_rows) {
        const ptr0 = passArrayF64ToWasm0(flat_columns, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmpca_transform(this.__wbg_ptr, ptr0, len0, n_rows);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v2;
    }
}
if (Symbol.dispose) WasmPca.prototype[Symbol.dispose] = WasmPca.prototype.free;

/**
 * Stateful PLS1 wrapper. The EDA `_partialLeastSquareRegression…`
 * adapter assembles the `WASM_OUTPUT_IDX` array (prediction,
 * regr_coeffs, T, U, P, q) from these getters; `regression_coefficients`
 * are already in raw space when stats are set.
 */
export class WasmPls {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmPlsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmpls_free(ptr, 0);
    }
    /**
     * Intercept `b₀`.
     * @returns {number}
     */
    bias() {
        const ret = wasm.wasmpls_bias(this.__wbg_ptr);
        return ret;
    }
    /**
     * Cumulative explained variance of `y`, length `A`.
     * @returns {Float64Array}
     */
    explainedVariance() {
        const ret = wasm.wasmpls_explainedVariance(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Fit on standardised predictors (`flat_x` concatenation) and
     * standardised response `y`.
     * @param {Float64Array} flat_x
     * @param {number} n_rows
     * @param {Float64Array} y
     */
    fit(flat_x, n_rows, y) {
        const ptr0 = passArrayF64ToWasm0(flat_x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(y, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmpls_fit(this.__wbg_ptr, ptr0, len0, n_rows, ptr1, len1);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * `components` is the number of latent factors `A` (1..=min(m, n)).
     * @param {number} components
     */
    constructor(components) {
        const ret = wasm.wasmpls_new(components);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        WasmPlsFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Full prediction `b₀ + Σ b_origⱼ·xⱼ` on raw predictors.
     * @param {Float64Array} flat_x
     * @param {number} n_rows
     * @returns {Float64Array}
     */
    predict(flat_x, n_rows) {
        const ptr0 = passArrayF64ToWasm0(flat_x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmpls_predict(this.__wbg_ptr, ptr0, len0, n_rows);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v2;
    }
    /**
     * Regression coefficients (raw space with stats set).
     * @returns {Float64Array}
     */
    regressionCoefficients() {
        const ret = wasm.wasmpls_regressionCoefficients(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Predictor μ/σ (one per feature) plus response μ/σ, for unwinding
     * the regression coefficients and intercept.
     * @param {Float64Array} x_means
     * @param {Float64Array} x_stds
     * @param {number} y_mean
     * @param {number} y_std
     */
    setFeatureStats(x_means, x_stds, y_mean, y_std) {
        const ptr0 = passArrayF64ToWasm0(x_means, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(x_stds, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmpls_setFeatureStats(this.__wbg_ptr, ptr0, len0, ptr1, len1, y_mean, y_std);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * X-scores `T`, flat `A × n_rows`.
     * @returns {Float64Array}
     */
    tScores() {
        const ret = wasm.wasmpls_tScores(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Y-scores `U`, flat `A × n_rows`.
     * @returns {Float64Array}
     */
    uScores() {
        const ret = wasm.wasmpls_uScores(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * X-loadings `P`, flat `A × m`.
     * @returns {Float64Array}
     */
    xLoadings() {
        const ret = wasm.wasmpls_xLoadings(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Y-loadings `q`, length `A`.
     * @returns {Float64Array}
     */
    yLoadings() {
        const ret = wasm.wasmpls_yLoadings(this.__wbg_ptr);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
}
if (Symbol.dispose) WasmPls.prototype[Symbol.dispose] = WasmPls.prototype.free;

/**
 * Stateful softmax-classifier wrapper. The EDA `_fitSoftmax` adapter
 * uses [`params`](WasmSoftmax::params) (`c × (n+1)`, `[W|B]`).
 */
export class WasmSoftmax {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        WasmSoftmaxFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_wasmsoftmax_free(ptr, 0);
    }
    /**
     * Train on standardised columns (`flat_x` concatenation, `n_rows`
     * samples each) with integer labels `0..c-1` (an `Int32Array`).
     * @param {Float64Array} flat_x
     * @param {number} n_rows
     * @param {Int32Array} targets
     */
    fit(flat_x, n_rows, targets) {
        const ptr0 = passArrayF64ToWasm0(flat_x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArray32ToWasm0(targets, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmsoftmax_fit(this.__wbg_ptr, ptr0, len0, n_rows, ptr1, len1);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
    /**
     * Cross-entropy loss per iteration.
     * @returns {Float64Array}
     */
    lossHistory() {
        const ret = wasm.wasmsoftmax_lossHistory(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Construct with the five hyperparameters the EDA UI exposes;
     * `class_weighting` (reference parity) and `seed` take their spec
     * defaults (`true`, `42`).
     * @param {number} classes
     * @param {number} learning_rate
     * @param {number} iterations
     * @param {number} penalty
     * @param {number} tolerance
     */
    constructor(classes, learning_rate, iterations, penalty, tolerance) {
        const ret = wasm.wasmsoftmax_new(classes, learning_rate, iterations, penalty, tolerance);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        WasmSoftmaxFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Model parameters `[W | B]`, flat `c × (n+1)` row-major. The EDA
     * adapter slices this into `c` columns of length `n+1`.
     * @returns {Float64Array}
     */
    params() {
        const ret = wasm.wasmsoftmax_params(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Predicted class indices as a `Uint32Array`.
     * @param {Float64Array} flat_x
     * @param {number} n_rows
     * @returns {Uint32Array}
     */
    predict(flat_x, n_rows) {
        const ptr0 = passArrayF64ToWasm0(flat_x, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.wasmsoftmax_predict(this.__wbg_ptr, ptr0, len0, n_rows);
        if (ret[3]) {
            throw takeFromExternrefTable0(ret[2]);
        }
        var v2 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v2;
    }
    /**
     * Optional μ/σ for normalising raw `predict` input. Does not affect
     * training.
     * @param {Float64Array} means
     * @param {Float64Array} stds
     */
    setFeatureStats(means, stds) {
        const ptr0 = passArrayF64ToWasm0(means, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArrayF64ToWasm0(stds, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.wasmsoftmax_setFeatureStats(this.__wbg_ptr, ptr0, len0, ptr1, len1);
        if (ret[1]) {
            throw takeFromExternrefTable0(ret[0]);
        }
    }
}
if (Symbol.dispose) WasmSoftmax.prototype[Symbol.dispose] = WasmSoftmax.prototype.free;
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg_Error_bce6d499ff0a4aff: function(arg0, arg1) {
            const ret = Error(getStringFromWasm0(arg0, arg1));
            return ret;
        },
        __wbg___wbindgen_throw_9c31b086c2b26051: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
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
        "./sci_comp_ml_bg.js": import0,
    };
}

const WasmElasticNetFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmelasticnet_free(ptr, 1));
const WasmPcaFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmpca_free(ptr, 1));
const WasmPlsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmpls_free(ptr, 1));
const WasmSoftmaxFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_wasmsoftmax_free(ptr, 1));

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

function getArrayU32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
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

function passArray32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getUint32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passArrayF64ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 8, 8) >>> 0;
    getFloat64ArrayMemory0().set(arg, ptr / 8);
    WASM_VECTOR_LEN = arg.length;
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

let WASM_VECTOR_LEN = 0;

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedFloat64ArrayMemory0 = null;
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
        module_or_path = new URL('sci_comp_ml_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
