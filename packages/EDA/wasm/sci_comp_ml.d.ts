/* tslint:disable */
/* eslint-disable */

/**
 * Stateful Elastic Net wrapper for the JS boundary.
 *
 * The EDA package uses this with `l1 = l2 = 0` (plain OLS, matching the
 * retired C++ `fitLinearRegressionParamsWithDataNormalizing`), but the
 * regularised path is available too. Columns cross the boundary as one
 * flat `Float64Array` (`[col0, col1, …]`) plus `n_rows`.
 */
export class WasmElasticNet {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Coefficients followed by the bias as the last element — a single
     * `Float64Array` of length `n_features + 1`. With μ/σ set, the
     * values are in the original feature space. This is exactly the
     * `[w₀…w_{m-1}, b]` layout the EDA `_fitLinearRegressionParams…`
     * adapter expects.
     */
    coefficientsWithBias(): Float64Array;
    /**
     * Train on standardised columns (`flat_columns` = concatenation,
     * `n_rows` rows each) against response `y`.
     */
    fit(flat_columns: Float64Array, n_rows: number, y: Float64Array): void;
    /**
     * Construct with the most-used hyperparameters; the rest take their
     * spec defaults (full-batch, `seed = 42`, `fit_intercept = true`).
     *
     * `tol` is the early-stopping threshold on |Δloss|. The EDA OLS path
     * passes a very small `tol` (with a large `epochs`) so full-batch GD
     * converges to the closed-form least-squares solution the retired
     * C++ produced; the interactive paths can pass a looser `tol`.
     */
    constructor(lr: number, epochs: number, l1: number, l2: number, tol: number);
    /**
     * Predict on standardised columns; returns a `Float64Array`.
     */
    predict(flat_columns: Float64Array, n_rows: number): Float64Array;
    /**
     * Optional μ/σ (one per feature) for unwinding `coefficients()` into
     * the original space. Does not affect `fit`/`predict`.
     */
    setFeatureStats(means: Float64Array, stds: Float64Array): void;
}

/**
 * Stateful PCA wrapper. The EDA `_principalComponentAnalysisNipals…`
 * adapter calls [`fit`](WasmPca::fit) then reads
 * [`scores`](WasmPca::scores) to build the `k` score columns.
 */
export class WasmPca {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Fraction of variance explained per component.
     */
    explainedVarianceRatio(): Float64Array;
    /**
     * Fit on already-centred columns (`flat_columns` concatenation,
     * `n_rows` rows each).
     */
    fit(flat_columns: Float64Array, n_rows: number): void;
    /**
     * Loadings `P`, flat `k × m` row-major.
     */
    loadings(): Float64Array;
    /**
     * Number of components actually extracted (≤ requested).
     */
    nComponents(): number;
    /**
     * `n_components = 0` means `min(n_rows−1, m)`. `tol`/`max_iter` take
     * the NIPALS convergence controls; `reorthogonalize`/`strict` use
     * their spec defaults (`false`).
     */
    constructor(n_components: number, tol: number, max_iter: number);
    /**
     * Cached training scores, flat `k × n_rows` row-major (component by
     * component). The EDA adapter slices this into `k` score columns of
     * length `n_rows`.
     */
    scores(): Float64Array;
    /**
     * Optional μ/σ for unwinding loadings (only σ used).
     */
    setFeatureStats(means: Float64Array, stds: Float64Array): void;
    /**
     * Project new centred columns; flat `k × n_rows` row-major.
     */
    transform(flat_columns: Float64Array, n_rows: number): Float64Array;
}

/**
 * Stateful PLS1 wrapper. The EDA `_partialLeastSquareRegression…`
 * adapter assembles the `WASM_OUTPUT_IDX` array (prediction,
 * regr_coeffs, T, U, P, q) from these getters; `regression_coefficients`
 * are already in raw space when stats are set.
 */
export class WasmPls {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Intercept `b₀`.
     */
    bias(): number;
    /**
     * Cumulative explained variance of `y`, length `A`.
     */
    explainedVariance(): Float64Array;
    /**
     * Fit on standardised predictors (`flat_x` concatenation) and
     * standardised response `y`.
     */
    fit(flat_x: Float64Array, n_rows: number, y: Float64Array): void;
    /**
     * `components` is the number of latent factors `A` (1..=min(m, n)).
     */
    constructor(components: number);
    /**
     * Full prediction `b₀ + Σ b_origⱼ·xⱼ` on raw predictors.
     */
    predict(flat_x: Float64Array, n_rows: number): Float64Array;
    /**
     * Regression coefficients (raw space with stats set).
     */
    regressionCoefficients(): Float64Array;
    /**
     * Predictor μ/σ (one per feature) plus response μ/σ, for unwinding
     * the regression coefficients and intercept.
     */
    setFeatureStats(x_means: Float64Array, x_stds: Float64Array, y_mean: number, y_std: number): void;
    /**
     * X-scores `T`, flat `A × n_rows`.
     */
    tScores(): Float64Array;
    /**
     * Y-scores `U`, flat `A × n_rows`.
     */
    uScores(): Float64Array;
    /**
     * X-loadings `P`, flat `A × m`.
     */
    xLoadings(): Float64Array;
    /**
     * Y-loadings `q`, length `A`.
     */
    yLoadings(): Float64Array;
}

/**
 * Stateful softmax-classifier wrapper. The EDA `_fitSoftmax` adapter
 * uses [`params`](WasmSoftmax::params) (`c × (n+1)`, `[W|B]`).
 */
export class WasmSoftmax {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Train on standardised columns (`flat_x` concatenation, `n_rows`
     * samples each) with integer labels `0..c-1` (an `Int32Array`).
     */
    fit(flat_x: Float64Array, n_rows: number, targets: Int32Array): void;
    /**
     * Cross-entropy loss per iteration.
     */
    lossHistory(): Float64Array;
    /**
     * Construct with the five hyperparameters the EDA UI exposes;
     * `class_weighting` (reference parity) and `seed` take their spec
     * defaults (`true`, `42`).
     */
    constructor(classes: number, learning_rate: number, iterations: number, penalty: number, tolerance: number);
    /**
     * Model parameters `[W | B]`, flat `c × (n+1)` row-major. The EDA
     * adapter slices this into `c` columns of length `n+1`.
     */
    params(): Float64Array;
    /**
     * Predicted class indices as a `Uint32Array`.
     */
    predict(flat_x: Float64Array, n_rows: number): Uint32Array;
    /**
     * Optional μ/σ for normalising raw `predict` input. Does not affect
     * training.
     */
    setFeatureStats(means: Float64Array, stds: Float64Array): void;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_wasmelasticnet_free: (a: number, b: number) => void;
    readonly __wbg_wasmpca_free: (a: number, b: number) => void;
    readonly __wbg_wasmpls_free: (a: number, b: number) => void;
    readonly __wbg_wasmsoftmax_free: (a: number, b: number) => void;
    readonly wasmelasticnet_coefficientsWithBias: (a: number) => [number, number];
    readonly wasmelasticnet_fit: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly wasmelasticnet_new: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly wasmelasticnet_predict: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly wasmelasticnet_setFeatureStats: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly wasmpca_explainedVarianceRatio: (a: number) => [number, number];
    readonly wasmpca_fit: (a: number, b: number, c: number, d: number) => [number, number];
    readonly wasmpca_loadings: (a: number) => [number, number, number, number];
    readonly wasmpca_nComponents: (a: number) => number;
    readonly wasmpca_new: (a: number, b: number, c: number) => [number, number, number];
    readonly wasmpca_scores: (a: number) => [number, number, number, number];
    readonly wasmpca_setFeatureStats: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly wasmpca_transform: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly wasmpls_bias: (a: number) => number;
    readonly wasmpls_explainedVariance: (a: number) => [number, number];
    readonly wasmpls_fit: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly wasmpls_new: (a: number) => [number, number, number];
    readonly wasmpls_predict: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly wasmpls_regressionCoefficients: (a: number) => [number, number];
    readonly wasmpls_setFeatureStats: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number];
    readonly wasmpls_tScores: (a: number) => [number, number, number, number];
    readonly wasmpls_uScores: (a: number) => [number, number, number, number];
    readonly wasmpls_xLoadings: (a: number) => [number, number, number, number];
    readonly wasmpls_yLoadings: (a: number) => [number, number, number, number];
    readonly wasmsoftmax_fit: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly wasmsoftmax_lossHistory: (a: number) => [number, number];
    readonly wasmsoftmax_new: (a: number, b: number, c: number, d: number, e: number) => [number, number, number];
    readonly wasmsoftmax_params: (a: number) => [number, number];
    readonly wasmsoftmax_predict: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly wasmsoftmax_setFeatureStats: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
