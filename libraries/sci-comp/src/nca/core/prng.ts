/**
 * Mulberry32 — small, fast, seedable PRNG with a 32-bit state.
 *
 * Period 2³². Quality is sufficient for bootstrap resampling at the
 * iteration counts typical for NCA (≤ 10⁵ resamples). Output is uniform on
 * `[0, 1)` by construction. Not suitable for cryptography.
 *
 * @example
 * const rand = mulberry32(42);
 * rand(); // 0.6011037519201636 — same on every run
 * rand(); // 0.4426399066392332
 *
 * @param seed - 32-bit seed. Coerced via `>>> 0`, so any safe integer works.
 * @returns A function that, on each call, returns the next pseudorandom
 *          float in `[0, 1)`.
 */
export function mulberry32(seed: number): () => number {
  let state = seed >>> 0;
  return () => {
    state = (state + 0x6D2B79F5) >>> 0;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Derive deterministic per-worker seeds from a single master seed.
 *
 * Used by the bootstrap worker pool (lives in `packages/NCA/`) to obtain a
 * sub-stream PRNG for each worker. Implementation (Variant A in the design
 * doc): run `mulberry32(masterSeed)` and convert the first `workerCount`
 * floats in `[0, 1)` to `uint32` seeds.
 *
 * Reproducibility guarantee: at fixed `masterSeed` AND fixed `workerCount`,
 * the returned seeds are identical across runs. **At different
 * `workerCount` the seeds differ** — this is an accepted limitation of
 * Mulberry32 (it is not a splittable PRNG). Pin the worker count when
 * bootstrap reproducibility is required.
 *
 * @param masterSeed - The user-supplied seed.
 * @param workerCount - Number of workers (`≥ 0`). `0` returns an empty array.
 * @returns A `Uint32Array` of length `workerCount` with deterministic seeds.
 */
export function deriveWorkerSeeds(masterSeed: number, workerCount: number): Uint32Array {
  const masterRand = mulberry32(masterSeed);
  const seeds = new Uint32Array(workerCount);
  for (let i = 0; i < workerCount; i++)
    seeds[i] = (masterRand() * 0x100000000) >>> 0;
  return seeds;
}
