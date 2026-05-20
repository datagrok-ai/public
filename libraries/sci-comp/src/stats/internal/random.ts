/**
 * Internal pseudo-random number generation.
 *
 * Used by methods that need reproducible randomness (e.g. permutation tests).
 * The default `Math.random` is non-seedable, so callers that need
 * reproducible output should pass an `rng` from `mulberry32` (or any other
 * `() => number` returning values in [0, 1)).
 */

/**
 * Mulberry32 — a tiny seedable PRNG with a period of 2^32.
 *
 * Output: `() => number` in `[0, 1)`. Suitable for permutation tests and
 * similar simulations where statistical quality is more important than
 * cryptographic strength.
 *
 * @param seed Any 32-bit integer (will be coerced via `>>> 0`).
 * @returns A function that produces uniform `[0, 1)` numbers.
 */
export function mulberry32(seed: number): () => number {
  let s = seed >>> 0;
  return function() {
    s = (s + 0x6D2B79F5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * In-place Fisher-Yates shuffle of a typed array using the supplied PRNG.
 *
 * @param arr Array to shuffle in place.
 * @param rng PRNG returning values in `[0, 1)`. Defaults to `Math.random`.
 */
export function shuffleInPlace(arr: Float64Array, rng: () => number = Math.random): void {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    const tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
  }
}
