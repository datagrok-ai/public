/**
 * Wraps an async function `fn` that fetches data from an API and returns a Promise. 
 * If `fn` succeeds, its result is returned. If it fails, the failure is logged and its error is rethrown.
 *
 * @template T - The type of the Promise returned by `fn`.
 * @param {() => Promise<T>} fn - The async function that fetches data from an API to wrap.
 * @returns {Promise<T>} - Resolves with `fn`'s result or rejects with `fn`'s error.
 */
export async function fetchWrapper<T>(fn: () => Promise<T>): Promise<T> {
  try {
    return await fn();
  } catch (error) {
    console.error('Promise rejected', error);
    throw error;
  }
}