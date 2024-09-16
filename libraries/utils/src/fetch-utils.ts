/**
 * Wraps an async function `fn` that fetches data from an API and returns a Promise. 
 * If `fn` succeeds, its result is returned. If it fails, a custom error with the message 
 * 'Promise rejected' is thrown.
 *
 * @template T - The type of the Promise returned by `fn`.
 * @param {() => Promise<T>} fn - The async function that fetches data from an API to wrap.
 * @returns {Promise<T>} - Resolves with `fn`'s result or rejects with 
 *                         'Promise rejected' if `fn` fails.
 * @throws {Error} - Throws 'Promise rejected' on failure.
 */
export async function fetchWrapper<T>(fn: () => Promise<T>): Promise<T> {
  try {
    return await fn();
  } catch (error) {
    console.error('Promise rejected');
    throw new Error('Promise rejected');
  }
}