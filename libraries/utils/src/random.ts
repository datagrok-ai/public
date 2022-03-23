/**
 * Generates single random float from 0 to range.
 *
 * @export
 * @param {number} range Max generating value.
 * @return {number} A random float generated.
 */
 export function randomFloat(range: number): number {
    return Math.random() * range;
 }
  
/**
 * Generates single random integer from 0 to range.
 *
 * @export
 * @param {number} range Max generating value.
 * @return {number} A random integer generated.
 */
 export function randomInt(range: number): number {
    return Math.floor(randomFloat(range));
 }
