/**
 * Denotes a vector of floating poit values.
 *
 * @export
 * @class Vector
 * @extends {Float32Array}
 */
export class Vector extends Float32Array {}

/**
 * Denotes a two-dimensional matrix.
 *
 * @export
 * @class Matrix
 * @extends {Array<Vector>}
 */
export class Matrix extends Array<Vector> {}

/**
 * Denotes cartesian coordinates.
 *
 * @export
 * @class Coordinates
 * @extends {Matrix}
 */
export class Coordinates extends Matrix {}

/**
 * Denotes an array of arbitrary-typed vectors.
 *
 * @export
 * @class Vectors
 * @extends {Array<any>}
 */
export class Vectors extends Array<any> {}

/**
 * Denotes a dictionary containing function options.
 *
 * @export
 * @type Options
 */
export type Options = {[name: string]: any};

/**
 * Denotes custom distance metric between the two given vectors.
 *
 * @export
 * @type DistanceMetric
 * @param {any} v1 The first vector.
 * @param {any} v2 The second vector.
 * @return {number} Distance between these two vectors.
 */
 export type DistanceMetric = (v1: any, v2: any) => (number);

/**
 * Denotes a simple string to string dictionary.
 *
 * @export
 * @type Options
 */
 export type StringDictionary = {[key: string]: string};