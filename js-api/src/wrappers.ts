/**
 * Wrapper functions for Dart-JavaScript interop.
 *
 * This file delegates to `DG.toJs()` and `DG.toDart()` rather than containing
 * the actual implementation. This is intentional for bootstrapping: during initial
 * load, the `DG` global is being set up and needs these simple delegates.
 * The actual implementation is in `wrappers_impl.ts`.
 *
 * @see {@link ./dart-interop.md} for full documentation on the interop system.
 * @module wrappers
 */

declare let DG: any;

/**
 * Converts a list of Dart objects to JavaScript objects.
 * Delegates to DG.paramsToJs() for bootstrapping compatibility.
 * @see wrappers_impl.ts for the actual implementation
 */
export function paramsToJs(params: any): any {
  return DG.paramsToJs(params);
}

/**
 * Converts a Dart object to its JavaScript wrapper.
 * Delegates to DG.toJs() for bootstrapping compatibility.
 *
 * @param dart - Dart object handle
 * @param check - When true, throws if the type cannot be converted
 * @returns JavaScript wrapper or native JS value
 * @see wrappers_impl.ts for the actual implementation
 * @see {@link ./dart-interop.md} for type conversion rules
 */
export function toJs(dart: any, check: boolean = false): any {
  return DG.toJs(dart);
}

/**
 * Extracts or converts a JavaScript object to its Dart representation.
 * Delegates to DG.toDart() for bootstrapping compatibility.
 *
 * @param x - JavaScript object (may have a `.dart` property or `.toDart()` method)
 * @returns Dart object handle
 * @see wrappers_impl.ts for the actual implementation
 * @see {@link ./dart-interop.md} for type conversion rules
 */
export function toDart(x: any): any {
  return DG.toDart(x);
}
