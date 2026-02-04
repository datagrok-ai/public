/**
 * Implementation of Dart-JavaScript interop conversion functions.
 *
 * This file contains the actual implementation of toJs() and toDart().
 * The wrappers.ts file delegates to these via the DG global for bootstrapping.
 *
 * @see {@link ./dart-interop.md} for full documentation on the interop system.
 * @module wrappers_impl
 */

import {Property} from "./entities";
import {FLOAT_NULL, TYPE, TYPES_SCALAR} from "./const";
import dayjs from "dayjs";
import { TypedEventArgs } from "./widgets";
let api: any = (typeof window !== 'undefined' ? window : global.window);

/**
 * Converts a list of Dart objects to JavaScript objects.
 * Iterates through the array and calls {@link toJs} on non-scalar types.
 *
 * @param params - Array of Dart objects
 * @returns Array of JavaScript objects/wrappers
 */
export function paramsToJs(params: any): any {
  let result = <any>[];
  for (let i = 0; i < params.length; i++) {
    let type = api.grok_GetType(params[i]);
    if (type !== null && (!TYPES_SCALAR.has(type) || type === TYPE.LIST || type === TYPE.MAP))
      result.push(toJs(params[i]));
    else
      result.push(params[i]);
  }

  return result;
}


/**
 * Converts a Dart object to its JavaScript wrapper or native JS value.
 *
 * Type conversion rules:
 * - `null`/`undefined` -> passed through
 * - `FLOAT_NULL` sentinel -> `null`
 * - `Map` -> plain object with recursively converted values
 * - `List` -> array via {@link paramsToJs}
 * - `DG.TypedEventArgs` -> {@link TypedEventArgs} wrapper
 * - `DateTime` -> `dayjs` object
 * - `ByteArray` -> passed through unchanged
 * - `BigInt` -> JavaScript `BigInt`
 * - Objects with existing wrappers -> returns the wrapper
 * - `Property` -> {@link Property} wrapper
 *
 * @param dart - Dart object handle
 * @param check - When true, throws an exception if the type cannot be converted
 * @returns JavaScript wrapper, native JS value, or the original Dart handle
 * @throws Error if `check` is true and type is not supported
 * @see {@link toDart} for the reverse operation
 * @see {@link ./dart-interop.md} for full documentation
 */
export function toJs(dart: any, check: boolean = false): any {
  if (dart === null)
    return null;
  if (dart == undefined)
    return undefined;
  let type = api.grok_GetType(dart);
  if (dart === FLOAT_NULL)
    return null;
  else if (type === TYPE.MAP) {
    let wrapper = api.grok_GetWrapper(dart);
    for (let key in wrapper) {
      if (wrapper.hasOwnProperty(key)) {
        let type = api.grok_GetType(wrapper[key]);
        if (type !== null && (!TYPES_SCALAR.has(type) || type === TYPE.LIST || type === TYPE.MAP))
          wrapper[key] = toJs(wrapper[key]);
      }
    }
    return wrapper;
  } else if (type === TYPE.LIST) {
    return paramsToJs(dart);
  } else if (type === 'DG.TypedEventArgs') {
    return new TypedEventArgs(dart);
  } else if (type === TYPE.DATE_TIME) {
    return dayjs(dart);
  } else if (type === TYPE.BYTE_ARRAY) {
    return dart;
  }
  else if (type === TYPE.BIG_INT)
    return BigInt(api.grok_BigInt_To_BigIntJs(dart));
  let wrapper = api.grok_GetWrapper(dart);
  if (wrapper != null)
    return wrapper;

  if (type === TYPE.PROPERTY)
    return new Property(dart);
  if (check)
    throw `Not supported type: ${type}`;

  return dart;
}

function isPlainObject(value: any) {
  return value?.constructor === Object;
}

/**
 * Converts a JavaScript object to its Dart representation.
 *
 * Conversion rules (checked in order):
 * 1. `null`/`undefined` -> passed through
 * 2. `dayjs` object -> Dart `DateTime`
 * 3. Object with `toDart()` method -> calls `x.toDart()`
 * 4. Object with `dart` property -> returns `x.dart`
 * 5. Plain object `{}` -> Dart `Map`
 * 6. JavaScript `BigInt` -> Dart `BigInt`
 * 7. Everything else -> passed through unchanged
 *
 * @param x - JavaScript object to convert
 * @returns Dart object handle
 * @see {@link toJs} for the reverse operation
 * @see {@link ./dart-interop.md} for full documentation
 */
export function toDart(x: any): any {
  if (x === undefined || x === null)
    return x;
  if (x instanceof dayjs)
    return api.grok_DayJs_To_DateTime(x.valueOf());
  if (typeof x.toDart === 'function')
    return x.toDart();
  if (typeof x.dart !== 'undefined')
    return x.dart;
  if (isPlainObject(x))
    return api.grok_JS_To_Map(x);
  if (typeof x === 'bigint')
    return api.grok_BigIntJs_To_BigInt(x.toString());
  return x;
}


