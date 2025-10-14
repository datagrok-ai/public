import {Property} from "./entities";
import {FLOAT_NULL, TYPE, TYPES_SCALAR} from "./const";
import dayjs from "dayjs";
import { TypedEventArgs } from "./widgets";
let api: any = (typeof window !== 'undefined' ? window : global.window);

/** Converts list of Dart objects to JavaScript objects by calling {@link toJs}
 * @param {object[]} params
 * @returns {object[]} - list of JavaScript objects */
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
 * Instantiates the corresponding JS handler for the Dart object [dart]. See also {@link toDart}
 * @param dart - Dart handle
 * @param {boolean} check - when true, throws an exception if the object can't be converted to JS.
 * @returns JavaScript wrapper for the Dart object
 * */
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

/** Extracts a Dart handle from the JavaScript wrapper. See also {@link toJs} */
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


