import {Property} from "./entities";
import {FLOAT_NULL, TYPE, TYPES_SCALAR} from "./const";
import {TypedEventArgs} from "./viewer";
import dayjs from "dayjs";

/** Converts list of Dart objects to JavaScript objects by calling {@link toJs}
 * @param {object[]} params
 * @returns {object[]} - list of JavaScript objects */
export function paramsToJs(params: any): any {
  let result = <any>[];
  for (let i = 0; i < params.length; i++) {
    let type = (<any>window).grok_GetType(params[i]);
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
  let type = (<any>window).grok_GetType(dart);
  if (dart === FLOAT_NULL)
    return null;
  else if (type === TYPE.MAP) {
    let wrapper = (<any>window).grok_GetWrapper(dart);
    for (let key in wrapper) {
      if (wrapper.hasOwnProperty(key)) {
        let type = (<any>window).grok_GetType(wrapper[key]);
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
  }

  let wrapper = (<any>window).grok_GetWrapper(dart);
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
  // if (x instanceof dayjs.Dayjs)
  //   return x.valueOf();
  if (typeof x.toDart === 'function')
    return x.toDart();
  if (typeof x.dart !== 'undefined')
    return x.dart;
  if (isPlainObject(x))
    return (<any>window).grok_JS_To_Map(x);
  return x;
}


