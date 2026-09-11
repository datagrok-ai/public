import {_toIterable} from "./utils_convert";
import {toDart, toJs} from "./wrappers";
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** A proxy to a Dart `List<T>`. */
export class DartList<T> implements Iterable<T> {
    dart: any;

    /** Creates a proxy to an existing Dart list. */
    static fromDart<T>(dart: any): DartList<T> {
        let list = new DartList();
        list.dart = dart;
        return list as DartList<T>;
    }

    [key: number]: any;

    /** Returns the number of objects in this list. */
    get length(): number { return api.grok_List_Get_Length(this.dart); }

    /** Removes all objects from this list; the length of the list becomes zero. */
    clear() { api.grok_List_Clear(this.dart); }

    /** Sorts this list. */
    sort() { api.grok_List_Sort(this.dart); }

    /** Adds [value] to the end of this list, extending the length by one. */
    push(value: T) { api.grok_List_Add(this.dart, value); }

    /** Returns the object at the given [index] in the list. */
    get(index: number): T { return api.grok_List_Get(this.dart, index); }

    /** Sets the value at the given [index] in the list to [value]. */
    set(index: number, value: T): T { return api.grok_List_Set(this.dart, index, value); }

    /** Removes the first occurrence of [value] from the list. */
    remove(value: T) { api.grok_List_Remove(this.dart, value); }

    includes(item: T, start?: number) {
        const length = this.length;
        for (let i = (start ? start : 0); i < length; i++) {
            if (this.get(i) === item)
                return true;
        }
        return false;
    }

    * [Symbol.iterator](): Iterator<T> {
        for (let i = 0; i < this.length; i++)
            yield this.get(i);
    }
}


/**
 * Proxies a Dart Map, API-compliant to ES 2015+
 */
/** Map methods of a {@link MapBag}, in addition to indexed access. */
export interface IMapProxy<V = any> {
  keys(): Iterable<string>;
  values(): Iterable<V>;
  entries(): Iterable<[string, V]>;
  forEach(callback: (key: string, value: V) => void): void;
  get(key: string): V;
  set(key: string, value: V): IMapProxy<V>;
  has(key: string): boolean;
  delete(key: string): boolean;
  clear(): void;
  size(): number;
  [Symbol.iterator](): Iterator<[string, V]>;
}

/** A Dart-backed string-keyed bag (`df.tags`, `column.temp`, `call.aux`, ...): read and write entries
 * by index, or use the map methods. Writes go straight to the Dart object. */
export type MapBag<V = any> = {[key: string]: V} & IMapProxy<V>;

/** A proxy to a Dart `Map<String, T>`; see {@link MapBag}. `valueType` restricts writes to one `typeof`. */
export const MapProxy: new <V = any>(dart: any, objectName?: string | null, valueType?: string | null) => MapBag<V> = new Proxy(class {
        dart: any;
        objectName: string | null;
        valueType: string | null;
        constructor(dart: any, objectName: string | null = null, valueType: string | null = null) {
            this.dart = dart;
            this.objectName = objectName;
            this.valueType = valueType;
        }
        keys(): Iterable<any> {
            return _toIterable(api.grok_Map_Keys(this.dart));
        }
        values() {
            return _toIterable(toJs(api.grok_Map_Values(this.dart)));
        }
        * [Symbol.iterator] () {
            for (let key of this.keys()) {
                const value = toJs(api.grok_Map_Get(this.dart, key));
                yield [key, toJs(value)];
            }
        }
        entries() {
            return this;
        }
        forEach(callback: (key: string, value: any) => void) {
            for (const [key, value] of this) {
                callback(key, toJs(value));
            }
        }
        delete(key: string) {
            return toJs(api.grok_Map_Delete(this.dart, toDart(key)));
        }
        get(key: any) {
            return toJs(api.grok_Map_Get(this.dart, key));
        }
        has(key: string) {
            return api.grok_Map_Has(this.dart, toDart(key));
        }
        set(key: string, value: any) {
            api.grok_Map_Set(this.dart, key, toDart(value));
            return this;
        }
        clear() {
            api.grok_Map_Clear(this.dart);
        }
        size() {
            return toJs(api.grok_Map_Size(this.dart));
        }
    }, {
        construct(target, args) {
            // @ts-ignore
            return new Proxy(new target(...args), {
                get: function (target: any, prop) {
                    const val = target[prop];
                    if (typeof val === 'function') {
                        return function (...args :string[]) {
                            return val.apply(target, args);
                        };
                    } else {
                        return toJs(api.grok_Map_Get(target.dart, prop));
                    }
                },
                set: function (target, prop, value) {
                    const valueType = typeof(value);
                    if (!target.valueType || target.valueType === valueType) {
                        api.grok_Map_Set(target.dart, prop, toDart(value));
                    } else {
                        throw new Error(`Entries of ${target.objectName} require type '${target.valueType}', passed '${valueType}'`);
                    }
                    return true;
                },
                deleteProperty: function (target, prop) {
                    api.grok_Map_Delete(target.dart, toDart(prop));
                    return true;
                },
                has: function (target, prop) {
                    return api.grok_Map_Has(target.dart, toDart(prop));
                },
                getOwnPropertyDescriptor(target, prop) {
                    return {
                        enumerable: true,
                        configurable: true
                    };
                },
                ownKeys: function (target) {
                    return Array.from(target.keys());
                }
            });
        }
    }
) as any;

// export class PropProxy {
//     constructor(dart) {
//         this.dart = dart;
//         return new Proxy({}, {
//             get: function(target, prop) { return DG.toJs(api.grok_PropMixin_Get(dart, prop)); },
//             set: function(target, prop, value) { api.grok_PropMixin_Set(dart, prop, DG.toDart(value)); }
//         })
//     }
// }
