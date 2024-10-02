import { toJs } from "./wrappers";
import { IDartApi } from "./api/grok_api.g";

const api: IDartApi = <any>window;

/**
 * Functionality to store key-value pairs with settings for further reuse.
 * Keys and values are stored in memory but are synchronized with the server every 10 seconds if data was changed.
 * Value can't be longer then 5000 symbols.
 * Use {@param isPrivate} to set what storage to use - private to the current user or shared with every user.
 * */
export class UserSettingsStorage {

  constructor() {
  }

  /** Adds a single value to storage
   * @param {string} name Storage name
   * @param {string} key
   * @param {string} value
   * @param {boolean} isPrivate Value will be available only for current user. If false, shared storage is used.
   * */
  add(name: string, key: string, value: string, isPrivate: boolean = true): void {
    api.grok_UserSettings_Add(name, key, value, isPrivate);
  }

  /**
   * Adds key-value pairs of object to storage
   * @param {string} name Storage name
   * @param {{[key: string]: string}} data
   * @param {boolean} isPrivate Value will be available only for current user. If false, shared storage is used.
   */
  addAll(name: string, data: {[key: string]: string}, isPrivate: boolean = true): void {
    api.grok_UserSettings_AddAll(name, data, isPrivate);
  }

  /**
   * Replaces existing data in storage with {@param data}
   * @param {string} name Storage name
   * @param {{[key: string]: string}} data
   * @param {boolean} isPrivate Data will be replaced in private or shared storage
   */
  put(name: string, data: {[key: string]: string}, isPrivate: boolean = true): void {
    api.grok_UserSettings_Put(name, data, isPrivate);
  }

  /**
   * Retrieves map from storage
   * @param {string} name Storage name
   * @param {boolean} isPrivate Map will be retrieved from private or shared storage.
   */
  get(name: string, isPrivate: boolean = true): {[key: string]: string} | undefined {
    return toJs(api.grok_UserSettings_Get(name, isPrivate));
  }

  /**
   * Retrieves value from storage
   * @param {string} name Storage name
   * @param {string} key
   * @param {boolean} isPrivate Value will be retrieved from private or shared storage.
   */
  getValue(name: string, key: string, isPrivate: boolean = true): string | undefined {
    return api.grok_UserSettings_GetValue(name, key, isPrivate);
  }

  /**
   * Deletes single value from storage
   * @param {string} name Storage name
   * @param {string} key
   * @param {boolean} isPrivate Value will be deleted from private or shared storage.
   */
  delete(name: string, key: string, isPrivate: boolean = true): void {
    return api.grok_UserSettings_Delete(name, key, isPrivate);
  }
}
