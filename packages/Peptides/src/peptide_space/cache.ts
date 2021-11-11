type CachingDict = {[key: string]: number};

// eslint-disable-next-line no-unused-vars
export class CachedCalculator {
  cache: CachingDict;
  constructor() {
    this.cache = {};
  }

  protected makeKey(...args: any[]): string {
    throw new Error('Method not implemented.');
  }

  /**
   * calcCached
   * @param {Function} func A function to call.
   * @param {any[]} args Arguments to path through.
   * @return {any} Result of calculation.
   */
  public calcCached(func: Function, ...args: any[]): any {
    const key = this.makeKey(...args);
    let result: any;

    if (!(key in this.cache)) {
      result = func(...args);
      this.cache[key] = result;
    }
    return result;
  }

  /**
   * resetCache
   */
  public resetCache() {
    this.cache = {};
  }
}
