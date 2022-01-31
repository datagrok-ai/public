export namespace ArrayUtils {

  /**
   * Finds indices of element(s) found with the given filter.
   *
   * @export
   * @template T Type of elements in the array.
   * @param {T[]} values Array to search in.
   * @param {(item: T) => boolean} filter Filtering function including an element if it fits a condition.
   * @return {number[]} Idices of elements corresponding the filter.
   */
  export function indexesOf<T>(values: T[], filter: (item: T) => boolean): number[] {
    const indexes = [];
    for (let i = 0; i < values.length; i++)
      if (filter(values[i]))
        indexes.push(i);
    return indexes;
  }

  /**
   * Generates array from a range [begin; end] or [begin; end) if endExclusive.
   *
   * @export
   * @param {number} begin Beginning of the range.
   * @param {number} end End of the range.
   * @param {boolean} [endExclusive=false] Whether to exclude the end of the range.
   * @return {Int32Array} The range between begin and end.
   */
  export function genRange(begin: number, end: number, endExclusive = false): Int32Array {
    const nItems = end - begin + (endExclusive ? 0 : 1);
    const series = new Int32Array(nItems);
    for (let i = 0; i < nItems; ++i) {
      series[i] = begin + i;
    }
    return series;
  }
  
  /**  
   * Returns order of values as if they are sorted.
   *
   * @export
   * @param {any[]} values Input array.
   * @param {boolean} [reverse=false] Whether to return reversed order.
   * @return {number[]} The order generated.
   */
  export function argSort(values: any[], reverse: boolean = false): number[] {
    const sortfn = reverse ? (a: any[], b: any[]) => (b[0] - a[0]) : (a: any[], b: any[]) => (a[0] - b[0]);
    const decor = (v: any, i: number) => [v, i]; // set index to value
    const undecor = (a: any[]) => a[1]; // leave only index
    const _argsort = (arr: any[]) => arr.map(decor).sort(sortfn).map(undecor);
    return _argsort(values);
  }
}