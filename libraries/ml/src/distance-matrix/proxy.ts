

/** Proxy for DistanceMatrix class. Allows to index matrix as matrix[i][j]
 * Note: much slower than direct indexing, but still much faster than recalculating distances.
 * will be used for T-SNE mainly.
 * @param {Float32Array}condensedArray - array of distances between all pairs of objects
 * @param {number}size - number of comparebles in matrix
 * @return {Float32Array} - proxy for condensedArray
*/
export function distanceMatrixProxy(condensedArray: Float32Array, size: number): Float32Array {
  const linearFunc = dmLinearIndex(size);
  function linearIndex(i: symbol | string, j: symbol | string) {
    const iNum = Number(i);
    const jNum = Number(j);
    return linearFunc(iNum, jNum);
  }

  function idx2Handler(idx1: symbol | string):ProxyHandler<Float32Array> {
    return (
      {
        get(target, idx2, _receiver) {
          if (idx1 === idx2) return 0;
          const linearIdx = Number(idx1) > Number(idx2) ? linearIndex(idx2, idx1) : linearIndex(idx1, idx2);
          return target[linearIdx];
        },
      }
    );
  }
  const idx1Handler: ProxyHandler<Float32Array> = {
    get(target, idx1, _receiver) {
      if (idx1 === 'length') return size;
      return new Proxy(target, idx2Handler(idx1));
    },
  };

  return new Proxy(condensedArray, idx1Handler);
}

export function dmLinearIndex(size: number): (i: number, j: number) => number {
  return (i: number, j: number) => size * i + j - Math.floor(((i + 2) * (i + 1)) / 2);
}
