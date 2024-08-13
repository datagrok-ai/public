

// this is kinda too much memory at one glance, but,
// we have to keep column and rows representation of the matrix, otherwise it will be to slow
const indexesMapRows = new Map<number, Map<number, number>>();
const indexesMapCols = new Map<number, Map<number, number>>();

let nRows = 0;
let pruneValue = 0;


onmessage = (e) => {
  // init step for the worker
  if (e.data.is) {
    const {is, js, ds} = e.data;
    nRows = e.data.nRows;
    pruneValue = e.data.pruneValue;
    for (let i = 0; i < is.length; i++) {
      if (ds[i] < pruneValue)
        continue;

      if (!indexesMapRows.has(is[i]))
        indexesMapRows.set(is[i], new Map<number, number>());

      indexesMapRows.get(is[i])!.set(js[i], ds[i]);

      // same for the columns
      if (!indexesMapCols.has(js[i]))
        indexesMapCols.set(js[i], new Map<number, number>());
        indexesMapCols.get(js[i])!.set(is[i], ds[i]);
    }
    postMessage('ready');
    return;
  }

  // main step for the worker
  const {blockStart, blockEnd}: {[_: string]: number} = e.data;

  const resI: number[] = [];
  const resJ: number[] = [];
  const resD: number[] = [];

  for (let i = blockStart; i < blockEnd; i++) {
    const iMap = indexesMapRows.get(i);
    if (!iMap)
      continue;

    for (let j = 0; j < nRows; j++) {
      const jMap = indexesMapCols.get(j); // important: get the columns map
      // jMap.get(i) and iMap.get(j) are the same, but we need to check both
      if (!jMap || !iMap.has(j) || !jMap.has(i))
        continue;

      let sum = 0;
      const [iMapPointer, jMapPointer] = iMap.size > jMap.size ? [jMap, iMap] : [iMap, jMap];
      for (const [k, v] of iMapPointer) {
        if (jMapPointer.has(k))
          sum += v * jMapPointer.get(k)!;
      }
      //   if (sum > pruneValue) {
      resI.push(i);
      resJ.push(j);
      resD.push(sum);
    //   }
    }
  }


  postMessage({i: new Uint32Array(resI), j: new Uint32Array(resJ), distance: new Float32Array(resD)});
};
