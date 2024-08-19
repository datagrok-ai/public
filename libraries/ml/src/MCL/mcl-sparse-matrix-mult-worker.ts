

const indexesMapRows = new Map<number, Map<number, number>>();

onmessage = (e) => {
  // init step for the worker
  if (e.data.is) {
    const {is, js, ds} = e.data;
    for (let i = 0; i < is.length; i++) {
      if (!indexesMapRows.has(is[i]))
        indexesMapRows.set(is[i], new Map<number, number>());

      indexesMapRows.get(is[i])!.set(js[i], ds[i]);
    }
    postMessage('ready');
    return;
  }

  // main step for the worker
  const {blockStart, blockEnd}: {[_: string]: number} = e.data;

  const resI: number[] = [];
  const resJ: number[] = [];
  const resD: number[] = [];

  // we do operations row by row
  // therefore, iterating over all columns is necessary
  // graphically speaking, one worker is working on a horizontal stripe in a matrix


  // blockstart and blockend are the indexes of the rows that the worker is responsible for
  for (let i = blockStart; i < blockEnd; i++) {
    const iMap = indexesMapRows.get(i);
    if (!iMap)
      continue;

    // this worker is responsible for every nonzero value in the row (or iMap)
    // since we are using iMap a lot, we can materialize it
    const iMapEntries = Array.from(iMap.entries());
    for (const [k, _v] of iMapEntries) {
      let sum = 0;
      // here we know that the value is not zero
      // therefore, we can iterate over the columns
      for (const [k1, v1] of iMapEntries) {
        // here we check the column
        const jVal = indexesMapRows.get(k1)?.get(k);
        if (jVal)
          sum += v1 * jVal;
      }

      resI.push(i);
      resJ.push(k);
      resD.push(sum);
    }
  }


  postMessage({i: new Uint32Array(resI), j: new Uint32Array(resJ), distance: new Float32Array(resD)});
};
