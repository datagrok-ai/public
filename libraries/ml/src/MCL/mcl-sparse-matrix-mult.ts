import {SparseMatrix} from '@datagrok-libraries/math';

// simple semaphore
let semaphore = Promise.resolve();

const lock = async () => {
  await semaphore;
  let release: () => void = () => {};
  const promise = new Promise<void>((resolve) => {
    release = resolve;
  });
  semaphore = promise;
  return release;
};


/** Multiplies sparse matrix by itself. note that it expects that sparse matrix is full
 * (i.e. not only the upper right corner) and it has self loops if any */
export async function multSparseMatrix(
  sparseMatrix: SparseMatrix, nRows: number, pruneValue: number
): Promise<SparseMatrix> {
  const workersNum = Math.min(Math.max(navigator.hardwareConcurrency - 2, 1), nRows);

  // number of horizontal strips that we will divide the matrix into
  // we divide the matrix into horizontal strips, because its faster for indexing and better for cache
  const numOfHorizontalStrinps = workersNum * 5;

  const indexStarts = new Uint32Array(numOfHorizontalStrinps);
  const indexEnds = new Uint32Array(numOfHorizontalStrinps);

  const blockHeight = nRows / numOfHorizontalStrinps;
  for (let i = 0; i < numOfHorizontalStrinps; i++) {
    indexStarts[i] = Math.floor(i * blockHeight);
    indexEnds[i] = i === numOfHorizontalStrinps - 1 ? nRows : Math.floor((i + 1) * blockHeight);
  }


  const workers: Worker[] = new Array(workersNum).fill(0)
    .map(() => new Worker(new URL('./mcl-sparse-matrix-mult-worker', import.meta.url)));


  const availableIndexes = new Set<number>();
  for (let i = 0; i < numOfHorizontalStrinps; i++)
    availableIndexes.add(i);

  const initPromises = workers.map((worker) => {
    return new Promise<void>((resolve) => {
      worker.postMessage({
        is: sparseMatrix.i,
        js: sparseMatrix.j,
        ds: sparseMatrix.distance,
        nRows,
        pruneValue
      });

      worker.onmessage = () => {
        resolve();
      };
    });
  });
  await Promise.all(initPromises);

  const res: SparseMatrix[] = [];

  const takeChunk = async (workerIdx: number) => {
    const release = await lock();
    const index = availableIndexes.values().next().value;
    if (index == null) {
      release();
      return;
    }
    availableIndexes.delete(index);
    release();
    await new Promise<void>((resolve) => {
      workers[workerIdx].postMessage({
        blockStart: indexStarts[index],
        blockEnd: indexEnds[index],
        workerIdx
      });

      workers[workerIdx].onmessage = (e) => {
        res.push(e.data);
        resolve();
      };
    });
    await takeChunk(workerIdx);
  };


  const promises = workers.map((_worker, i) => {
    return takeChunk(i);
  });
  await Promise.all(promises);
  workers.forEach((worker) => worker.terminate());


  let totLen = 0;
  for (const r of res)
    totLen += r.i.length;

  const i = new Uint32Array(totLen);
  const j = new Uint32Array(totLen);
  const d = new Float32Array(totLen);
  let offset = 0;
  for (const r of res) {
    i.set(r.i, offset);
    j.set(r.j, offset);
    d.set(r.distance, offset);
    offset += r.i.length;
  }

  // don't forget to get the updated values for the diagonal, their rows need to be basically squared and summed
  // also, the diagonal values are not pruned

  return {i, j, distance: d};
}
