export function kFold(rowCount: number, folds: number, seed: number): Int32Array {
  const random = mulberry32(seed);
  const order = new Int32Array(rowCount);
  for (let i = 0; i < rowCount; i++)
    order[i] = i;
  for (let i = rowCount - 1; i > 0; i--) {
    const j = Math.floor(random() * (i + 1));
    [order[i], order[j]] = [order[j], order[i]];
  }
  const fold = new Int32Array(rowCount);
  for (let j = 0; j < rowCount; j++)
    fold[order[j]] = j % folds;
  return fold;
}

function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
