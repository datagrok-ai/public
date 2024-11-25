export type Linkage = {
  fChain: number,
  sChain: number,
  /** Continuous 1-based numbering */ fMonomer: number,
  /** Continuous 1-based numbering */ sMonomer: number,
  fR: number,
  sR: number
}

/** Gets 0-based in-index (simple polymer) of out-index (continuous) {@link idx} */
export function getInnerIdx(outIdx: number, monomers: string[][]): [number, number] {
  // let prevSpCount = 0;
  // for (let spI = 0; spI < monomers.length && idx >= (prevSpCount + monomers[spI].length); ++spI)
  //   prevSpCount += monomers[spI].length;
  // return idx - prevSpCount;
  let inIdx = outIdx;
  let spIdx: number;
  for (spIdx = 0; spIdx < monomers.length && inIdx >= monomers[spIdx].length; ++spIdx)
    inIdx -= monomers[spIdx].length;
  return [inIdx, spIdx];
}

/** Gets 0-based out-index of 0-based in-index {@link inIdx} monomer of simple polymer {@link spIdx} */
export function getOuterIdx(inIdx: number, spIdx: number, monomers: string[][]): number {
  let outIdx = 0;
  for (let i = 0; i < spIdx; ++i)
    outIdx += monomers[i].length;
  return outIdx + inIdx;
}
