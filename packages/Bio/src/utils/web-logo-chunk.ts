/** Residue frequencies per position of an MSA header chunk [start, end): for each position, a Map from residue to
 * its share of the rows ({@link oneOverCount} per occurrence). The result is the one the plain loop
 * `map.set(residue, (map.get(residue) || 0) + oneOverCount)` over the rows gives: every position of the chunk has a
 * Map, residues in first-seen order, and each share is {@link oneOverCount} added once per occurrence in row order,
 * so the floating-point sums are the same. Residues are counted by id in a typed array instead (a char-code table
 * for one-character residues), the Maps are built once at the end.
 * @param {Function} nextRow the next row after the given one (-1 to start), -1 at the end
 * @param {Function} getRegion the row's residues in [start, end) (an array or a string)
 * @param {number} start first position of the chunk
 * @param {number} end position after the chunk
 * @param {number} oneOverCount the share of one row
 * @return {Map<number, Map<string, number>>} */
export function getWebLogoChunkFrequencies(
  nextRow: (prevRowIdx: number) => number, getRegion: (rowIdx: number) => ArrayLike<string>,
  start: number, end: number, oneOverCount: number,
): Map<number, Map<string, number>> {
  const width = Math.max(0, end - start);
  const residues: string[] = [];
  const idOfResidue = new Map<string, number>();
  const idOfCode = new Int32Array(128).fill(-1);
  const order: number[][] = new Array(width);
  for (let pos = 0; pos < width; ++pos) order[pos] = [];
  let cap = 32;
  let counts = new Int32Array(width * cap);

  const newId = (residue: string): number => {
    const id = residues.length;
    residues.push(residue);
    if (id >= cap) { // more distinct residues than columns per position: widen the table
      const newCap = cap * 2;
      const newCounts = new Int32Array(width * newCap);
      for (let pos = 0; pos < width; ++pos)
        newCounts.set(counts.subarray(pos * cap, pos * cap + cap), pos * newCap);
      counts = newCounts;
      cap = newCap;
    }
    return id;
  };

  for (let rowIdx = -1; (rowIdx = nextRow(rowIdx)) !== -1;) {
    const region = getRegion(rowIdx);
    const n = region.length;
    if (n === 0) continue;
    if (n > width) // not a region of the chunk: the plain loop's result (or its error)
      return getWebLogoChunkFrequenciesPlain(nextRow, getRegion, start, end, oneOverCount);
    for (let pos = 0; pos < n; ++pos) {
      const residue = region[pos];
      let id: number;
      const code = typeof residue === 'string' && residue.length === 1 ? residue.charCodeAt(0) : 128;
      if (code < 128) {
        id = idOfCode[code];
        if (id < 0)
          idOfCode[code] = id = newId(residue);
      } else {
        const known = idOfResidue.get(residue);
        if (known === undefined)
          idOfResidue.set(residue, id = newId(residue));
        else
          id = known;
      }
      if (counts[pos * cap + id]++ === 0)
        order[pos].push(id);
    }
  }

  const res = new Map<number, Map<string, number>>();
  for (let pos = 0; pos < width; ++pos) {
    const residueMap = new Map<string, number>();
    for (const id of order[pos]) {
      const count = counts[pos * cap + id];
      let share = 0;
      for (let i = 0; i < count; ++i) share += oneOverCount;
      residueMap.set(residues[id], share);
    }
    res.set(start + pos, residueMap);
  }
  return res;
}

/** The plain Map loop {@link getWebLogoChunkFrequencies} stands for. */
export function getWebLogoChunkFrequenciesPlain(
  nextRow: (prevRowIdx: number) => number, getRegion: (rowIdx: number) => ArrayLike<string>,
  start: number, end: number, oneOverCount: number,
): Map<number, Map<string, number>> {
  const webLogoData: Map<number, Map<string, number>> = new Map();
  for (let i = start; i < end; i++)
    webLogoData.set(i, new Map());
  for (let i = -1; (i = nextRow(i)) !== -1;) {
    const seqChunk = getRegion(i);
    if (seqChunk.length === 0) continue;
    for (let pos = 0; pos < seqChunk.length; pos++) {
      const residue = seqChunk[pos];
      const residueMap = webLogoData.get(start + pos)!;
      residueMap.set(residue, (residueMap.get(residue) || 0) + oneOverCount);
    }
  }
  return webLogoData;
}
