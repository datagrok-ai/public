import {ISeqSplitted, MonomerFreqs, SeqColStats} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

/** The monomer frequencies and the same-length flag of a sequence column, {@link SeqHandler.stats}: every row's
 * canonical monomers counted in row order. The counts go to a prototype-less object first (no `in` test that walks
 * the prototype chain per monomer); the result is the object the plain loop `if (!(m in freq)) freq[m] = 0;
 * freq[m] += 1` makes, with its keys in the same order.
 * @param {number} rowCount
 * @param {Function} getSplitted the row's split sequence
 * @return {SeqColStats} */
export function getSeqColStats(rowCount: number, getSplitted: (rowIdx: number) => ISeqSplitted): SeqColStats {
  const counts: MonomerFreqs = Object.create(null);
  let sameLength = true;
  let firstLength: number | null = null;
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    const mSeq: ISeqSplitted = getSplitted(rowIdx);
    const length = mSeq.length;
    if (firstLength == null)
      firstLength = length;
    else if (length !== firstLength)
      sameLength = false;

    for (let posIdx = 0; posIdx < length; ++posIdx) {
      const cm = mSeq.getCanonical(posIdx);
      const count = counts[cm];
      counts[cm] = count === undefined ? 1 : count + 1;
    }
  }
  return {freq: toFreqObject(counts), sameLength: sameLength};
}

/** A plain object with the counts, keys in first-seen order. A monomer named like an Object.prototype member
 * (`toString`, ...) gets what the plain loop gave it: its `in` test found the inherited member, so each count was
 * added to that. */
function toFreqObject(counts: MonomerFreqs): MonomerFreqs {
  const freq: MonomerFreqs = {};
  for (const m of Object.keys(counts)) {
    const count = counts[m];
    if (m in Object.prototype) {
      if (!(m in freq)) freq[m] = 0;
      for (let i = 0; i < count; ++i) freq[m] += 1;
    } else
      freq[m] = count;
  }
  return freq;
}
