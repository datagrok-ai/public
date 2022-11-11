import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';

import {Vector} from '@datagrok-libraries/utils/src/type-declarations';
import {vectorLength, vectorDotProduct} from '@datagrok-libraries/utils/src/vector-operations';
import {SeqPalette} from '../seq-palettes';
import {Aminoacids, AminoacidsPalettes} from '../aminoacids';
import {Nucleotides, NucleotidesPalettes} from '../nucleotides';
import {UnknownSeqPalettes} from '../unknown';
import {UnitsHandler} from '../utils/units-handler';

/** enum type to simplify setting "user-friendly" notation if necessary */
export const enum NOTATION {
  FASTA = 'fasta',
  SEPARATOR = 'separator',
  HELM = 'helm',
}

export const enum ALIGNMENT {
  SEQ_MSA = 'SEQ.MSA',
  SEQ = 'SEQ',
}

export const enum ALPHABET {
  DNA = 'DNA',
  RNA = 'RNA',
  PT = 'PT',
  UN = 'UN',
}

export const enum TAGS {
  aligned = 'aligned',
  alphabet = 'alphabet',
  alphabetSize = '.alphabetSize',
  alphabetIsMultichar = '.alphabetIsMultichar',
  separator = 'separator',
};

export type SeqColStats = { freq: MonomerFreqs, sameLength: boolean }
export type SplitterFunc = (seq: string) => string[];
export type MonomerFreqs = { [m: string]: number };

/** Stats of sequences with specified splitter func, returns { freq, sameLength }.
 * @param {DG.Column} seqCol
 * @param {number} minLength
 * @param {SplitterFunc} splitter
 * @return { SeqColStats }, sameLength: boolean } stats of column sequences
 */
export function getStats(seqCol: DG.Column, minLength: number, splitter: SplitterFunc): SeqColStats {
  const freq: { [m: string]: number } = {};
  let sameLength = true;
  let firstLength = null;

  for (const seq of seqCol.categories) {
    const mSeq = splitter(seq);

    if (firstLength == null)
      firstLength = mSeq.length;
    else if (mSeq.length !== firstLength)
      sameLength = false;

    if (mSeq.length >= minLength) {
      for (const m of mSeq) {
        if (!(m in freq))
          freq[m] = 0;
        freq[m] += 1;
      }
    }
  }
  return {freq: freq, sameLength: sameLength};
}

export const monomerRe: RegExp = /\[(\w+)\]|(\w)|(-)/g;

/** Split sequence for single character monomers, square brackets multichar monomer names or gap symbol.
 * @param {any} seq object with sequence
 * @return {string[]} array of monomers
 */
export function splitterAsFasta(seq: any): string[] {
  return wu<RegExpMatchArray>(seq.toString().matchAll(monomerRe))
    .map((ma: RegExpMatchArray) => {
      let mRes: string;
      const m: string = ma[0];
      if (m.length > 1) {
        mRes = ma[1];
      } else {
        mRes = m;
      }
      return mRes;
    }).toArray();
}

/** Gets method to split sequence by separator
 * @param {string} separator
 * @param limit
 * @return {SplitterFunc}
 */
export function getSplitterWithSeparator(separator: string, limit: number | undefined = undefined): SplitterFunc {
  return (seq: string) => {
    return seq.split(separator, limit);
  };
}

const helmRe: RegExp = /(PEPTIDE1|DNA1|RNA1)\{([^}]+)}/g;
const helmPp1Re: RegExp = /\[([^\[\]]+)]/g;


/** Splits Helm string to monomers, but does not replace monomer names to other notation (e.g. for RNA).
 * Only for linear polymers, does not split RNA for ribose and phosphate monomers.
 * @param {string} seq Source string of HELM notation
 * @return {string[]}
 */
export function splitterAsHelm(seq: any): string[] {
  helmRe.lastIndex = 0;
  const ea: RegExpExecArray | null = helmRe.exec(seq.toString());
  const inSeq: string | null = ea ? ea[2] : null;

  const mmPostProcess = (mm: string): string => {
    helmPp1Re.lastIndex = 0;
    const pp1M = helmPp1Re.exec(mm);
    if (pp1M && pp1M.length >= 2)
      return pp1M[1];
    else
      return mm;
  };

  const mmList: string[] = inSeq ? inSeq.split('.') : [];
  return mmList.map(mmPostProcess);
}

/** Get splitter method to split sequences to monomers
 * @param {string} units
 * @param {string} separator
 * @param limit
 * @return {SplitterFunc}
 */
export function getSplitter(units: string, separator: string, limit: number | undefined = undefined): SplitterFunc {
  if (units.toLowerCase().startsWith(NOTATION.FASTA))
    return splitterAsFasta;
  else if (units.toLowerCase().startsWith(NOTATION.SEPARATOR))
    return getSplitterWithSeparator(separator, limit);
  else if (units.toLowerCase().startsWith(NOTATION.HELM))
    return splitterAsHelm;
  else
    throw new Error(`Unexpected units ${units} .`);

  // TODO: Splitter for HELM
}

/** Generate splitter function for sequence column
 * @param {DG.Column} col
 * @return {SplitterFunc} Splitter function
 */
export function getSplitterForColumn(col: DG.Column): SplitterFunc {
  if (col.semType !== DG.SEMTYPE.MACROMOLECULE)
    throw new Error(`Get splitter for semType "${DG.SEMTYPE.MACROMOLECULE}" only.`);

  const units = col.getTag(DG.TAGS.UNITS);
  const separator = col.getTag(TAGS.separator);
  return getSplitter(units, separator);
}

const longMonomerPartRe: RegExp = /(\w+)/g;

/** Convert long monomer names to short ones */
export function monomerToShort(amino: string, maxLengthOfMonomer: number): string {
  const shortAminoMatch: RegExpMatchArray | null = amino.match(longMonomerPartRe);
  const needAddDots: boolean = amino.length > maxLengthOfMonomer || (shortAminoMatch?.length ?? 0) > 1;
  const shortAmino = shortAminoMatch?.[0] ?? ' ';
  return !needAddDots ? shortAmino : shortAmino.substring(0, maxLengthOfMonomer) + '…';
}

/** Calculate similarity in current sequence and alphabet.
 * @param {MonomerFreqs} freq
 * @param {Set<string>} alphabet
 * @param {string} gapSymbol
 * @return {number} Cosine similarity
 */
export function getAlphabetSimilarity(freq: MonomerFreqs, alphabet: Set<string>, gapSymbol: string = '-'): number {
  const keys = new Set<string>([...new Set(Object.keys(freq)), ...alphabet]);
  keys.delete(gapSymbol);

  const freqA: number[] = [];
  const alphabetA: number[] = [];
  for (const m of keys) {
    freqA.push(m in freq ? freq[m] : 0);
    alphabetA.push(alphabet.has(m) ? 1 : 0);
  }
  /* There were a few ideas: chi-squared, pearson correlation (variance?), scalar product */
  const freqV: Vector = new Vector(freqA);
  const alphabetV: Vector = new Vector(alphabetA);
  return vectorDotProduct(freqV, alphabetV) / (vectorLength(freqV) * vectorLength(alphabetV));
}

export function detectAlphabet(stats: SeqColStats): string {
  const alphabetCandidates: [string, Set<string>][] = [
    [ALPHABET.PT, UnitsHandler.PeptideFastaAlphabet],
    [ALPHABET.DNA, UnitsHandler.DnaFastaAlphabet],
    [ALPHABET.RNA, UnitsHandler.RnaFastaAlphabet],
  ];

  // Calculate likelihoods for alphabet_candidates
  const alphabetCandidatesSim: number[] = alphabetCandidates.map(
    (c) => getAlphabetSimilarity(stats.freq, c[1]));
  const maxCos = Math.max(...alphabetCandidatesSim);
  const alphabet = maxCos > 0.65 ? alphabetCandidates[alphabetCandidatesSim.indexOf(maxCos)][0] : 'UN';
  return alphabet;
}

/** Selects a suitable palette based on column data
 * @param {DG.Column} seqCol Column to look for a palette
 * @param {number}  minLength minimum length of sequence to detect palette (empty strings are allowed)
 * @return {SeqPalette} Palette corresponding to the alphabet of the sequences in the column
 */
export function pickUpPalette(seqCol: DG.Column, minLength: number = 5): SeqPalette {
  let alphabet: string;
  if (seqCol.semType == DG.SEMTYPE.MACROMOLECULE) {
    const uh: UnitsHandler = new UnitsHandler(seqCol);
    alphabet = uh.alphabet;
  } else {
    const stats: SeqColStats = getStats(seqCol, minLength, splitterAsFasta);
    alphabet = detectAlphabet(stats);
  }

  const res = getPaletteByType(alphabet);
  return res;
}

export function getPaletteByType(paletteType: string): SeqPalette {
  switch (paletteType) {
  case 'PT':
    return AminoacidsPalettes.GrokGroups;
  case 'NT':
  case 'DNA':
  case 'RNA':
    return NucleotidesPalettes.Chromatogram;
    // other
  default:
    return UnknownSeqPalettes.Color;
  }
}

export function pickUpSeqCol(df: DG.DataFrame): DG.Column | null {
  const semTypeColList = df.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  let resCol: DG.Column | null = semTypeColList.find((col) => {
    const units = col.getTag(DG.TAGS.UNITS);
    return units ? units.indexOf('MSA') !== -1 : false;
  }) ?? null;
  if (!resCol && semTypeColList.length > 0)
    resCol = semTypeColList[0];
  return resCol;
}
