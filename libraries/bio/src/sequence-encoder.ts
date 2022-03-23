import {assert} from '@datagrok-libraries/utils/src/vector-operations';

type SideChainScale = {[name: string]: number};
type SideChainScaleCollection = {[name: string]: SideChainScale};

class SideChainScales {
  static scales: SideChainScaleCollection = {
    // Wimley-White interfacial hydrophobicity scale
    'WimleyWhite': {
      '-': 0,
      'A': 0.17,
      'C': -0.24,
      'D': -0.07, // Asp-: 1.23
      'E': -0.01, // Glu-: 2.02
      'F': -1.13, //
      'G': 0.01,
      'H': 0.17, // His+: 0.96
      'I': -0.31,
      'K': 0.99, // Lys+
      'L': -0.56,
      'M': -0.23,
      'N': 0.42,
      'P': 0.45,
      'Q': 0.58,
      'R': 0.81, // Arg+
      'S': 0.13,
      'T': 0.14,
      'V': 0.07,
      'W': -1.85,
      'Y': -0.94,
    },
    'categorial': {
      '-': 0,
      'A': 1,
      'C': 2,
      'D': 3,
      'E': 4,
      'F': 5,
      'G': 6,
      'H': 7,
      'I': 8,
      'K': 9,
      'L': 10,
      'M': 11,
      'N': 12,
      'P': 13,
      'Q': 14,
      'R': 15,
      'S': 16,
      'T': 17,
      'V': 18,
      'W': 19,
      'Y': 20,
    },
  };

  static getAvailableScales(): string[] {
    return Object.entries(this.scales).map(([k, _]) => k);
  }

  static getScale(name: string): SideChainScale {
    assert(!(this.scales[name] === undefined), `Scale '${name}' was not found.`);
    return this.scales[name];
  }
}

/**
 * Class to categorial encode/decode aligned amino acid residues sequence.
 *
 * @export
 * @class AlignedSequenceEncoder
 */
export class AlignedSequenceEncoder {
  protected aa2num: SideChainScale;
  protected num2aa: {[code: number]: string};

  constructor(scale: string = 'categorial') {
    this.aa2num = SideChainScales.getScale(scale);
    this.num2aa = {};
    Object.entries(this.aa2num).forEach(([k, v]) => (this.num2aa[v] = k));
  }

  /**
     * Truncate NH2 and -COOH terminals of the given sequence.
     *
     * @static
     * @param {string} seq The sequence provided.
     * @return {string} Truncated sequence.
     * @memberof AlignedSequenceEncoder
     */
  static _truncateSequence(seq: string): string {
    let start = 0;
    let end = seq.length;
    const termina = ['NH2', 'COOH'];

    if (seq.startsWith(termina[0])) {
      const l = termina[0].length; // Cut only 'NH2' without following '-'.
      assert(seq[l] == '-', `Wrong sequence format: ${termina[0]} without following '-' in '${seq}'.`);
      start = l;
    }
    if (seq.endsWith(termina[1])) {
      const l = termina[1].length+1; // Cut both 'COOH' and precending '-'.
      assert(seq[end-l] == '-', `Wrong sequence format: ${termina[1]} without '-' precending in '${seq}'.`);
      end -= l;
    }
    return seq.substring(start, end);
  }

  /**
     * Cuts auxiliary defises before a residue.
     *
     * @static
     * @param {string} seq The sequence to process.
     * @return {string} Processed sequence.
     * @memberof AlignedSequenceEncoder
     */
  static _dropDefises(seq: string): string {
    return seq.replace(/(-)([^-]+)/g, '$2');
  }

  /**
     * Performs truncation and cutting auxiliary defises.
     *
     * @static
     * @param {string} sequence The sequence work under process.
     * @return {string} Result of cleaning.
     * @memberof AlignedSequenceEncoder
     */
  static clean(sequence: string): string {
    return AlignedSequenceEncoder._dropDefises(AlignedSequenceEncoder._truncateSequence(sequence));
  }

  /**
     * Categorial encode of the sequence provided.
     *
     * @param {string} sequence The sequence.
     * @return {number[]} Encoded vector.
     * @memberof AlignedSequenceEncoder
     */
  public encode(sequence: string): number[] {
    const nItems = sequence.length;
    const values = new Array(nItems).fill(0);

    for (let i = 0; i < nItems; ++i) {
      const char = sequence[i];

      assert(char in this.aa2num, `Unknown char '${char}' found in sequence '${sequence}'`);

      values[i] = this.encodeLettter(char);
    }
    return values;
  }

  public encodeLettter(letter: string): number {
    return this.aa2num[letter];
  }

  /**
     * Decode the encoded vector into the sequence back.
     *
     * @param {number[]} value The vector encoded.
     * @return {string} Decoded sequence.
     * @memberof AlignedSequenceEncoder
     */
  public decode(value: number[]): string {
    let s: string = '';

    for (let i = 0; i < value.length; ++i) {
      const code = value[i];

      assert(code in this.num2aa, `Unknown code '${code}' found in vector '${value}'`);

      s += this.num2aa[code];
    }
    return s;
  }
}
