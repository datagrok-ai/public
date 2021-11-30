import {assert} from './operations';

/**
 * Class to categorial encode/decode aligned amino acid residues sequence.
 *
 * @export
 * @class AlignedSequenceEncoder
 */
export class AlignedSequenceEncoder {
    protected aa2num: {[name: string]: number};
    protected num2aa: {[code: number]: string};

    constructor () {
        this.aa2num = {
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
        };
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
            assert(seq[l] == '-', `Wrong sequence format: ${termina[0]} without following '-' in '${seq}'.`)
            start = l;
        }
        if (seq.endsWith(termina[1])) {
            const l = termina[1].length+1; // Cut both 'COOH' and precending '-'.
            assert(seq[end-l] == '-', `Wrong sequence format: ${termina[1]} without '-' precending in '${seq}'.`)
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
        return seq.replace(/(-)([^-]+)/g, '$2')
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
        let values = new Array(nItems).fill(0);
    
        for (let i = 0; i < nItems; ++i) {
            let char = sequence[i];
    
            assert(char in this.aa2num, `Unknown char '${char}' found in sequence '${sequence}'`);
    
            values[i] = this.aa2num[char];
        }
        return values;
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
            let code = value[i];

            assert(code in this.num2aa, `Unknown code '${code}' found in vector '${value}'`);

            s += this.num2aa[code];
        }
        return s;
    }
}
