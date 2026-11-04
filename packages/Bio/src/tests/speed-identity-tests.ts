import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {StringListSeqSplitted, splitterAsBiln} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {GAP_SYMBOL} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

import {DrawStyle, getPrintColumnTemps, printLeftOrCentered} from '@datagrok-libraries/bio/src/utils/cell-renderer';

import {getSeqColStats} from '../utils/seq-helper/seq-stats';
import {getWebLogoChunkFrequencies, getWebLogoChunkFrequenciesPlain} from '../utils/web-logo-chunk';

// The speed rewrites of 2026-11-02 against the plain loops they replace: the same results, key order included.

/** Monomers named like Object.prototype members, integer-like, non-ASCII, multi-character, empty, gaps */
const TOKENS = ['A', 'toString', '__proto__', 'constructor', 'valueOf', '1', '10', '2', '01', 'é', 'Ω', 'ab', '', '-',
  'x y', '\u{1F600}', 'C', 'G', ...Array.from({length: 40}, (_, i) => `m${i}`)];

function randomRows(n: number, seed: number): string[][] {
  const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
  return Array.from({length: n}, () =>
    Array.from({length: Math.floor(rnd() * 30)}, () => TOKENS[Math.floor(rnd() * TOKENS.length)]));
}

function plainStats(splitted: ISeqSplitted[]): { freq: { [m: string]: number }, sameLength: boolean } {
  const freq: { [m: string]: number } = {};
  let sameLength = true;
  let firstLength = null;
  for (const mSeq of splitted) {
    if (firstLength == null) firstLength = mSeq.length;
    else if (mSeq.length !== firstLength) sameLength = false;
    for (let posIdx = 0; posIdx < mSeq.length; ++posIdx) {
      const cm = mSeq.getCanonical(posIdx);
      if (!(cm in freq)) freq[cm] = 0;
      freq[cm] += 1;
    }
  }
  return {freq, sameLength};
}

const chunkToString = (m: Map<number, Map<string, number>>): string =>
  JSON.stringify([...m.entries()].map(([pos, mm]) => [pos, [...mm.entries()]]));

category('SeqHandler: stats', () => {
  test('as the plain loop', async () => {
    for (const [rows, seed] of [[0, 1], [1, 2], [300, 3]]) {
      const splitted = randomRows(rows, seed).map((ms) => new StringListSeqSplitted(ms, '-'));
      const exp = plainStats(splitted);
      const res = getSeqColStats(splitted.length, (i) => splitted[i]);
      expect(JSON.stringify(res), JSON.stringify(exp));
      expect(Object.keys(res.freq).join('|'), Object.keys(exp.freq).join('|'));
    }
  });
});

category('WebLogo.chunk', () => {
  test('frequencies as the plain loop', async () => {
    const rows = randomRows(500, 7);
    for (const [start, end] of [[0, 10], [5, 25], [28, 40]]) {
      // regions as the splitted gives them, every other row filtered out
      const regions = rows.map((ms) => ms.slice(start, end));
      const nextRow = (i: number) => {const j = i + 2 - ((i + 2) % 2); return j < regions.length ? j : -1;};
      const getRegion = (i: number) => regions[i];
      const share = 1 / 250;
      expect(chunkToString(getWebLogoChunkFrequencies(nextRow, getRegion, start, end, share)),
        chunkToString(getWebLogoChunkFrequenciesPlain(nextRow, getRegion, start, end, share)));
      // one-character residues read from a string region (FASTA)
      const strRegions = regions.map((r) => r.filter((m) => m.length === 1).join(''));
      expect(chunkToString(getWebLogoChunkFrequencies(nextRow, (i) => strRegions[i], start, end, share)),
        chunkToString(getWebLogoChunkFrequenciesPlain(nextRow, (i) => strRegions[i], start, end, share)));
    }
  });
});

category('splitters: BILN', () => {
  test('originals and canonicals', async () => {
    const reOriginal = /^\[(.*)\]/;
    const reMark = /\(\d{1,2},\d{1,2}\)$/;
    for (const biln of ['A-C(1,3)-G.H-C(1,3)', '[D-Arg]-[Lys(Me)2](2,3)-[x-y-z]-C(2,3)', '[a](1,2)(3,3)-X)-Q(123,3)',
      '[[A-B]]-[]-[ab]c-a[b]-z(1,1))', 'A-"-"-\'-\'-C']) {
      const ss = splitterAsBiln(biln);
      const raw = ss.getOriginalRegion(0, ss.length);
      for (let k = 0; k < ss.length; ++k) {
        const original = raw[k].replace(reOriginal, '$1');
        let canonical = original;
        while (reMark.test(canonical)) canonical = canonical.replace(reMark, '');
        expect(ss.getOriginal(k), original);
        expect(ss.getCanonical(k), ss.isGap(k) ? GAP_SYMBOL : canonical);
      }
    }
  });
});

category('renderers: print options', () => {
  test('column temps given once draw as read per call', async () => {
    /** A context recording what is drawn and set */
    const recorder = (log: string[]): CanvasRenderingContext2D => new Proxy({} as any, {
      get: (t, k) => k === 'measureText' ? (s: string) => ({width: String(s).length * 7, fontBoundingBoxAscent: 9,
        fontBoundingBoxDescent: 3}) : typeof k === 'string' && ['fillText', 'fillRect', 'save', 'restore'].includes(k) ?
        (...a: any[]) => {log.push(`${k}(${a.join(',')})`);} : t[k],
      set: (t, k, v) => {log.push(`${String(k)}=${v}`); t[k] = v; return true;},
    });
    for (const temps of [{}, {'color-code': false}, {'compare-with-current': false}, {'highlight-difference': 'equal'}]) {
      const gridCell = {cell: {column: {temp: temps}}} as any;
      for (const [s, wordIdx, ref] of [['A', 0, ['A', 'C']], ['C', 1, ['A', 'C']], ['Lys_Me', 2, []]] as const) {
        const opts = {color: 'rgb(10,20,30)', left: true, transparencyRate: 0, separator: '-', last: false,
          drawStyle: DrawStyle.classic, maxWord: [0, 10, 20], wordIdx, gridCell, referenceSequence: [...ref] as string[],
          maxLengthOfMonomer: 3, monomerTextSizeMap: {}, selectedPosition: 2};
        const perCall: string[] = [];
        const once: string[] = [];
        const x1 = printLeftOrCentered(recorder(perCall), s, 5, 5, 100, 20, {...opts});
        const x2 = printLeftOrCentered(recorder(once), s, 5, 5, 100, 20,
          {...opts, monomerTextSizeMap: {}, columnTemps: getPrintColumnTemps(gridCell)});
        expect(x2, x1);
        expect(once.join('|'), perCall.join('|'));
      }
    }
  });
});
