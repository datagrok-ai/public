import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {
  AsmRole, AsmSequence, packBatch, unpackRow,
} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import {buildSeqBatch, toAsmSequence} from '@datagrok-libraries/bio/src/monomer-works/assembly-input';
import {
  getEffectiveMonomers, getMonomerLibJson, getMonomerLibKey, getSeqColMonomerLib,
} from '@datagrok-libraries/bio/src/monomer-works/library-json';
import {getMonomerLibHelper, IMonomerLib, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {splitterAsBiln} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

/** The input layer of the macromolecule assembler (Bio on Crux): any notation to the parsed sequence with exact
 * library keys, per-row errors, the packed batch and the effective library JSON. No assembler is called. */
category('assemblyInput', () => {
  let seqHelper: ISeqHelper;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let monomerLib: IMonomerLib;

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests(); // HELMCoreLibrary.json + polytool-lib.json
    monomerLib = monomerLibHelper.getMonomerLib();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  function seqCol(list: string[], units: NOTATION, alphabet?: ALPHABET, separator?: string): DG.Column<string> {
    const col = DG.Column.fromList(DG.TYPE.STRING, 'seq', list) as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = units;
    col.setTag(bioTAGS.aligned, ALIGNMENT.SEQ);
    if (alphabet) col.setTag(bioTAGS.alphabet, alphabet);
    if (separator) col.setTag(bioTAGS.separator, separator);
    return col;
  }

  function rowOf(col: DG.Column<string>, rowIdx: number = 0): AsmSequence | string {
    const sh = seqHelper.getSeqHandler(col);
    const helm = sh.notation === NOTATION.HELM ? col.get(rowIdx)! : undefined;
    const res = toAsmSequence(sh.getSplitted(rowIdx), sh, monomerLib, {helm});
    return res.seq ?? res.error;
  }

  function seqOf(col: DG.Column<string>, rowIdx: number = 0): AsmSequence {
    const res = rowOf(col, rowIdx);
    if (typeof res === 'string') throw new Error(`Unexpected row error: ${res}`);
    return res;
  }

  function errorOf(col: DG.Column<string>, rowIdx: number = 0): string {
    const res = rowOf(col, rowIdx);
    if (typeof res !== 'string') throw new Error(`Row error expected, got ${JSON.stringify(res)}`);
    return res;
  }

  const symbolsOf = (s: AsmSequence) => s.monomers.map((m) => m.gap ? '-' : (m.symbol || `{${m.inlineSmiles}}`));
  const S = AsmRole.SUGAR;
  const B = AsmRole.BASE;
  const L = AsmRole.LINKER;

  test('packRoundTrip', async () => {
    const s: AsmSequence = {monomers: [
      {symbol: 'A', polymerType: 'PEPTIDE', role: AsmRole.AA, position: 0},
      {symbol: '', polymerType: 'PEPTIDE', role: AsmRole.AA, gap: true, position: 1},
      {symbol: 'C', polymerType: 'PEPTIDE', role: AsmRole.AA, position: 2},
    ], graphInfo: {disjointSeqStarts: [0], connections: [
      {seqIndex1: 0, monomerIndex1: 0, rGroup1: 1, seqIndex2: 0, monomerIndex2: 2, rGroup2: 10}]}};
    const batch = packBatch([s, null, s], [null, null, 'unknown monomer PEPTIDE:Xyz']);
    expectArray(batch.rowError, [null, 'no sequence', 'unknown monomer PEPTIDE:Xyz']);
    expect(batch.symbols.length, 2);
    const u = unpackRow(batch, 0);
    expectArray(symbolsOf(u), ['A', '-', 'C']);
    expect(u.graphInfo.connections[0].rGroup2, 10);
    expect(unpackRow(batch, 2).monomers.length, 0);
  });

  test('fastaPeptide', async () => {
    const s = seqOf(seqCol(['AC-D'], NOTATION.FASTA, ALPHABET.PT));
    expectArray(symbolsOf(s), ['A', 'C', '-', 'D']);
    expectArray(s.monomers.map((m) => m.position), [0, 1, 2, 3]);
    expectArray(s.monomers.map((m) => m.polymerType), ['PEPTIDE', 'PEPTIDE', 'PEPTIDE', 'PEPTIDE']);
  });

  test('fastaDnaExpanded', async () => {
    // Every base gets the default sugar (d for DNA) and a phosphate before its sugar after the first nucleotide;
    // a gap adds no sugar or phosphate (the old linear path left a dangling phosphate after a leading gap)
    const s = seqOf(seqCol(['-AC-G'], NOTATION.FASTA, ALPHABET.DNA));
    expectArray(symbolsOf(s), ['-', 'd', 'A', 'p', 'd', 'C', '-', 'p', 'd', 'G']);
    expectArray(s.monomers.map((m) => m.position), [0, 1, 1, 2, 2, 2, 3, 4, 4, 4]);
    expectArray(s.monomers.filter((m) => !m.gap).map((m) => m.role), [S, B, L, S, B, L, S, B]);
  });

  test('fastaRnaSugar', async () => {
    const s = seqOf(seqCol(['AU'], NOTATION.FASTA, ALPHABET.RNA));
    expectArray(symbolsOf(s), ['r', 'A', 'p', 'r', 'U']);
  });

  test('separatorNotABase', async () => {
    // The old linear path attached r/p tokens to a default sugar as bases; now it is the row's error
    const err = errorOf(seqCol(['A.p.C'], NOTATION.SEPARATOR, ALPHABET.DNA, '.'));
    expect(err.includes('not a nucleobase'), true, err);
  });

  test('helmRnaKeys', async () => {
    // D15: RNA R/P resolve to the library's r/p
    const s = seqOf(seqCol(['RNA1{R(A)P.R(C)P.R(G)}$$$$V2.0'], NOTATION.HELM));
    expectArray(symbolsOf(s), ['r', 'A', 'p', 'r', 'C', 'p', 'r', 'G']);
    expectArray(s.monomers.map((m) => m.role), [S, B, L, S, B, L, S, B]);
    expectArray(s.monomers.map((m) => m.position), [0, 1, 2, 3, 4, 5, 6, 7]);
  });

  test('helmCyclic', async () => {
    const s = seqOf(seqCol(['PEPTIDE1{C.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$V2.0'], NOTATION.HELM));
    expect(s.graphInfo.connections.length, 1);
    const c = s.graphInfo.connections[0];
    expectArray([c.seqIndex1, c.monomerIndex1, c.rGroup1, c.seqIndex2, c.monomerIndex2, c.rGroup2], [0, 0, 3, 0, 3, 3]);
  });

  test('helmDuplexPairs', async () => {
    // Hydrogen bonds add no covalent bond: two separate chains, no error
    const s = seqOf(seqCol(['RNA1{R(A)P.R(C)}|RNA2{R(G)P.R(U)}$RNA1,RNA2,2:pair-5:pair$$$V2.0'], NOTATION.HELM));
    expectArray(s.graphInfo.disjointSeqStarts, [0, 5]);
  });

  test('helmInlineSmiles', async () => {
    const s = seqOf(seqCol(['PEPTIDE1{[C[C@H](N[*:1])C(=O)[*:2]].G}$$$$V2.0'], NOTATION.HELM));
    expect(s.monomers[0].symbol, '');
    expect(s.monomers[0].inlineSmiles, 'C[C@H](N[*:1])C(=O)[*:2]');
  });

  test('helmGap', async () => {
    const s = seqOf(seqCol(['PEPTIDE1{A.*.C}$$$$V2.0'], NOTATION.HELM));
    expectArray(symbolsOf(s), ['A', '-', 'C']);
    expectArray(s.monomers.map((m) => m.position), [0, 1, 2]);
  });

  test('rowErrors', async () => {
    const col = seqCol([
      'PEPTIDE1{A.[Xyz].C}$$$$V2.0',
      'PEPTIDE1{A.C}$PEPTIDE1,PEPTIDE1,?:R1-2:R2$$$V2.0',
      'BLOB1{Bead}|PEPTIDE1{A}$$$$V2.0',
      'PEPTIDE1{A.(G,C).A}$$$$V2.0',
      'PEPTIDE1{A\'2-5\'.G}$$$$V2.0',
      'PEPTIDE1{' + new Array(120).fill('A').join('.') + '}$PEPTIDE1,PEPTIDE1,120:R2-1:R1$$$V2.0',
      'PEPTIDE1{A.C.D}$$$$V2.0#',
      'PEPTIDE1{Ac(1).A.A}$$$$',
    ], NOTATION.HELM);
    expect(errorOf(col, 0), 'unknown monomer PEPTIDE:Xyz');
    expect(errorOf(col, 1).includes('?'), true);
    expect(errorOf(col, 2), 'BLOB polymers are not supported');
    expect(errorOf(col, 3).includes('ambiguous'), true);
    expect(errorOf(col, 4).includes('range'), true);
    // Positions >= 100 (context/01 B3): the hwe splitter keeps the connection, so the 120-mer cycle builds
    const cyclic = seqOf(col, 5);
    const conns = cyclic.graphInfo.connections
      .map((c) => [c.seqIndex1, c.monomerIndex1, c.rGroup1, c.seqIndex2, c.monomerIndex2, c.rGroup2].join(','));
    expectArray(conns, ['0,119,2,0,0,1']);
    // HELM hwe cannot read (no monomers), or reads only partly: hwe's error
    expect(errorOf(col, 6).startsWith('HELM unexpected-char'), true, errorOf(col, 6));
    expect(errorOf(col, 7).startsWith('HELM unterminated-brace'), true, errorOf(col, 7));
  });

  // bio-edge 2026-11-02: an empty polymer gave the internal 'invalid chain start 0 at 1' (or '2 at 1' last)
  test('emptyPolymer', async () => {
    const col = seqCol([
      'PEPTIDE1{}|PEPTIDE2{A.C}$$$$V2.0',
      'PEPTIDE1{A.C}|PEPTIDE2{}$$$$V2.0',
      'PEPTIDE1{A}|RNA1{}|PEPTIDE2{C}$$$$V2.0',
      'PEPTIDE1{A}|PEPTIDE2{C}$$$$V2.0',
    ], NOTATION.HELM);
    expect(errorOf(col, 0), 'polymer PEPTIDE1 has no monomers');
    expect(errorOf(col, 1), 'polymer PEPTIDE2 has no monomers');
    expect(errorOf(col, 2), 'polymer RNA1 has no monomers');
    expectArray(seqOf(col, 3).graphInfo.disjointSeqStarts, [0, 1]);
  });

  // bio-edge 2026-11-02, pinned (not fixed): a BILN link number with one end (`C(1,3)-A-A`) or three
  // (`C(1,3)-A-C(1,3)-C(1,3)`) is dropped silently: the row builds without that bond. An error would change rows that
  // build today (the BILN splitter, which other agents tonight change too): listed for a decision.
  test('bilnUnmatchedLinks-pinned', async () => {
    const bilnHandler = {notation: NOTATION.BILN} as ISeqHandler;
    for (const [biln, conns] of [['C(1,3)-A-A', 0], ['C(1,3)-A-C(1,3)-C(1,3)', 1], ['C(1,3)-A.C(1,3).C(1,3)', 1]] as
      [string, number][]) {
      const res = toAsmSequence(splitterAsBiln(biln), bilnHandler, monomerLib);
      expect(res.error, null, biln);
      expect(res.seq!.graphInfo.connections.length, conns, biln);
    }
  });

  test('biln', async () => {
    // Samples of the BILN detector tests (detection needs realistic sequences)
    const col = DG.Column.fromList(DG.TYPE.STRING, 'seq', [
      'A-C(1,3)-G-G-H-A-V-E-A-K-Y-L-V-C(3,3)-S.G-I-V-E-A-C(2,3)-C(1,3)-T-S-I-C(2,3)-S-L-Y-Q-L-E-N-Y-C(3,3)-Y',
      'C-C(1,3)-S-W-P-A-R-C(2,3)-L-H-Q-D-L-C(3,3)-NH2.C(1,1)(2,2)(3,3)',
      'D-T-H-F-P-I-C(1,3)-I-F-C(2,3)-C(3,3)-G-C(2,3)-C(4,3)-H-R-S-K-C(3,3)-G-M-C(4,3)-C(1,3)-K-T',
    ]) as DG.Column<string>;
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
    if (semType) col.semType = semType;
    await seqHelper.getSeqHandler(col).refinerPromise;
    expect(col.meta.units, NOTATION.BILN);
    const s1 = seqOf(col, 1);
    expect(symbolsOf(s1).join('-'), 'C-C-S-W-P-A-R-C-L-H-Q-D-L-C-NH2-C');
    expectArray(s1.graphInfo.disjointSeqStarts, [0, 15]);
    const conns = s1.graphInfo.connections
      .map((c) => [c.seqIndex1, c.monomerIndex1, c.rGroup1, c.seqIndex2, c.monomerIndex2, c.rGroup2].join(','));
    expectArray(conns, ['0,13,3,1,0,3', '0,7,3,1,0,2', '0,1,3,1,0,1']);
    expect(seqOf(col, 0).graphInfo.connections.length, 3);
    expect(seqOf(col, 2).graphInfo.connections.length, 4);
  });

  test('buildSeqBatch', async () => {
    // Empty cells give empty sequences, one unknown monomer fails its row only (not the column)
    const col = seqCol(['PEPTIDE1{A.C}$$$$V2.0', '', 'PEPTIDE1{A.[Xyz]}$$$$V2.0', 'PEPTIDE1{G}$$$$V2.0'],
      NOTATION.HELM);
    const sh = seqHelper.getSeqHandler(col);
    const batch = buildSeqBatch(sh, null, monomerLib);
    expectArray(batch.rowError, [null, null, 'unknown monomer PEPTIDE:Xyz', null]);
    expectArray(symbolsOf(unpackRow(batch, 0)), ['A', 'C']);
    expect(unpackRow(batch, 1).monomers.length, 0);
    const sub = buildSeqBatch(sh, [3, 0], monomerLib);
    expectArray(symbolsOf(unpackRow(sub, 0)), ['G']);
  });

  test('libraryJson', async () => {
    const lib1 = getMonomerLibJson(monomerLib);
    const lib2 = getMonomerLibJson(monomerLib);
    expect(lib1.key, lib2.key);
    expect(lib1.json === lib2.json, true, 'JSON is cached per key');
    const monomers = JSON.parse(lib1.json) as any[];
    expect(monomers.length, getEffectiveMonomers(monomerLib).length);
    expect(monomers.some((m) => m.polymerType === 'PEPTIDE' && m.symbol === 'meI'), true);
    expect(monomers.every((m) => !('lib' in m) && !('wem' in m)), true);
  });

  test('libraryOverride', async () => {
    const a = monomerLib.getMonomer('PEPTIDE', 'A')!;
    const ovr = {...a, symbol: 'Ovr', name: 'override test', wem: undefined};
    const override = monomerLib.override({PEPTIDE: {Ovr: ovr}}, 'test');
    const col = seqCol(['PEPTIDE1{A.[Ovr]}$$$$V2.0'], NOTATION.HELM);
    col.temp[MmcrTemps.overriddenLibrary] = override;
    const lib = getSeqColMonomerLib(col, monomerLib);
    expect(getMonomerLibKey(lib) !== getMonomerLibKey(monomerLib), true, 'override gets its own key');
    const sh = seqHelper.getSeqHandler(col);
    const res = toAsmSequence(sh.getSplitted(0), sh, lib, {helm: col.get(0)!});
    expectArray(symbolsOf(res.seq!), ['A', 'Ovr']);
    expect(JSON.parse(getMonomerLibJson(lib).json).some((m: any) => m.symbol === 'Ovr'), true);
  });
});
