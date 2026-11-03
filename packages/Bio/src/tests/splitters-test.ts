import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {after, before, category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {TAGS as bioTAGS, splitterAsFasta} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {splitterAsHelm} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {ISeqGraphInfoHelm, ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {parseHelm} from '@datagrok-libraries/hwe';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {splitToMonomersUI} from '../utils/split-to-monomers';
import {awaitGrid} from './utils';
import * as C from '../utils/constants';

import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

category('splitters', async () => {
  before(async () => {
  });

  after(async () => {
  });

  const _helm1 = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$';

  const _helm2 = 'PEPTIDE1{meI.hHis.Hcy.Q.T.W.Q.Phe_4NH2.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.N.meK}$$$';

  const data: { [key: string]: [string, string[]] } = {
    fastaMulti: [
      'M[MeI]YKETLL[MeF]PKTDFPMRGGL[MeA]',
      ['M', 'MeI', 'Y', 'K', 'E', 'T', 'L', 'L', 'MeF', 'P',
        'K', 'T', 'D', 'F', 'P', 'M', 'R', 'G', 'G', 'L', 'MeA'],
    ],
    fastaFromHelm: [
      '[meI][Pip][dK][Thr_PO3H2][L-hArg(Et,Et)][D-Tyr_Et][Tyr_ab-dehydroMe][dV]EN[D-Orn][D-aThr][Phe_4Me]',
      ['meI', 'Pip', 'dK', 'Thr_PO3H2', 'L-hArg(Et,Et)', 'D-Tyr_Et', 'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn',
        'D-aThr', 'Phe_4Me'],
    ],

    helm1: [
      'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dE', 'Thr_PO3H2', 'Aca', 'D-Tyr_Et',
        'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn', 'D-aThr', 'Phe_4Me'],
    ],
    helm2: [
      'PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.D-Dap.dV.E.N.pnG.Phe_4Me}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dK', 'Thr_PO3H2', 'Aca',
        'D-Tyr_Et', 'D-Dap', 'dV', 'E', 'N', 'pnG', 'Phe_4Me'],
    ],
    // HELM editor dialog returns HELM string with multichar monomer names in square brackets
    helm3: [
      'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dK].[Thr_PO3H2].[Aca].[D-Tyr_Et].[D-Dap].[dV].E.N.[pnG].[Phe_4Me]}$$$',
      ['meI', 'hHis', 'Aca', 'N', 'T', 'dK', 'Thr_PO3H2', 'Aca',
        'D-Tyr_Et', 'D-Dap', 'dV', 'E', 'N', 'pnG', 'Phe_4Me'],
    ],

    // splitterAsHelm gives one position per hwe monomer: an RNA nucleotide `sugar(base)phosphate` is 3 positions,
    // `sugar(base)` 2, and a lone linker (`P`) one.
    testHelm1: [
      'RNA1{R(U)P.R(T)P.R(G)P.R(C)P.R(A)}$$$$',
      ['R', 'U', 'P', 'R', 'T', 'P', 'R', 'G', 'P', 'R', 'C', 'P', 'R', 'A'],
    ],

    testHelm2: [
      'RNA1{P.R(U)P.R(T)}$$$$',
      ['P', 'R', 'U', 'P', 'R', 'T'],
    ],
    testHelm3: [
      'RNA1{P.R(U).P.R(T)}$$$$',
      ['P', 'R', 'U', 'P', 'R', 'T'],
    ],
  };

  test('fastaMulti', async () => { await _testFastaSplitter(data.fastaMulti[0], data.fastaMulti[1]); });
  test('fastaFromHelm', async () => { await _testFastaSplitter(data.fastaFromHelm[0], data.fastaFromHelm[1]); });

  test('helm1', async () => { await _testHelmSplitter(data.helm1[0], data.helm1[1]); });
  test('helm2', async () => { await _testHelmSplitter(data.helm2[0], data.helm2[1]); });
  test('helm3-multichar', async () => { await _testHelmSplitter(data.helm3[0], data.helm3[1]); });

  // examples from Helm's former tests/test.csv file
  test('testHelm1', async () => { await _testHelmSplitter(data.testHelm1[0], data.testHelm1[1]); });
  test('testHelm2', async () => { await _testHelmSplitter(data.testHelm2[0], data.testHelm2[1]); });
  test('testHelm3', async () => { await _testHelmSplitter(data.testHelm3[0], data.testHelm3[1]); });
});

// splitterAsHelm on the hwe parser (Bio on Crux, PLAN §5.2): bugs of the regex splitter (context/01 §9)
category('splitters', () => {
  const originals = (ss: ISeqSplitted): string[] =>
    wu.count(0).take(ss.length).map((i) => ss.getOriginal(i)).toArray();
  const conns = (ss: ISeqSplitted): string[] => (ss.graphInfo?.connections ?? []).map((c) =>
    [c.seqIndex1, c.monomerIndex1, c.rGroup1, c.seqIndex2, c.monomerIndex2, c.rGroup2].join(','));
  const giOf = (ss: ISeqSplitted): ISeqGraphInfoHelm => ss.graphInfo as ISeqGraphInfoHelm;

  test('helm-B1-mixed-types', async () => {
    // The chain index comes from the chain's order in the string, not from the number in its id
    const ss = splitterAsHelm('PEPTIDE1{A.C.D}|CHEM1{[PEG2]}$PEPTIDE1,CHEM1,3:R2-1:R1$$$V2.0');
    expectArray(originals(ss), ['A', 'C', 'D', 'PEG2']);
    expectArray(conns(ss), ['0,2,2,1,0,1']);
    expectArray(ss.graphInfo!.polymerTypes!, ['PEPTIDE', 'PEPTIDE', 'PEPTIDE', 'CHEM']);
    expectArray(giOf(ss).chainIds!, ['PEPTIDE1', 'CHEM1']);
  });

  test('helm-B1-chem-first', async () => {
    // HELM_CHEMS.csv row 10 shape: the legacy splitter dropped this connection
    const ss = splitterAsHelm('CHEM1{[SMCC]}|PEPTIDE1{G.G.C}$PEPTIDE1,CHEM1,3:R3-1:R1$$$V2.0');
    expectArray(originals(ss), ['SMCC', 'G', 'G', 'C']);
    expectArray(ss.graphInfo!.disjointSeqStarts, [0, 1]);
    expectArray(conns(ss), ['1,2,3,0,0,1']);
    expect(ss.graphInfo!.connections[0].seq1Type, 'PEPTIDE');
    expect(ss.graphInfo!.connections[0].seq2Type, 'CHEM');
  });

  test('helm-B2-annotated-chains', async () => {
    // A polymer annotation no longer merges two chains
    const ss = splitterAsHelm('PEPTIDE1{A.C.D}"HC"|PEPTIDE2{G.C}$$$$V2.0');
    expectArray(originals(ss), ['A', 'C', 'D', 'G', 'C']);
    expectArray(ss.graphInfo!.disjointSeqStarts, [0, 3]);
    expectArray(ss.graphInfo!.polymerTypes!, ['PEPTIDE', 'PEPTIDE', 'PEPTIDE', 'PEPTIDE', 'PEPTIDE']);
  });

  test('helm-B3-position-100', async () => {
    // A 120-mer head-to-tail cycle keeps its connection
    const ss = splitterAsHelm(`PEPTIDE1{${new Array(120).fill('A').join('.')}}$PEPTIDE1,PEPTIDE1,120:R2-1:R1$$$V2.0`);
    expect(ss.length, 120);
    expectArray(conns(ss), ['0,119,2,0,0,1']);
  });

  test('helm-B3-rgroup-10', async () => {
    const ss = splitterAsHelm('PEPTIDE1{A.C}$PEPTIDE1,PEPTIDE1,1:R10-2:R3$$$V2.0');
    expectArray(conns(ss), ['0,0,10,0,1,3']);
  });

  test('helm-inline-smiles-dot', async () => {
    // Inline SMILES with '.' or nested brackets is one monomer
    const ss = splitterAsHelm('CHEM1{[C(=O)O.N[*:1]]}|PEPTIDE1{A.G.C}$PEPTIDE1,CHEM1,3:R3-1:R1$$$V2.0');
    expectArray(originals(ss), ['C(=O)O.N[*:1]', 'A', 'G', 'C']);
    expectArray(conns(ss), ['1,2,3,0,0,1']);
  });

  test('helm-inline-smiles-rna-base', async () => {
    const ss = splitterAsHelm('RNA1{R([NC1=NC(N([*:1])C=C1)=O])P.R(A)}$$$$');
    expectArray(originals(ss), ['R', 'NC1=NC(N([*:1])C=C1)=O', 'P', 'R', 'A']);
    expectArray(giOf(ss).biotypes!, ['HELM_SUGAR', 'HELM_BASE', 'HELM_LINKER', 'HELM_SUGAR', 'HELM_BASE']);
  });

  test('helm-monomer-annotation', async () => {
    const ss = splitterAsHelm('PEPTIDE1{A.A.C"mutation".D}$$$$V2.0');
    expectArray(originals(ss), ['A', 'A', 'C', 'D']);
  });

  test('helm-hbonds-out-of-connections', async () => {
    // Old consumers (joinToHelm, joinToBiln, the POM choice) read connections as covalent bonds
    const ss = splitterAsHelm('RNA1{R(A)P.R(C)}|RNA2{R(G)P.R(U)}$RNA1,RNA2,2:pair-5:pair$$$V2.0');
    expect(ss.graphInfo!.connections.length, 0);
    const hb = giOf(ss).hbonds!;
    expect(hb.length, 1);
    expect(hb[0].kind, 'hbond');
    expectArray([hb[0].seqIndex1, hb[0].monomerIndex1, hb[0].seqIndex2, hb[0].monomerIndex2], [0, 1, 1, 4]);
  });

  const repeatsOf = (ss: ISeqSplitted): string[] =>
    (giOf(ss).repeats ?? []).map((r) => [r.start, r.end, r.min, r.max].join(','));
  const errorsOf = (helm: string): string[] => giOf(splitterAsHelm(helm)).errors ?? [];

  test('helm-repeats', async () => {
    // Listed once, as written (position = hwe continuousId - 1): renderers and WebLogo see the source tokens.
    // graphInfo.repeats holds the fixed repeats for the assembly input, which expands them.
    const ss = splitterAsHelm('PEPTIDE1{A.(G.A.C)\'3\'.A}$$$$V2.0');
    expectArray(originals(ss), ['A', 'G', 'A', 'C', 'A']);
    expectArray(repeatsOf(ss), ['1,4,3,3']);
    expect(giOf(ss).errors!.length, 0);
    // One monomer, and the end of a chain followed by another chain
    const one = splitterAsHelm('PEPTIDE1{A.G.A.C\'5\'}|PEPTIDE2{C.C}$PEPTIDE1,PEPTIDE2,4:R3-1:R3$$$V2.0');
    expectArray(repeatsOf(one), ['3,4,5,5']);
    expectArray(conns(one), ['0,3,3,1,0,3']);
  });

  test('helm-repeats-rna', async () => {
    // An RNA nucleotide repeat covers sugar, base and linker; a sugar with its base alone takes the base in
    const ss = splitterAsHelm('RNA1{R(C)P.(R(U)P)\'3\'.R(A)P}$$$$V2.0');
    expectArray(originals(ss), ['R', 'C', 'P', 'R', 'U', 'P', 'R', 'A', 'P']);
    expectArray(repeatsOf(ss), ['3,6,3,3']);
    expectArray(repeatsOf(splitterAsHelm('RNA1{(R(U))\'2\'.R(A)P}$$$$V2.0')), ['0,2,2,2']);
  });

  test('helm-repeats-errors', async () => {
    // D9: ranges, CHEM repeats, nested repeats and zero counts are per-row errors, not repeats
    const range = splitterAsHelm('RNA1{R(A)P.(R(N)P)\'4\'.(R(G)P)\'3-7\'}$$$$V2.0');
    expectArray(repeatsOf(range), ['3,6,4,4']);
    expectArray(giOf(range).errors!, ['repeat range \'3-7\' at position 7 is not supported']);
    expect(errorsOf('CHEM1{[Az]\'3\'}$$$$V2.0').some((e) => e.includes('repeat \'3\'')), true);
    expect(errorsOf('PEPTIDE1{A.((G.A)\'2\'.C)\'3\'.A}$$$$V2.0').some((e) => e.startsWith('nested repeats')), true);
    expect(errorsOf('PEPTIDE1{A\'0\'.C}$$$$V2.0').some((e) => e.includes('repeat-zero-count')), true);
  });

  test('helm-ambiguous-connections', async () => {
    // Multi-site positions, ranges, '?' and group endpoints: hwe draws them at the first candidate; the adapter
    // leaves them out of connections (old consumers would take them for definite bonds) and errors the row
    for (const [helm, text] of [
      ['PEPTIDE1{A.C.K.D}|CHEM1{[sDBL]}$PEPTIDE1,CHEM1,(1,3):R3-1:R1$$$V2.0', 'PEPTIDE1,CHEM1,(1,3):R3-1:R1'],
      ['PEPTIDE1{A.C.K.D}|CHEM1{[sDBL]}$PEPTIDE1,CHEM1,2-4:R3-1:R1$$$V2.0', 'PEPTIDE1,CHEM1,2-4:R3-1:R1'],
      ['PEPTIDE1{A.C}$PEPTIDE1,PEPTIDE1,?:R1-2:R2$$$V2.0', 'PEPTIDE1,PEPTIDE1,?:R1-2:R2'],
      ['BLOB1{Bead}|PEPTIDE1{A.K}$G1,PEPTIDE1,K:R3-1:R1$G1(BLOB1+PEPTIDE1)$$V2.0', 'G1,PEPTIDE1,K:R3-1:R1'],
    ]) {
      const ss = splitterAsHelm(helm);
      expect(ss.graphInfo!.connections.length, 0, helm);
      expectArray(giOf(ss).errors!, [`ambiguous connection '${text}' is not supported`]);
    }
    // A definite connection next to an ambiguous one stays
    const mixed = splitterAsHelm(
      'PEPTIDE1{A.C}|PEPTIDE2{G}$PEPTIDE1,PEPTIDE2,(1,2):R3-1:R1|PEPTIDE1,PEPTIDE2,2:R2-1:R1$$$V2.0');
    expectArray(conns(mixed), ['0,1,2,1,0,1']);
    // A symbol position whose symbol occurs once is one monomer: a connection, no error
    const sym = splitterAsHelm('PEPTIDE1{A.A.C.D}|CHEM1{[PEG2]}$PEPTIDE1,CHEM1,C:R3-1:R1$$$V2.0');
    expectArray(conns(sym), ['0,2,3,1,0,1']);
    expect(giOf(sym).errors!.length, 0);
    expect(errorsOf('PEPTIDE1{C.A.C}|CHEM1{[PEG2]}$PEPTIDE1,CHEM1,C:R3-1:R1$$$V2.0')[0]
      .startsWith('ambiguous connection'), true);
  });

  test('helm-errors-for-assembly', async () => {
    expect(errorsOf('PEPTIDE1{A.(G,C).A}$$$$V2.0')[0], 'ambiguous monomer \'(G,C)\' at position 2 is not supported');
    expect(errorsOf('PEPTIDE1{A.C}$PEPTIDE1,PEPTIDE1,?:R1-2:R2$$$V2.0')[0].startsWith('ambiguous connection'), true);
    expect(errorsOf('PEPTIDE1{A.C.D}$$$$V2.0').length, 0);
    // Two polymers with one id: two chains, and an error (hwe sends connections to the last one)
    const dup = splitterAsHelm('PEPTIDE1{A.C}|PEPTIDE1{G.G}$$$$V2.0');
    expectArray(dup.graphInfo!.disjointSeqStarts, [0, 2]);
    expect(giOf(dup).errors![0].includes('duplicate-polymer-id'), true);
  });

  test('helm-whitespace', async () => {
    // context/01 B10: hwe 1.1.0 skips whitespace outside brackets (a warning), so these rows parse in hwe
    for (const helm of ['PEPTIDE1{A.C}$$$$V2.0 ', 'PEPTIDE1{A . C}$$$$V2.0', ' PEPTIDE1{A.C}\t$$$$V2.0']) {
      const ss = splitterAsHelm(helm);
      expectArray(originals(ss), ['A', 'C']);
      expect(giOf(ss).errors!.length, 0);
    }
  });

  test('helm-blob-with-spaces', async () => {
    // hwe 1.0.5 failed the whole parse on the space; the regex splitter glued the rest of the string to the BLOB
    const ss = splitterAsHelm('BLOB1{Gold Particle}$$$$V2.0');
    expectArray(originals(ss), ['Gold Particle']);
    expectArray(ss.graphInfo!.polymerTypes!, ['BLOB']);
    expect(giOf(ss).errors!.length, 0);
    const conj = splitterAsHelm('PEPTIDE1{C.C}|BLOB1{Gold Particle}"Au10, Diameter:10nm"$PEPTIDE1,BLOB1,C:R3-?:?' +
      '$G1(PEPTIDE1:20-34+BLOB1)${"Name":"Gold particle conjugated with peptides","Load":26}$V2.0');
    expectArray(originals(conj), ['C', 'C', 'Gold Particle']);
    expectArray(conj.graphInfo!.disjointSeqStarts, [0, 2]);
    expectArray(giOf(conj).errors!, ['ambiguous connection \'PEPTIDE1,BLOB1,C:R3-?:?\' is not supported']);
  });

  test('helm-position-invariant', async () => {
    // Hover and PolyTool map hwe atoms to positions by continuousId - 1 (context/01 §12)
    for (const helm of [
      'RNA1{R(U)P.R(T)P}|RNA2{P.R(A)P.R(A)}$RNA1,RNA2,2:pair-6:pair$$$',
      'CHEM1{[SMCC]}|PEPTIDE1{G.[meI].C}$PEPTIDE1,CHEM1,3:R3-1:R1$$$V2.0',
      'RNA1{[dR](A)p.r([NC1=NC(N([*:1])C=C1)=O])p.r(G)}$$$$V2.0',
      'RNA1{R(C)P.(R(U)P)\'3\'.R(A)P}$$$$V2.0',
      'PEPTIDE1{A . C}|PEPTIDE2{G}$$$$V2.0 ',
    ]) {
      const ss = splitterAsHelm(helm);
      const atoms = parseHelm(helm, null, {mode: 'lenient'}).mol.atoms;
      expect(ss.length, atoms.length);
      for (const a of atoms)
        expect(ss.getOriginal(a.bio.continuousId - 1), a.symbol);
    }
  });

  test('helm-no-sections', async () => {
    // The regex splitter returned no monomers without the '$' sections
    expectArray(originals(splitterAsHelm('PEPTIDE1{A.C}')), ['A', 'C']);
  });

  test('helm-unparsable', async () => {
    // A string hwe's lexer rejects (an unexpected character): no monomers, hwe's error, unbuildable at atomic level
    const ss = splitterAsHelm('PEPTIDE1{A.C.D}$$$$V2.0#');
    expect(ss.length, 0);
    expect(giOf(ss).recovered, true);
    expect(giOf(ss).errors![0].startsWith('HELM unexpected-char'), true);
  });

  test('helm-recovered', async () => {
    // Invalid HELM (a peptide branch): hwe recovers after 'Ac' and skips the rest of the polymer. The row keeps the
    // monomers hwe read, with hwe's error, and stays unbuildable at atomic level.
    const ss = splitterAsHelm('PEPTIDE1{Ac(1).A.A.C(1).G}$$$$');
    expectArray(originals(ss), ['Ac']);
    expect(giOf(ss).recovered, true);
    expect(giOf(ss).errors![0].startsWith('HELM unterminated-brace'), true);
    // A recovery that loses no monomer
    const rep = splitterAsHelm('CHEM1{[Az]\'3\'}$$$$V2.0');
    expectArray(originals(rep), ['Az']);
    expect(giOf(rep).recovered, true);
    expect(giOf(rep).errors!.some((e) => e.startsWith('HELM repeat-on-non-iterable-polymer')), true);
  });
});

category('splitters', () => {
  let seqHelper: ISeqHelper;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  test('splitToMonomers', async () => {
    const df: DG.DataFrame = await grok.dapi.files.readCsv('System:AppData/Bio/samples/MSA.csv');

    const seqCol = df.getCol('MSA');
    const semType = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
    if (semType)
      seqCol.semType = semType;
    seqCol.setTag(bioTAGS.aligned, C.MSA);

    const newDf: DG.DataFrame = await grok.functions.call('Bio:splitToMonomersTopMenu', {table: df, sequence: seqCol});
    expect(newDf.columns.names().includes('17'), true);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(newDf);

    // TODO: Check cell.renderer for columns of monomers
    const tv: DG.TableView = grok.shell.addTableView(newDf);
    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);
  });

  test('getHelmMonomers', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(
      `HELM,Activity
PEPTIDE1{hHis.N.T}$$$,5.30751
PEPTIDE1{hHis.Aca.Cys_SEt}$$$,5.72388
`);
    await grok.data.detectSemanticTypes(df);
    const expectedMonomerList = ['hHis', 'Aca', 'Cys_SEt', 'N', 'T'];

    const helmCol: DG.Column = df.getCol('HELM');
    const res = await grok.functions.call('Bio:getHelmMonomers', {sequence: helmCol}) as string[];

    const missed = expectedMonomerList.filter((m) => !res.includes(m));
    const unexpected = res.filter((m) => !expectedMonomerList.includes(m));
    if (missed.length > 0 || unexpected.length) {
      const msgs = [];
      if (missed.length > 0)
        msgs.push(`Missed monomers ${JSON.stringify(missed)}.`);
      if (unexpected.length > 0)
        msgs.push(`Unexpected monomers ${JSON.stringify(unexpected)}.`);

      throw new Error(msgs.join(' '));
    }
  });

  // test('helmAsFasta', async () => {
  //   // The columns can't be empty for SeqHandler
  //   /* eslint-disable max-len */
  //   const srcSeq = '[meI][Pip][dK][Thr_PO3H2][L-hArg(Et,Et)][D-Tyr_Et][Tyr_ab-dehydroMe][dV]EN[D-Orn][D-aThr][Phe_4Me]';
  //   const tgtSeqA = ['meI', 'Pip', 'dK', 'Thr_PO3H2', 'L-hArg(Et,Et)', 'D-Tyr_Et', 'Tyr_ab-dehydroMe', 'dV', 'E', 'N', 'D-Orn', 'D-aThr', 'Phe_4Me'];
  //   /* eslint-enable max-len */
  //   const resSeqA = splitterAsFasta(srcSeq);
  //   expectArray(resSeqA, tgtSeqA);
  // });
});

export async function _testFastaSplitter(src: string, tgt: string[]) {
  const resSS: ISeqSplitted = splitterAsFasta(src);
  console.debug(`Bio: tests: splitters: src=${JSON.stringify(src)}, res=${JSON.stringify(resSS)} .`);
  expectArray(wu.count(0).take(resSS.length).map((p) => resSS.getOriginal(p)).toArray(), tgt);
}

export async function _testHelmSplitter(src: string, tgt: string[]) {
  const resSS: ISeqSplitted = splitterAsHelm(src);
  console.debug(`Bio: tests: splitters: src=${JSON.stringify(src)}, res=${JSON.stringify(resSS)} .`);
  expectArray(wu.count(0).take(resSS.length).map((p) => resSS.getOriginal(p)).toArray(), tgt);
}
