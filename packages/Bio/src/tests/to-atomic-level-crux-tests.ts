/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {getMonomerLibHelper, IMonomerLib, IMonomerLibBase, IMonomerLibHelper, Monomer}
  from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {assembleRowWithMap, AsmLayoutName, SEQUENCE_SRC_LAYOUT_TAG}
  from '@datagrok-libraries/bio/src/monomer-works/assembler';

import {sequenceToMolfile} from '../utils/sequence-to-mol';
import * as mg from './utils/molblock-graph';

/** Reproducing tests of the bugs the assembler (Bio on Crux, PLAN §5.4) fixes: connections the old splitter or
 * converter dropped, per-row errors instead of whole-column failures, stereo the linear path lost, double bonds,
 * chains, caps. The chemistry is read from the molblocks (./utils/molblock-graph) or compared as canonical SMILES
 * through Chem. */
category('toAtomicLevelCrux', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let seqHelper: ISeqHelper;
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

  function seqCol(list: string[], units: NOTATION, alphabet?: ALPHABET): DG.Column<string> {
    const col = DG.Column.fromList(DG.TYPE.STRING, 'seq', list) as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = units;
    col.setTag(bioTAGS.aligned, ALIGNMENT.SEQ);
    if (alphabet) col.setTag(bioTAGS.alphabet, alphabet);
    return col;
  }

  const helmCol = (list: string[]) => seqCol(list, NOTATION.HELM);

  /** Converts the column the way To Atomic Level does; one molblock per row ('' for a failed row). */
  async function convert(col: DG.Column<string>, nonlinear: boolean, lib?: IMonomerLibBase):
    Promise<{ mols: string[], warnings: string[], molCol: DG.Column<string> }> {
    const df = DG.DataFrame.fromColumns([col]);
    const res = await sequenceToMolfile(df, col, nonlinear, false, lib ?? monomerLib, seqHelper);
    expect(res.molCol != null, true, `no molecule column: ${res.warnings.join('; ')}`);
    const molCol = res.molCol!;
    const mols = Array.from({length: col.length}, (_, i) => molCol.get(i) ?? '');
    return {mols, warnings: res.warnings, molCol};
  }

  /** One row that must build, as a graph. */
  async function buildOne(seq: string, nonlinear: boolean, lib?: IMonomerLibBase): Promise<mg.MolblockGraph> {
    const {mols, warnings} = await convert(helmCol([seq]), nonlinear, lib);
    expect(warnings.length, 0, `unexpected row errors: ${warnings.join('; ')}`);
    expect(!!mols[0], true, `no molecule for ${seq}`);
    return mg.parseMolblockGraph(mols[0]);
  }

  /** Atoms with specified stereo in a canonical SMILES */
  const stereoCentres = (smiles: string): number => (smiles.match(/\[[^\]]*@[^\]]*\]/g) ?? []).length;

  const warningOfRow = (warnings: string[], rowIdx: number): string | undefined =>
    warnings.find((w) => w.startsWith(`Cannot get molfile of row #${rowIdx}:`));

  // B3: connections at positions >= 100 were dropped by the old splitter (the 100-mer came out linear)
  test('connections-100', async () => {
    const helm = 'PEPTIDE1{' + new Array(100).fill('G').join('.') + '}$PEPTIDE1,PEPTIDE1,100:R2-1:R1$$$V2.0';
    for (const nonlinear of [false, true]) { // linear falls back to the non-linear layout for the cycle
      const g = await buildOne(helm, nonlinear);
      // cyclo-(Gly)100: 4 heavy atoms per residue, one ring (a linear chain would have 401 atoms, 400 bonds)
      expect(mg.heavyAtomCount(g), 400, 'heavy atoms');
      expect(mg.heavyBondCount(g), 400, 'heavy bonds');
      expect(mg.fragmentCount(g), 1, 'fragments');
    }
  });

  // B1: polymer ids out of numeric order and mixed polymer types (the old splitter pointed CHEM1 at chain 0)
  test('mixed-type-connection', async () => {
    const helm = 'CHEM1{[*N[C@H](C(=O)*)Cc1ccc(cc1)OP(=O)(O)O |$_R1;;;;;_R2;;;;;;;;;;;;$|]}|PEPTIDE1{[C[C@H](N[*])C(=O)C[*] |$;;;_R1;;;;_R2$|].G.G.G.C.C.K.K.K.K}$PEPTIDE1,CHEM1,10:R3-1:R1$$$V2.0';
    const g = await buildOne(helm, true);
    expect(mg.fragmentCount(g), 1, 'the PEPTIDE1-CHEM1 connection must join the chains');
    expect(mg.countElement(g, 'P'), 1, 'the CHEM monomer\'s phosphate');
  });

  // B2: an annotated chain was merged with the next one or lost
  test('annotated-chain', async () => {
    const g = await buildOne('PEPTIDE1{A.C.G}"chain annotation"$$$$V2.0', false);
    expect(mg.fragmentCount(g), 1, 'fragments');
    const {mols} = await convert(helmCol(['PEPTIDE1{A.C.G}"chain annotation"$$$$V2.0']), false);
    expect(await mg.canonicalSmiles(mols[0]), await mg.canonicalSmiles('C[C@H](N)C(=O)N[C@@H](CS)C(=O)NCC(=O)O'));
  });

  // siRNA duplex: hydrogen bonds ('pair') add no covalent bond; the old non-linear converter threw on them
  test('sirna-pair', async () => {
    const helm = 'RNA1{r(A)p.r(C)}|RNA2{r(G)p.r(U)}$RNA1,RNA2,2:pair-5:pair$$$V2.0';
    for (const nonlinear of [false, true]) {
      const g = await buildOne(helm, nonlinear);
      expect(mg.fragmentCount(g), 2, 'two strands, no bond between them');
      expect(mg.furanoseRings(g).length, 4, 'four riboses');
      expect(mg.countElement(g, 'P'), 2, 'one linker per strand');
    }
  });

  // '?' and repeated-symbol endpoints: an error for that row only
  test('unbuildable-endpoints-per-row', async () => {
    const {mols, warnings} = await convert(helmCol([
      'PEPTIDE1{A.C}$PEPTIDE1,PEPTIDE1,?:R1-2:R2$$$V2.0',
      'PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,C:R3-3:R3$$$V2.0',
      'PEPTIDE1{A.C}$$$$V2.0',
    ]), true);
    expect(mols[0], '', 'row 0 (?)');
    expect(mols[1], '', 'row 1 (C:R3)');
    expect(warningOfRow(warnings, 0) != null && warningOfRow(warnings, 1) != null, true, warnings.join('; '));
    expect(!!mols[2] && warningOfRow(warnings, 2) == null, true, 'row 2 builds');
  });

  const cysLike = (symbol: string, sideR: number, sideCap: string, sideCapName: string): Monomer => ({
    symbol: symbol, name: `${symbol} (test monomer)`, molfile: '',
    smiles: `O=C([C@H](CS[*:${sideR}])N[*:1])[*:2]`,
    polymerType: 'PEPTIDE', monomerType: 'Backbone', naturalAnalog: 'C', id: 0, author: 'Test', createDate: null,
    rgroups: [
      {alternateId: 'R1-H', capGroupName: 'H', capGroupSMILES: '[*:1][H]', label: 'R1'},
      {alternateId: 'R2-OH', capGroupName: 'OH', capGroupSMILES: 'O[*:2]', label: 'R2'},
      {alternateId: `R${sideR}-${sideCapName}`, capGroupName: sideCapName, capGroupSMILES: sideCap, label: `R${sideR}`},
    ],
    meta: {},
  } as Monomer);

  const libWith = (...monomers: Monomer[]): IMonomerLibBase => monomerLib.override(
    {'PEPTIDE': Object.fromEntries(monomers.map((m) => [m.symbol, m]))}, 'toAtomicLevelCrux test monomers');

  // R10+ attachment points (the old converter read one digit of the label: R10 -> R1)
  test('r10-connection', async () => {
    const lib = libWith(cysLike('Cys_R10', 10, '[*:10][H]', 'H'));
    const g = await buildOne('PEPTIDE1{C.[Cys_R10].A}$PEPTIDE1,PEPTIDE1,1:R3-2:R10$$$V2.0', true, lib);
    expect(mg.bondsBetween(g, 'S', 'S'), 1, 'the R3-R10 disulfide');
    expect(mg.fragmentCount(g), 1, 'fragments');
    expect(mg.heavyBondCount(g), mg.heavyAtomCount(g), 'one ring');
  });

  // 5'-terminal monomers with R1 only (Bio, DBCO): the old non-linear converter could not join them
  test('rna-5p-terminal-r1-only-nonlinear', async () => {
    const g = await buildOne('RNA1{[Bio].r(T)}$$$$V2.0', true);
    expect(mg.fragmentCount(g), 1, 'Bio joined to the chain');
    expect(mg.cyclicUreaRings(g) >= 1, true, 'biotin\'s ureido ring');
    expect(mg.furanoseRings(g).length, 1, 'the r(T) ribose');
  });

  // An 'X' cap: the old converter wrote a literal X atom into the molfile
  test('x-cap', async () => {
    const lib = libWith(cysLike('Cys_XCap', 3, 'X[*:3]', 'X'));
    const joined = await buildOne('PEPTIDE1{[Cys_XCap].C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0', true, lib);
    expect(mg.countElement(joined, 'X'), 0, 'no X atom');
    expect(mg.bondsBetween(joined, 'S', 'S'), 1, 'the R3-R3 disulfide');
    // Left unjoined, the X-capped R-group has no cap: a row error (never a literal X atom)
    const {mols, warnings} = await convert(helmCol(['PEPTIDE1{[Cys_XCap]}$$$$V2.0']), true, lib);
    expect(mols[0], '', 'no molecule when the X-capped R-group is unjoined');
    expect(warningOfRow(warnings, 0) != null, true, `a row error: ${warnings.join('; ')}`);
  });

  // B1 of the linear path: the sugars' C1' and C3' lost their stereo (only C4' was kept)
  test('sugar-stereo', async () => {
    const cases: [string, ALPHABET, number][] = [
      ['A', ALPHABET.DNA, 3], ['AC', ALPHABET.DNA, 6], ['A', ALPHABET.RNA, 4],
    ];
    for (const [seq, alphabet, centres] of cases) {
      const {mols, warnings} = await convert(seqCol([seq], NOTATION.FASTA, alphabet), false);
      expect(warnings.length, 0, warnings.join('; '));
      expect(stereoCentres(await mg.canonicalSmiles(mols[0])), centres, `${alphabet} ${seq}: specified sugar centres`);
    }
  });

  // B2 of the linear path: wedges moved to the wrong bonds, meI and meT lost their Cbeta centre
  test('wedge-shift-monomers', async () => {
    for (const symbol of ['meI', 'meT']) {
      const {mols, warnings} = await convert(helmCol([`PEPTIDE1{[${symbol}]}$$$$V2.0`]), false);
      expect(warnings.length, 0, warnings.join('; '));
      expect(stereoCentres(await mg.canonicalSmiles(mols[0])), 2, `${symbol}: Calpha and Cbeta`);
    }
  });

  // The old linear path lost the non-aromatic double bonds of dehydro/alkene monomers (a wrong formula on 513/540
  // rows of HELM.csv)
  test('dehydro-double-bonds', async () => {
    const symbols = ['Dha', 'Tyr_ab-dehydroMe', 'Phe_ab-dehydro', 'Bmt'];
    let expected = 0;
    for (const symbol of symbols) {
      const m = monomerLib.getMonomer('PEPTIDE', symbol);
      expect(!!m?.molfile, true, `${symbol} in the library`);
      expected += mg.chainDoubleBonds(mg.parseMolblockGraph(m!.molfile));
    }
    for (const nonlinear of [false, true]) {
      const g = await buildOne(`PEPTIDE1{${symbols.map((s) => `[${s}]`).join('.')}}$$$$V2.0`, nonlinear);
      expect(mg.chainDoubleBonds(g), expected, `double bonds outside rings (nonlinear=${nonlinear})`);
    }
  });

  /** The row's monomer map keys and atom lists, through the hover's call */
  async function mapOf(col: DG.Column<string>, layout: AsmLayoutName) {
    const res = await assembleRowWithMap(col, 0, monomerLib, seqHelper, layout);
    expect(res.error, null, `row error: ${res.error}`);
    const atomCount = mg.parseMolblockGraph(res.molfile).elements.length;
    const keys = [...res.monomers.keys()].sort((a, b) => a - b);
    return {keys, monomers: res.monomers, atomCount};
  }

  // The gap map: gaps own no atoms and the other monomers keep their positions (gaps counted)
  test('gap-map', async () => {
    const fasta = await mapOf(seqCol(['MD-YK'], NOTATION.FASTA, ALPHABET.PT), 'linear');
    expectArray(fasta.keys, [0, 1, 3, 4]);
    const seen = new Set<number>();
    for (const v of fasta.monomers.values()) {
      for (const a of v.atoms) {
        expect(seen.has(a), false, `atom ${a} in two monomers`);
        expect(a < fasta.atomCount, true, `atom ${a} out of the molblock`);
        seen.add(a);
      }
    }
    const helm = await mapOf(helmCol(['PEPTIDE1{A.*.C}$$$$V2.0']), 'nonlinear');
    expectArray(helm.keys, [0, 2]);
  });

  // B20: multi-chain peptides in linear were fused into one chain (and the hover map disagreed)
  test('multi-chain-peptide-linear', async () => {
    const col = helmCol(['PEPTIDE1{A.C}|PEPTIDE2{G.K}$$$$V2.0']);
    const {mols, warnings, molCol} = await convert(col, false);
    expect(warnings.length, 0, warnings.join('; '));
    const g = mg.parseMolblockGraph(mols[0]);
    expect(mg.fragmentCount(g), 2, 'two chains, not fused');
    expect(mg.heavyAtomCount(g), 26, 'Ala-Cys (12) + Gly-Lys (14)');
    expect(molCol.getTag(SEQUENCE_SRC_LAYOUT_TAG), 'linear');
    const map = await mapOf(col, 'linear');
    expectArray(map.keys, [0, 1, 2, 3]);
  });

  // B5: one unknown monomer failed the whole column
  test('unknown-monomer-per-row', async () => {
    const {mols, warnings} = await convert(helmCol(['PEPTIDE1{A.[Xyz].C}$$$$V2.0', 'PEPTIDE1{A.C}$$$$V2.0']), false);
    expect(mols[0], '', 'row 0');
    expect(warningOfRow(warnings, 0)?.includes('Xyz') ?? false, true, warnings.join('; '));
    expect(!!mols[1], true, 'row 1 builds');
  });

  // BLOB, ambiguity and repeat ranges: an error for that row only
  test('blob-ambiguity-range-per-row', async () => {
    const {mols, warnings} = await convert(helmCol([
      'BLOB1{Bead}|PEPTIDE1{A}$$$$V2.0',
      'PEPTIDE1{A.(G,C).A}$$$$V2.0',
      'PEPTIDE1{A\'2-5\'.G}$$$$V2.0',
      'PEPTIDE1{A.G}$$$$V2.0',
    ]), true);
    for (let rowIdx = 0; rowIdx < 3; ++rowIdx) {
      expect(mols[rowIdx], '', `row ${rowIdx}`);
      expect(warningOfRow(warnings, rowIdx) != null, true, `row ${rowIdx} warning: ${warnings.join('; ')}`);
    }
    expect(!!mols[3], true, 'row 3 builds');
  });

  // context/04 surprise 5: the non-linear conversion of a non-HELM column renamed the user's column and added a
  // temporary HELM column during the call
  test('no-temporary-columns', async () => {
    const col = seqCol(['MDYKETLLMPK'], NOTATION.FASTA, ALPHABET.PT);
    const df = DG.DataFrame.fromColumns([col]);
    const names: string[][] = [];
    const sub = df.onColumnsAdded.subscribe(() => names.push(df.columns.names()));
    try {
      await sequenceToMolfile(df, col, true, false, monomerLib, seqHelper);
    } finally {
      sub.unsubscribe();
    }
    expectArray(df.columns.names(), ['seq', 'molfile(seq)']);
    expect(names.every((n) => n[0] === 'seq'), true, `the sequence column was renamed: ${JSON.stringify(names)}`);
  });

  // helmToAtomicLevel ignored `highlight` (PolyTool's "Highlight monomers" did nothing)
  test('helm-to-atomic-level-highlight', async () => {
    for (const highlight of [true, false]) {
      const res = await seqHelper.helmToAtomicLevel(helmCol(['PEPTIDE1{A.C.G}$$$$V2.0']), false, highlight);
      expect(res.molCol?.getTag(ChemTags.SEQUENCE_SRC_HL_MONOMERS), String(highlight));
      expect(res.molCol?.getTag(SEQUENCE_SRC_LAYOUT_TAG), 'nonlinear');
    }
  });

  // The legacy synchronous converter (D14) runs on the assembler's synchronous twins after its async factory readied
  // them: the same molecules as the workers' column
  test('sync-converter', async () => {
    const helms = ['PEPTIDE1{A.C.G}$$$$V2.0', '', 'PEPTIDE1{A.[Xyz]}$$$$V2.0'];
    const converter = await seqHelper.getHelmToMolfileConverter(monomerLib);
    const withMaps = converter.convertToMolfileV3K(helms);
    expect(withMaps[1].molfile, '', 'an empty HELM');
    expect(withMaps[2].molfile, '', 'an unknown monomer');
    expectArray([...withMaps[0].monomers.keys()].sort((a, b) => a - b), [0, 1, 2]);
    expect(await mg.canonicalSmiles(withMaps[0].molfile),
      await mg.canonicalSmiles('C[C@H](N)C(=O)N[C@@H](CS)C(=O)NCC(=O)O'), 'the molecule');
    expect(seqHelper.helmToAtomicLevelSingle(helms[0], converter).molfile, withMaps[0].molfile,
      'helmToAtomicLevelSingle');
    const syncCol = converter.convertToMolfileV3KColumn(helmCol(helms));
    const workers = await seqHelper.helmToAtomicLevel(helmCol(helms), false, false);
    expect(syncCol.get(0), workers.molCol!.get(0), 'the main thread\'s column as the workers\'');
    expect(syncCol.get(2), '', 'the column\'s failed row');
  });

  // bio-edge 2026-11-02: HELM to Molecule assembled the whole column on the main thread (the legacy synchronous
  // converter): the page froze on large columns, and a row too large for memory switched Crux off for the session.
  // It runs on the workers now: the same column as the synchronous converter's
  test('helm-to-molecule-workers', async () => {
    const col = helmCol(['PEPTIDE1{A.C.G}$$$$V2.0', '', 'PEPTIDE1{A.[Xyz]}$$$$V2.0',
      'PEPTIDE1{C.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$V2.0', 'RNA1{R(A)P.R(C)P.R(G)}$$$$V2.0',
      'PEPTIDE1{A.C}|PEPTIDE2{G.K}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0', 'PEPTIDE1{[meI].[hHis].[Aca].N.T.[dK]}$$$$V2.0']);
    const df = DG.DataFrame.fromColumns([col]);
    const converter = await seqHelper.getHelmToMolfileConverter(monomerLib);
    const sync = converter.convertToRdKitBeautifiedMolfileColumn(col, true, null, monomerLib);
    // the registered function (the package's own init has run for it: getMolColumnFromHelm reads _package.seqHelper)
    const workers: DG.Column<string> = await grok.functions.call('Bio:getMolFromHelm', {df, helmCol: col});
    expect(workers.name, sync.name);
    expect(workers.semType, DG.SEMTYPE.MOLECULE);
    expect(workers.meta.units, sync.meta.units);
    for (const tag of [ChemTags.SEQUENCE_SRC_COL, SEQUENCE_SRC_LAYOUT_TAG])
      expect(workers.getTag(tag), sync.getTag(tag), tag);
    for (let rowIdx = 0; rowIdx < col.length; ++rowIdx)
      expect(workers.get(rowIdx) ?? '', sync.get(rowIdx) ?? '', `row ${rowIdx}`);
    expect(!!workers.get(3) && !!workers.get(4), true, 'the cyclic and RNA rows build');
  });

  // bio-edge 2026-11-02: a connection to a repeated monomer went to its first copy, R2 too: `A.(G)'3'` closed head to
  // tail failed with "R-group R2 of PEPTIDE:G used twice" (the first copy's R2 joins the second copy). R2 goes to the
  // last copy now: the same molecule as the sequence written out
  test('repeat-r2-connection-last-copy', async () => {
    for (const [repeated, writtenOut] of [
      [`PEPTIDE1{A.(G)'3'}$PEPTIDE1,PEPTIDE1,2:R2-1:R1$$$V2.0`, 'PEPTIDE1{A.G.G.G}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0'],
      [`PEPTIDE1{C.(A.G)'2'}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0`, 'PEPTIDE1{C.A.G.A.G}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$V2.0'],
      [`PEPTIDE1{C.(A.C)'2'.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0`, 'PEPTIDE1{C.A.C.A.C.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0'],
    ]) {
      for (const nonlinear of [false, true]) {
        const {mols, warnings} = await convert(helmCol([repeated, writtenOut]), nonlinear);
        expect(warnings.length, 0, `${repeated}: ${warnings.join('; ')}`);
        expect(await mg.canonicalSmiles(mols[0]), await mg.canonicalSmiles(mols[1]), `${repeated} nonlinear=${nonlinear}`);
      }
    }
  });

  // Fixed repeats expand; the copies keep the source token's position
  test('fixed-repeat-expanded', async () => {
    const col = helmCol(['PEPTIDE1{A\'3\'.G}$$$$V2.0']);
    const {mols, warnings} = await convert(col, false);
    expect(warnings.length, 0, warnings.join('; '));
    expect(await mg.canonicalSmiles(mols[0]),
      await mg.canonicalSmiles('C[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)NCC(=O)O'));
    const map = await mapOf(col, 'linear');
    expectArray(map.keys, [0, 1]);
  });
});
