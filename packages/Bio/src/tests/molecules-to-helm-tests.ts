/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getMonomerLibHelper, IMonomerLib, IMonomerLibBase, IMonomerLibHelper, MonomerLibData}
  from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {getMacromoleculeAssembler} from '@datagrok-libraries/bio/src/monomer-works/assembler';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {IMacromoleculeAssembler} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';

import {helmListColumn} from '../utils/helm-to-molfile/converter/converter';
import {REGENERATED_COL_NAME, convertMoleculesToHelm} from '../utils/molecules-to-helm';
import {_package} from '../package-test';

/** Molecules to HELM (utils/molecules-to-helm.ts: Chem's macromolecule assembler `toHelm`, Crux on Chem's workers).
 * Known peptides of every topology to their exact HELM, with the R numbers of every connection and the same molecule
 * back; Bio's public samples taken to molecules (To Atomic Level's path) and back, as molecules (canonical SMILES with
 * stereo) and as normalized HELM, against the counts measured on localhost (2026-11-18); per-row errors; the three
 * functions called as the platform calls them (the package bundle), their column, end to end on a sample table; the
 * library key shared with To Atomic Level. The library is Bio's test library (HELMCore + polytool-lib). The timing test
 * also times every library the stand's providers list, combined in memory: timing only. */
category('moleculesToHelm', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let seqHelper: ISeqHelper;
  let monomerLib: IMonomerLib;
  let rdKit: RDModule;

  before(async () => {
    seqHelper = await getSeqHelper();
    rdKit = await getRdKitModule();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests(); // HELMCoreLibrary.json + polytool-lib.json
    monomerLib = monomerLibHelper.getMonomerLib();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  /** Canonical SMILES with stereo ('' for none, 'unread' for a molecule that does not read). */
  function canonical(molecule: string | null | undefined): string {
    if (!molecule) return '';
    let mol = null;
    try {
      mol = rdKit.get_mol(molecule);
      return mol ? mol.get_smiles() : 'unread';
    } catch {
      return 'unread';
    } finally {
      mol?.delete();
    }
  }

  /** An inline monomer as a molecule with its attachment points: Datagrok's CXSMILES (`_Rn` labels) and the `[*:n]`
   * SMILES of Bio's samples compare equal when they are one piece with the same R groups. */
  const inlineKey = (token: string): string => `smiles:${canonical(inlineWithMaps(token))}`;

  /** HELM to molecules as To Atomic Level does (non-linear); '' for a row that does not build. */
  async function toMolecules(helms: string[], lib: IMonomerLibBase = monomerLib): Promise<string[]> {
    const col = helmListColumn(helms);
    DG.DataFrame.fromColumns([col]);
    const res = await seqHelper.helmToAtomicLevel(col, false, false, lib);
    expect(res.molCol != null, true, `no molecule column: ${res.warnings.join('; ')}`);
    return Array.from({length: helms.length}, (_, i) => res.molCol!.get(i) ?? '');
  }

  async function readSample(file: string): Promise<DG.DataFrame> {
    const text = await _package.files.readAsText(`samples/${file}`);
    return DG.DataFrame.fromCsv(text, file.endsWith('.tsv') ? {delimiter: '\t'} : {});
  }

  const sampleHelms = async (file: string, colName: string): Promise<string[]> =>
    Array.from((await readSample(file)).col(colName)!.toList(), (h) => h ?? '');

  /** A Molecule column of `values` in a table of its own. */
  function moleculeTable(values: (string | null)[]): { df: DG.DataFrame, molCol: DG.Column<string> } {
    const molCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'molecules', values) as DG.Column<string>;
    molCol.semType = DG.SEMTYPE.MOLECULE;
    return {df: DG.DataFrame.fromColumns([molCol]), molCol};
  }

  // -- known peptides: the exact HELM, each connection's R numbers, the same molecule back --

  /** Inputs written in Datagrok's HELM (public monomers), and the HELM Molecules to HELM writes for their molecules
   * (crux-core f92fd728c), each checked by hand: its monomers, connections and R numbers, an inline piece's `_R`
   * labels on the right atoms; the test builds each back to the input's molecule. */
  const KNOWN: { name: string, helm: string, expected: string, connections: string[] }[] = [
    {name: 'linear', helm: 'PEPTIDE1{A.G.C}$$$$V2.0', expected: 'PEPTIDE1{A.G.C}$$$$V2.0', connections: []},
    {name: 'linear, multi-letter monomers', helm: 'PEPTIDE1{[meI].[hHis].[Aca].N.T}$$$$V2.0',
      expected: 'PEPTIDE1{[meI].[hHis].[Aca].N.T}$$$$V2.0', connections: []},
    {name: 'head-to-tail ring', helm: 'PEPTIDE1{A.G.C.G.F}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$V2.0',
      expected: 'PEPTIDE1{A.G.C.G.F}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$V2.0', connections: ['PEPTIDE1:5:R2-PEPTIDE1:1:R1']},
    {name: 'disulfide', helm: 'PEPTIDE1{C.A.G.C}$PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$V2.0',
      expected: 'PEPTIDE1{C.A.G.C}$PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$V2.0', connections: ['PEPTIDE1:1:R3-PEPTIDE1:4:R3']},
    {name: 'bicyclic: a head-to-tail ring and a disulfide across it',
      helm: 'PEPTIDE1{C.A.G.C.G.F}$PEPTIDE1,PEPTIDE1,6:R2-1:R1|PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$V2.0',
      expected: 'PEPTIDE1{C.A.G.C.G.F}$PEPTIDE1,PEPTIDE1,6:R2-1:R1|PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$V2.0',
      connections: ['PEPTIDE1:6:R2-PEPTIDE1:1:R1', 'PEPTIDE1:1:R3-PEPTIDE1:4:R3']},
    {name: 'two chains joined by a disulfide', helm: 'PEPTIDE1{A.C.G}|PEPTIDE2{G.C.A}$PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$V2.0',
      expected: 'PEPTIDE1{A.C.G}|PEPTIDE2{G.C.A}$PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$V2.0',
      connections: ['PEPTIDE1:2:R3-PEPTIDE2:2:R3']},
    // an S-acetyl on C written as a CHEM entity: HELMCore has the acetyl as its PEPTIDE monomer [ac] (CC(=O)[*:2]),
    // and a library monomer is preferred to an unmatched piece, so it comes back as [ac], its R2 (the carbonyl C) on
    // C's R3 (the S): the same thioester
    {name: 'an S-acetyl on a side chain: HELMCore\'s [ac] on its R2',
      helm: 'PEPTIDE1{A.C.G}|CHEM1{[*C(C)=O |$_R1;;;$|]}$PEPTIDE1,CHEM1,2:R3-1:R1$$$V2.0',
      expected: 'PEPTIDE1{A.C.G}|PEPTIDE2{[ac]}$PEPTIDE1,PEPTIDE2,2:R3-1:R2$$$V2.0',
      connections: ['PEPTIDE1:2:R3-PEPTIDE2:1:R2']},
    // cyclooctylalanine, in no library: inline CXSMILES, atom 0 `*` labelled _R1 on the N, atom 13 `*` labelled _R2 on
    // the carbonyl C
    {name: 'an inline-SMILES monomer on the backbone (cyclooctylalanine, in no library)',
      helm: 'PEPTIDE1{A.G.[[*:1]N[C@@H](CC1CCCCCCC1)C(=O)[*:2]].G}$$$$V2.0',
      expected: 'PEPTIDE1{A.G.[*N[C@@H](CC1CCCCCCC1)C(*)=O |$_R1;;;;;;;;;;;;;_R2;$|].G}$$$$V2.0', connections: []},
    // a 4-(trifluoromethyl)benzoyl on C's S, in no library: a CHEM entity, atom 0 `*` labelled _R1 on the carbonyl C
    {name: 'a CHEM entity on a side chain (4-CF3-benzoyl on C, in no library)',
      helm: 'PEPTIDE1{A.C.G}|CHEM1{[*C(=O)c1ccc(cc1)C(F)(F)F |$_R1;;;;;;;;;;;;$|]}$PEPTIDE1,CHEM1,2:R3-1:R1$$$V2.0',
      expected: 'PEPTIDE1{A.C.G}|CHEM1{[*C(=O)c1ccc(C(F)(F)F)cc1 |$_R1;;;;;;;;;;;;$|]}$PEPTIDE1,CHEM1,2:R3-1:R1$$$V2.0',
      connections: ['PEPTIDE1:2:R3-CHEM1:1:R1']},
    // crosslinks through side chains, their linker in no library (a CHEM entity, its R groups numbered by the partner
    // placed first): an m-xylene bis-thioether staple between C1 and C4 (atom 0 `*` _R1, atom 8 `*` _R2, meta on the
    // ring), a TBMB bicycle, 1,3,5-tris(methylene)benzene on C1, C3 and C5 (atoms 0, 6 and 10: _R1, _R2, _R3), and
    // the staple between the chain's two ends (C1 and C5). Each comes back as written. Before crux-core 3a47d8fce the
    // C1/C4 staple came back as the same molecule written as ONE inline residue (C1, the xylene and C4) in a ring.
    {name: 'a bis-thioether staple between C1 and C4 (m-xylene, in no library)',
      helm: 'PEPTIDE1{C.A.G.C}|CHEM1{[*Cc1cccc(C*)c1 |$_R1;;;;;;;;_R2;$|]}$PEPTIDE1,CHEM1,1:R3-1:R1|PEPTIDE1,CHEM1,4:R3-1:R2$$$V2.0',
      expected: 'PEPTIDE1{C.A.G.C}|CHEM1{[*Cc1cccc(C*)c1 |$_R1;;;;;;;;_R2;$|]}$PEPTIDE1,CHEM1,1:R3-1:R1|PEPTIDE1,CHEM1,4:R3-1:R2$$$V2.0',
      connections: ['PEPTIDE1:1:R3-CHEM1:1:R1', 'PEPTIDE1:4:R3-CHEM1:1:R2']},
    {name: 'a TBMB bicycle: a 1,3,5-tris(methylene)benzene on C1, C3 and C5',
      helm: 'PEPTIDE1{C.A.C.G.C}|CHEM1{[*Cc1cc(C*)cc(C*)c1 |$_R1;;;;;;_R2;;;;_R3;$|]}$PEPTIDE1,CHEM1,1:R3-1:R1|PEPTIDE1,CHEM1,3:R3-1:R2|PEPTIDE1,CHEM1,5:R3-1:R3$$$V2.0',
      expected: 'PEPTIDE1{C.A.C.G.C}|CHEM1{[*Cc1cc(C*)cc(C*)c1 |$_R1;;;;;;_R2;;;;_R3;$|]}$PEPTIDE1,CHEM1,1:R3-1:R1|PEPTIDE1,CHEM1,3:R3-1:R2|PEPTIDE1,CHEM1,5:R3-1:R3$$$V2.0',
      connections: ['PEPTIDE1:1:R3-CHEM1:1:R1', 'PEPTIDE1:3:R3-CHEM1:1:R2', 'PEPTIDE1:5:R3-CHEM1:1:R3']},
    {name: 'a staple between the chain\'s two ends (m-xylene on C1 and C5)',
      helm: 'PEPTIDE1{C.A.G.G.C}|CHEM1{[*Cc1cccc(C*)c1 |$_R1;;;;;;;;_R2;$|]}$PEPTIDE1,CHEM1,1:R3-1:R1|PEPTIDE1,CHEM1,5:R3-1:R2$$$V2.0',
      expected: 'PEPTIDE1{C.A.G.G.C}|CHEM1{[*Cc1cccc(C*)c1 |$_R1;;;;;;;;_R2;$|]}$PEPTIDE1,CHEM1,1:R3-1:R1|PEPTIDE1,CHEM1,5:R3-1:R2$$$V2.0',
      connections: ['PEPTIDE1:1:R3-CHEM1:1:R1', 'PEPTIDE1:5:R3-CHEM1:1:R2']},
    // a CHEM linker attached through a specified centre that bears both R groups: CH(CH3) between C2's S (R1) and K4's
    // N (R2), written [*:1][C@@H](C)[*:2]. Each `*` of such a centre is written [*:n] beside its _Rn so the centre
    // keeps its tag (crux-core 03c1fcf02; before, the tag was lost and the row refused): C[C@H]([*:1])[*:2] is the
    // input's centre (one swap of neighbours, @@ to @), atoms 2 and 3 labelled _R1 and _R2
    {name: 'a CHEM linker attached through a specified centre bearing both R groups',
      helm: 'PEPTIDE1{A.C.G.K.A}|CHEM1{[[*:1][C@@H](C)[*:2]]}$PEPTIDE1,CHEM1,2:R3-1:R1|PEPTIDE1,CHEM1,4:R3-1:R2$$$V2.0',
      expected: 'PEPTIDE1{A.C.G.K.A}|CHEM1{[C[C@H]([*:1])[*:2] |$;;_R1;_R2$|]}$PEPTIDE1,CHEM1,2:R3-1:R1|PEPTIDE1,CHEM1,4:R3-1:R2$$$V2.0',
      connections: ['PEPTIDE1:2:R3-CHEM1:1:R1', 'PEPTIDE1:4:R3-CHEM1:1:R2']},
  ];

  test('known peptides: the exact HELM, each connection\'s R numbers, the same molecule back', async () => {
    const mols = await toMolecules(KNOWN.map((k) => k.helm));
    mols.forEach((m, i) => expect(!!m, true, `${KNOWN[i].name}: To Atomic Level built no molecule`));
    const res = await convertMoleculesToHelm(mols, monomerLib);
    const back = await toMolecules(res.helms.map((h) => h ?? ''));
    KNOWN.forEach((k, i) => {
      expect(res.errors[i], null, `${k.name}: error`);
      expect(res.helms[i], k.expected, k.name);
      expect(connectionsOf(res.helms[i]).join(' | '), k.connections.join(' | '), `${k.name}: connections and R numbers`);
      expect(canonical(back[i]), canonical(mols[i]), `${k.name}: the regenerated HELM's molecule`);
    });
    // the counts: library monomers, inline-SMILES monomers, CHEM entities
    const counts = (i: number): string => [res.matched[i], res.inline[i], res.chem[i]].join();
    for (let i = 0; i < 6; i++) {
      const monomers = monomersOf(res.helms[i], 0).length + monomersOf(res.helms[i], 1).length;
      expect(counts(i), `${monomers},0,0`, `${KNOWN[i].name}: every monomer from the library`);
    }
    expect(counts(6), '4,0,0', 'the S-acetyl\'s row: A, C, G and [ac] from the library');
    expect(counts(7), '3,1,0', 'the inline monomer\'s row: matched, inline, CHEM');
    expect(counts(8), '3,0,1', 'the CHEM entity\'s row: matched, inline, CHEM');
    expect(counts(9), '4,0,1', 'the C1/C4 staple: C, A, G, C matched, the linker a CHEM entity');
    expect(counts(10), '5,0,1', 'the TBMB bicycle: five residues matched, the linker a CHEM entity');
    expect(counts(11), '5,0,1', 'the end-to-end staple: five residues matched, the linker a CHEM entity');
    expect(counts(12), '5,0,1', 'the chiral linker: five residues matched, the linker a CHEM entity');
    const linker = monomersOf(res.helms[12], 1);
    expect(linker.length === 1 && canonical(inlineWithMaps(linker[0])) === canonical('[*:1][C@@H](C)[*:2]'), true,
      `the chiral linker keeps its centre, R1 and R2 where they were: ${linker}`);
    // the pieces no library monomer is, as molecules with their attachment points (R groups as [*:n])
    const entity = monomersOf(res.helms[8], 1);
    expect(entity.length === 1 &&
      canonical(inlineWithMaps(entity[0])) === canonical('[*:1]C(=O)c1ccc(C(F)(F)F)cc1'), true,
    `the CHEM entity is the 4-CF3-benzoyl, its R1 on the carbonyl C: ${entity}`);
    const chain = monomersOf(res.helms[7], 0);
    expect(chain.length === 4 && chain[0] === 'A' && chain[1] === 'G' && chain[3] === 'G', true, `the chain: ${chain}`);
    expect(canonical(inlineWithMaps(chain[2])), canonical('[*:1]N[C@@H](CC1CCCCCCC1)C(=O)[*:2]'),
      `the inline monomer, its R1 on the N and its R2 on the carbonyl C: ${chain[2]}`);
  }, {timeout: 120000});

  // -- Bio's public samples: HELM -> molecule -> HELM -> molecule --

  /** The round trip of `helms`: how many rows build, convert back, give the same molecule and the same HELM. */
  async function roundTrip(helms: string[]):
    Promise<{ built: number, regenerated: number, molEqual: number, helmEqual: number, ms: number,
      errorRows: number[], errors: string[] }> {
    const mols = await toMolecules(helms);
    const rows = mols.map((m, i) => m ? i : -1).filter((i) => i >= 0);
    const t0 = performance.now();
    const res = await convertMoleculesToHelm(rows.map((i) => mols[i]), monomerLib);
    const ms = performance.now() - t0;
    const regen = res.helms.map((h) => h ?? '');
    const back = await toMolecules(regen);
    let molEqual = 0;
    let helmEqual = 0;
    rows.forEach((row, k) => {
      const a = canonical(mols[row]);
      if (a && a !== 'unread' && a === canonical(back[k])) molEqual++;
      if (regen[k] && normalizeHelm(regen[k], inlineKey) === normalizeHelm(helms[row], inlineKey)) helmEqual++;
    });
    const failed = rows.map((row, k) => res.errors[k] ? row : -1).filter((row) => row >= 0);
    return {built: rows.length, regenerated: regen.filter((h) => h).length, molEqual, helmEqual, ms,
      errorRows: failed, errors: res.errors.filter((e): e is string => !!e)};
  }

  /** Bio's public samples and what was measured on localhost (2026-11-17, Chem with crux-core f92fd728c, this test
   * library): rows To Atomic Level builds, rows that come back as the same molecule, rows whose HELM is the same after
   * normalizing (inline monomers compared as molecules). The measured values are the floors: a count under them is a
   * regression. Where the same molecule comes back as another HELM, it is under another symbol of one structure
   * (`Bmt` for `Bmt_E`, `am` for `NH2`, `gGlu` for an E on R3, `meA` and `D-aThr` for inline pieces the library has)
   * or a CHEM entity for an inline piece with one backbone R group. The 3 RNA rows of HELM-non-linear.csv are out of
   * the converter's scope (PEPTIDE and CHEM): each is its row's error ("joins none"), never a HELM of another
   * molecule, and every other row converts (`errorRows`: rows of the file). Bio's old Python script on the same rows,
   * for scale: molecule equal 539, 2, 0, 41, 2, 1, the RNA rows HELM of other molecules. */
  type Sample = { file: string, col: string, built: number, molEqual: number, helmEqual: number, errorRows: number[] };
  const SAMPLES: Sample[] = [
    {file: 'HELM.csv', col: 'HELM', built: 540, molEqual: 540, helmEqual: 535, errorRows: []},
    {file: 'HELM-non-linear.csv', col: 'HELM', built: 7, molEqual: 4, helmEqual: 4, errorRows: [4, 5, 6]},
    {file: 'HELM_BI_CYCLIC.csv', col: 'HELM', built: 6, molEqual: 6, helmEqual: 5, errorRows: []},
    {file: 'HELM_cyclic.tsv', col: 'sequence', built: 41, molEqual: 41, helmEqual: 37, errorRows: []},
    {file: 'HELM_CHEMS.csv', col: 'HELM', built: 9, molEqual: 9, helmEqual: 3, errorRows: []},
    {file: 'BILN_W_HELM.csv', col: 'helm(BILN)', built: 9, molEqual: 9, helmEqual: 7, errorRows: []},
  ];

  for (const sample of SAMPLES) {
    test(`round trip: ${sample.file}`, async () => {
      const r = await roundTrip(await sampleHelms(sample.file, sample.col));
      const line = `${sample.file}: built ${r.built}, regenerated ${r.regenerated}, molecule equal ${r.molEqual}, ` +
        `HELM equal ${r.helmEqual}, ${r.ms.toFixed(0)} ms`;
      console.log(line);
      expect(r.built, sample.built, `${line}: rows built`);
      expect(r.molEqual >= sample.molEqual, true, `${line}: the same molecule in at least ${sample.molEqual}`);
      expect(r.helmEqual >= sample.helmEqual, true, `${line}: the same HELM in at least ${sample.helmEqual}`);
      expect(r.errorRows.join(), sample.errorRows.join(), `${line}: the rows that are errors`);
      expect(r.errors.every((e) => e.includes('joins none')), true,
        `${line}: the errors: ${r.errors.join('; ').slice(0, 300)}`);
      return line;
    }, {timeout: 120000});
  }

  // -- per-row errors, the three functions as the platform calls them, the column --

  test('per-row errors: empty and bad molecules, nothing thrown', async () => {
    const values = ['', 'C[C@H](N)C(=O)NCC(=O)O', 'not a molecule', 'C[C@H](N)C(=O)NCC(=O)O.Cl',
      'NCC(=O)N[C@@H](CS)C(=O)O', null];
    const res = await convertMoleculesToHelm(values, monomerLib);
    expect(res.helms.length, values.length, 'a value per row');
    expect(res.helms[0] === '' && res.errors[0] === null, true, 'an empty row: no HELM, no error');
    expect(res.helms[1], 'PEPTIDE1{A.G}$$$$V2.0', 'A-G');
    expect(res.helms[2] === null && String(res.errors[2]).includes('cannot be read'), true,
      `no molecule: an error (${res.errors[2]})`);
    expect(res.helms[3] === null && String(res.errors[3]).includes('joins none'), true,
      `a counter-ion that joins nothing: an error, not a silent drop (${res.errors[3]})`);
    expect(res.helms[4], 'PEPTIDE1{G.C}$$$$V2.0', 'G-C');
    expect(res.helms[5] === '' && res.errors[5] === null, true, 'a null row: no HELM, no error');
    // the column function, as the platform calls it: the failed rows' cells are empty, the others hold their HELM
    const {df, molCol} = moleculeTable(values);
    const col: DG.Column<string> = await grok.functions.call('Bio:moleculesToHelmColumn', {table: df, molecules: molCol});
    expect(Array.from({length: col.length}, (_, i) => col.get(i) ?? '').join(' | '),
      ['', 'PEPTIDE1{A.G}$$$$V2.0', '', '', 'PEPTIDE1{G.C}$$$$V2.0', ''].join(' | '), 'the column\'s cells');
  });

  test('out of scope: an RNA strand is the row\'s error, never a HELM of another molecule', async () => {
    const [mol] = await toMolecules(['RNA1{d(A)p.d(G)p.d(C)p}$$$$V2.0']);
    expect(!!mol, true, 'To Atomic Level built the strand');
    const res = await convertMoleculesToHelm([mol, 'C[C@H](N)C(=O)NCC(=O)O'], monomerLib);
    expect(res.helms[0] === null && String(res.errors[0]).startsWith('a fragment matches no monomer and joins none'),
      true, `the strand: ${res.helms[0]} ${String(res.errors[0]).slice(0, 120)}`);
    expect([res.matched[0], res.inline[0], res.chem[0]].join(), '0,0,0', 'the strand\'s counts');
    expect(res.helms[1], 'PEPTIDE1{A.G}$$$$V2.0', 'the next row converts');
  });

  test('Molecule to HELM Single: linear, head-to-tail, empty and bad molecules', async () => {
    const single = (molecule: string): Promise<string> =>
      grok.functions.call('Bio:moleculeToHelmSingle', {molecule});
    expect(await single('C[C@H](N)C(=O)NCC(=O)O'), 'PEPTIDE1{A.G}$$$$V2.0', 'A-G from SMILES');
    const cyclic = 'PEPTIDE1{A.G.C.G.F}$PEPTIDE1,PEPTIDE1,5:R2-1:R1$$$V2.0';
    const [mol] = await toMolecules([cyclic]);
    expect(await single(mol), cyclic, 'head-to-tail from its MOL block');
    expect(await single(''), '', 'an empty molecule');
    expect(await single('not a molecule'), '', 'no molecule');
  });

  test('the column: name, semantic type, units, renderer; a second run', async () => {
    const helms = ['PEPTIDE1{A.G.C}$$$$V2.0', 'PEPTIDE1{F.F}$$$$V2.0'];
    const {df, molCol} = moleculeTable(await toMolecules(helms));
    molCol.meta.units = DG.UNITS.Molecule.MOLBLOCK;
    await grok.functions.call('Bio:moleculesToHelmTopMenu', {table: df, molecules: molCol});
    const first = df.col(REGENERATED_COL_NAME)!;
    expect(first != null, true, `no ${REGENERATED_COL_NAME} column: ${df.columns.names()}`);
    expect(first.semType, DG.SEMTYPE.MACROMOLECULE, 'semType');
    expect(first.meta.units, NOTATION.HELM, 'units');
    expect(first.getTag('cell.renderer'), 'helm', 'renderer');
    expect(first.get(0), helms[0], 'row 0');
    expect(first.get(1), helms[1], 'row 1');
    const second: DG.Column<string> = await grok.functions.call('Bio:moleculesToHelmColumn', {table: df, molecules: molCol});
    expect(second.name, `${REGENERATED_COL_NAME} (2)`, 'the second column\'s name');
    expect(df.columns.length, 3, 'one column a run');
    expect(second.get(1), helms[1], 'the same HELM');
  });

  test('end to end on a sample table: To Atomic Level, then Molecules to HELM, the same molecules', async () => {
    const df = await readSample('HELM_cyclic.tsv');
    const seqCol = df.col('sequence')!;
    await grok.data.detectSemanticTypes(df);
    expect(seqCol.semType, DG.SEMTYPE.MACROMOLECULE, 'the sequences detected');
    await grok.functions.call('Bio:toAtomicLevel', {table: df, seqCol, nonlinear: true, highlight: false});
    const molCol = df.columns.toList().find((c) => c.semType === DG.SEMTYPE.MOLECULE) as DG.Column<string>;
    expect(molCol != null, true, `no molecule column: ${df.columns.names()}`);
    await grok.functions.call('Bio:moleculesToHelmTopMenu', {table: df, molecules: molCol});
    const helmCol = df.col(REGENERATED_COL_NAME) as DG.Column<string>;
    expect(helmCol != null && helmCol.semType === DG.SEMTYPE.MACROMOLECULE, true, `no HELM column: ${df.columns.names()}`);
    const back = await toMolecules(Array.from({length: df.rowCount}, (_, i) => helmCol.get(i) ?? ''));
    let same = 0;
    for (let i = 0; i < df.rowCount; i++)
      if (canonical(molCol.get(i)) && canonical(molCol.get(i)) === canonical(back[i])) same++;
    expect(same, df.rowCount, 'the same molecule in every row');
  }, {timeout: 120000});

  test('the library: the key and JSON To Atomic Level registers (one compile per worker for both directions)', async () => {
    // both directions as the platform calls them (the package bundle, its library cache), the assembler's
    // setLibrary watched
    const assembler: IMacromoleculeAssembler = await getMacromoleculeAssembler();
    const calls: { key: string, json: string }[] = [];
    const setLibrary = assembler.setLibrary;
    assembler.setLibrary = async (key: string, json: string) => {
      calls.push({key, json});
      return setLibrary.call(assembler, key, json);
    };
    try {
      const seqCol = helmListColumn(['PEPTIDE1{A.G.C}$$$$V2.0'], 'sequence');
      const df = DG.DataFrame.fromColumns([seqCol]);
      await grok.functions.call('Bio:toAtomicLevel', {table: df, seqCol, nonlinear: true, highlight: false});
      const toAtomic = calls.length;
      const molCol = df.columns.toList().find((c) => c.semType === DG.SEMTYPE.MOLECULE)!;
      const helmCol: DG.Column<string> = await grok.functions.call('Bio:moleculesToHelmColumn', {table: df, molecules: molCol});
      expect(helmCol.get(0), 'PEPTIDE1{A.G.C}$$$$V2.0', 'converted');
      expect(toAtomic > 0 && calls.length > toAtomic, true, `setLibrary calls: To Atomic Level ${toAtomic}, in all ${calls.length}`);
      const toHelmCall = calls[calls.length - 1];
      expect(toHelmCall.key, calls[toAtomic - 1].key, 'Molecules to HELM\'s key is To Atomic Level\'s');
      expect(toHelmCall.json === calls[toAtomic - 1].json, true, 'the same JSON, not another copy of the library');
    } finally {
      assembler.setLibrary = setLibrary;
    }
  });

  // -- timing --

  /** Every library the providers list, combined in memory (HELMCore first, a symbol's first monomer kept): the
   * libraries a user who enabled them all would have. Timing only. */
  async function everyLibrary(): Promise<{ lib: IMonomerLibBase, libraries: number, monomers: number }> {
    const names: [string, string][] = [];
    for (const provider of await monomerLibHelper.getProviders()) {
      for (const name of await provider.listLibraries())
        names.push([provider.name, name]);
    }
    names.sort((a, b) => Number(b[1] === 'HELMCoreLibrary.json') - Number(a[1] === 'HELMCoreLibrary.json'));
    const data: MonomerLibData = {};
    let libraries = 0;
    let monomers = 0;
    for (const [provider, name] of names) {
      const one = await monomerLibHelper.readSingleLibrary(provider, name).catch(() => null);
      if (!one) continue;
      libraries++;
      for (const m of one.toJSON()) {
        const byType = data[m.polymerType] ??= {};
        if (byType[m.symbol]) continue;
        byType[m.symbol] = JSON.parse(JSON.stringify(m)); // override() sets its own fields on what it gets
        monomers++;
      }
    }
    return {lib: monomerLib.override(data, 'moleculesToHelm timing: every library'), libraries, monomers};
  }

  test('timing: 1,000 molecules, the test library and every library of the stand', async () => {
    const mols = (await toMolecules(await sampleHelms('HELM.csv', 'HELM'))).filter((m) => m);
    const input = Array.from({length: 1000}, (_, i) => mols[i % mols.length]);
    const time = async (lib: IMonomerLibBase): Promise<{ ms: number, ok: number }> => {
      const t0 = performance.now();
      const res = await convertMoleculesToHelm(input, lib);
      return {ms: performance.now() - t0, ok: res.helms.filter((h) => h).length};
    };
    const rate = (r: { ms: number }): string => (1000 * input.length / r.ms).toFixed(0);
    const first = await time(monomerLib); // the library's first column: each worker compiles it and its index
    const again = await time(monomerLib);
    const every = await everyLibrary();
    const everyFirst = await time(every.lib);
    const everyAgain = await time(every.lib);
    const line = `1,000 molecules: the test library ${first.ms.toFixed(0)} ms first, ${again.ms.toFixed(0)} ms again ` +
      `(${rate(again)} molecules/s, ${again.ok} converted); every library of the stand (${every.libraries} libraries, ` +
      `${every.monomers} monomers; timing only) ${everyFirst.ms.toFixed(0)} ms first, ${everyAgain.ms.toFixed(0)} ms ` +
      `again (${rate(everyAgain)} molecules/s; ${(everyAgain.ms / again.ms).toFixed(1)}x the test library's time)`;
    console.log(line);
    expect(again.ok, input.length, `${line}: every molecule converted with the test library`);
    expect(everyAgain.ok, input.length, `${line}: every molecule converted with every library`);
    // Measured on localhost (2026-11-17), these 1,000 rows warm: Crux about 90 ms with the test library and 410 ms
    // with every library of the stand (9,783 monomers); Bio's old Python script took 7,135 ms and 14,785 ms. The bounds
    // are the owner's (crux-core decision mol-to-helm, 8): at least 1,000 molecules a second with the test library,
    // and every library of the stand within 3 s.
    expect(again.ms < 1000, true, `${line}: under 1,000 molecules/s with the test library`);
    expect(everyAgain.ms < 3000, true, `${line}: over 3 s with every library`);
    return line;
  }, {timeout: 300000});
});

// -- HELM helpers --

/** HELM's top-level parts split at `sep`, not inside brackets, braces, parentheses or quotes (an inline monomer's
 * CXSMILES holds `$` and `|`). */
function splitHelm(text: string, sep: string): string[] {
  const parts: string[] = [];
  let depth = 0;
  let quoted = false;
  let cur = '';
  for (const ch of text) {
    if (ch === '"') quoted = !quoted;
    else if (!quoted && '[({'.includes(ch)) depth++;
    else if (!quoted && '])}'.includes(ch)) depth--;
    else if (!quoted && depth === 0 && ch === sep) {
      parts.push(cur);
      cur = '';
      continue;
    }
    cur += ch;
  }
  parts.push(cur);
  return parts;
}

/** The connections of `helm` as written, `PEPTIDE1:4:R2-PEPTIDE1:1:R1` each, in order. */
function connectionsOf(helm: string | null): string[] {
  const section = splitHelm(helm ?? '', '$')[1] ?? '';
  return section ? splitHelm(section, '|').map((c) => {
    const m = /^(\w+),(\w+),(\d+):(R\d+)-(\d+):(R\d+)$/.exec(c);
    return m ? `${m[1]}:${m[3]}:${m[4]}-${m[2]}:${m[5]}:${m[6]}` : `unread: ${c}`;
  }) : [];
}

/** The monomers of chain `chain` (0-based) of `helm`: symbols without their brackets, inline SMILES as written. */
function monomersOf(helm: string | null, chain: number): string[] {
  const polymer = splitHelm(splitHelm(helm ?? '', '$')[0], '|')[chain] ?? '';
  const body = /^\w+\{(.*)\}$/s.exec(polymer)?.[1];
  if (body === undefined) return []; // no such chain
  return splitHelm(body, '.').map((t) => t.startsWith('[') && t.endsWith(']') ? t.slice(1, -1) : t);
}

/** An inline monomer as Datagrok writes it (CXSMILES, its R groups `*` atoms labelled `_R1`, `_R2`, ... in `|$...$|`)
 * as plain SMILES with those atoms written `[*:1]`, `[*:2]`, ...: a molecule any engine reads, its attachment points
 * kept, so a monomer compares by canonical SMILES. Inline SMILES already in the `[*:n]` form are returned as they
 * are.
 * @param {string} inline
 * @return {string} */
function inlineWithMaps(inline: string): string {
  const m = /^(.*?)\s*\|\$([^$]*)\$\|\s*$/s.exec(inline);
  if (!m) return inline;
  const labels = m[2].split(';');
  const relabel = (token: string, label: string | undefined): string => {
    const r = /^_R(\d+)$/.exec(label ?? '');
    return r && (token === '*' || token === '[*]') ? `[*:${r[1]}]` : token;
  };
  const smiles = m[1];
  let atom = -1;
  let out = '';
  for (let i = 0; i < smiles.length;) {
    const ch = smiles[i];
    if (ch === '[') { // a bracket atom
      const j = smiles.indexOf(']', i);
      out += relabel(smiles.slice(i, j + 1), labels[++atom]);
      i = j + 1;
    } else if (smiles.startsWith('Cl', i) || smiles.startsWith('Br', i)) {
      ++atom;
      out += smiles.slice(i, i + 2);
      i += 2;
    } else if (/[BCNOPSFIbcnops*]/.test(ch)) { // an organic-subset atom
      out += relabel(ch, labels[++atom]);
      i++;
    } else { // bonds, branches, ring closures
      out += ch;
      i++;
    }
  }
  return out;
}

/** HELM compared modulo the writer's choices: chains renumbered in a canonical order, a head-to-tail ring rotated to
 * its smallest rotation, connections as unordered pairs, brackets, annotations and the version left out. An inline
 * SMILES monomer (a token with `*` atoms) is replaced by `inline(token)` (as text without it).
 * @param {string} helm
 * @param {function} inline
 * @return {string} */
function normalizeHelm(helm: string, inline?: (token: string) => string): string {
  if (!helm) return '';
  type End = { chain: string, idx: number, r: string };
  const sections = splitHelm(helm.trim(), '$');
  const chains: { id: string, type: string, tokens: string[] }[] = [];
  for (const p of splitHelm(sections[0] ?? '', '|')) {
    const m = /^([A-Z]+)(\d+)\{(.*)\}(".*")?$/s.exec(p.trim());
    if (!m) return `unparsed: ${helm}`;
    const tokens = splitHelm(m[3], '.').map((t) => {
      t = t.trim();
      const bare = t.startsWith('[') && t.endsWith(']') && splitHelm(t.slice(1, -1), '.').length === 1 ? t.slice(1, -1) : t;
      return inline && bare.includes('*') ? inline(bare) : bare;
    });
    chains.push({id: `${m[1]}${m[2]}`, type: m[1], tokens});
  }
  const conns: [End, End][] = [];
  for (const c of splitHelm(sections[1] ?? '', '|')) {
    if (!c.trim()) continue;
    const m = /^([A-Z]+\d+),([A-Z]+\d+),(\d+):([A-Za-z0-9?]+)-(\d+):([A-Za-z0-9?]+)$/.exec(c.trim());
    if (!m) return `unparsed: ${helm}`;
    conns.push([{chain: m[1], idx: +m[3], r: m[4]}, {chain: m[2], idx: +m[5], r: m[6]}]);
  }
  for (const ch of chains) {
    const n = ch.tokens.length;
    const closure = conns.find(([a, b]) => a.chain === ch.id && b.chain === ch.id &&
      ((a.idx === n && a.r === 'R2' && b.idx === 1 && b.r === 'R1') || (b.idx === n && b.r === 'R2' && a.idx === 1 && a.r === 'R1')));
    if (!closure || n < 2) continue;
    let best = 0;
    let bestKey: string | null = null;
    for (let s = 0; s < n; s++) {
      const key = [...ch.tokens.slice(s), ...ch.tokens.slice(0, s)].join('\u0001');
      if (bestKey === null || key < bestKey) {
        bestKey = key;
        best = s;
      }
    }
    if (best === 0) continue;
    ch.tokens = [...ch.tokens.slice(best), ...ch.tokens.slice(0, best)];
    for (const pair of conns) {
      for (const e of pair) {
        if (e.chain === ch.id && e !== closure[0] && e !== closure[1])
          e.idx = ((e.idx - 1 - best + n) % n) + 1;
      }
    }
  }
  const order = [...chains].sort((a, b) => (a.type + a.tokens.join('.')).localeCompare(b.type + b.tokens.join('.')));
  const rename = new Map<string, string>();
  const counts: { [type: string]: number } = {};
  for (const ch of order) rename.set(ch.id, `${ch.type}${(counts[ch.type] = (counts[ch.type] ?? 0) + 1)}`);
  const ends = conns.map((pair) => pair.map((e) => `${rename.get(e.chain) ?? e.chain}:${e.idx}:${e.r}`).sort().join('-')).sort();
  return order.map((ch) => `${rename.get(ch.id)}{${ch.tokens.join('.')}}`).join('|') + '$' + ends.join('|');
}
