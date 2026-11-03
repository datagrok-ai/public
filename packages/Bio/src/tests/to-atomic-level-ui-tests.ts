import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/monomer-library';
import {sequenceToMolfile} from '../utils/sequence-to-mol';
import {EMBED_3D_FAILED, EMBED_3D_TOO_LONG, embedMolecule3d, showMoleculeDialog, toAtomicLevelSingle}
  from '../widgets/to-atomic-level-widget';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {heavyAtomBondCounts} from './utils/molblock-graph';

type TestDataTargetType = { atomCount: number, bondCount: number };
type TestDataType = {
  src: { seq: string, units: NOTATION },
  tgt: TestDataTargetType,
};

category('toAtomicLevel-ui', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let seqHelper: ISeqHelper;
  let monomerLib: IMonomerLib;

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries

    monomerLib = monomerLibHelper.getMonomerLib();
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  const tests: { [testName: string]: TestDataType } = {
    'fasta': {
      src: {seq: 'MDYKETLLMPK', units: NOTATION.FASTA,},
      tgt: {atomCount: 94, bondCount: 95,},
    },
    'fasta-with-gap': {
      src: {seq: 'MD-YKETLLMPK', units: NOTATION.FASTA,},
      tgt: {atomCount: 94, bondCount: 95,},
    },
    'helm': {
      src: {seq: 'PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2}$$$$', units: NOTATION.HELM,},
      tgt: {atomCount: 68, bondCount: 68,},
    },
    'helm-with-gap': {
      src: {seq: 'PEPTIDE1{meI.hHis.*.Aca.N.T.dK.Thr_PO3H2}$$$$', units: NOTATION.HELM,},
      tgt: {atomCount: 68, bondCount: 68,},
    },
  };

  //   const fastaCsv = `seq
  // MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW
  // MIEVFLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL
  // MMELVLKTIIGPIVVGVVLRIVDKWLNKDK
  // `;
  //   const helmCsv = `seq
  // PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.Aze.dV.E.N.dV.Phe_4Me}$$$$
  // PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.meM.D-Chg.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$$
  // PEPTIDE1{meI.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.D-Cit.N.D-Orn.D-aThr.Phe_4Me}$$$$
  // `;

  const getSeqCol = async (testData: TestDataType): Promise<DG.Column<string>> => {
    const seq = testData.src.seq;
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [seq])]);
    await grok.data.detectSemanticTypes(df);
    return df.getCol('seq');
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}-linear`, async () => {
      const seqCol = await getSeqCol(testData);
      await _testToAtomicLevelFunc(seqCol, false, testData.tgt);
    });
  }

  // The non-linear layout no longer needs the Helm package (the assembler builds both layouts)
  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}-nonlinear`, async () => {
      const seqCol = await getSeqCol(testData);
      await _testToAtomicLevelFunc(seqCol, true, testData.tgt);
    });
  }

  /** The conversion's warning balloons that were not on screen in {@link shown}. */
  const newConversionBalloons = (shown: Set<Element>): Element[] =>
    Array.from(document.querySelectorAll('.d4-balloon'))
      .filter((b) => !shown.has(b) && (b.textContent ?? '').includes('Molfile conversion resulted in'));
  const balloonsNow = (): Set<Element> => new Set(document.querySelectorAll('.d4-balloon'));
  const seqColOf = async (seqs: string[]): Promise<DG.Column<string>> => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', seqs)]);
    await grok.data.detectSemanticTypes(df);
    return df.getCol('seq');
  };
  const panelValue = async (seq: string): Promise<DG.SemanticValue> =>
    DG.SemanticValue.fromTableCell((await seqColOf([seq])).dataFrame.cell(0, 'seq'));
  const UNKNOWN_MONOMER_HELM = 'PEPTIDE1{A.Unknown_Monomer_X.C}$$$$';

  // The panels (Molecular Structure, Molecular 3D Structure) show the error themselves: no balloon per row
  test('panel-failed-sequence-no-balloon', async () => {
    const shown = balloonsNow();
    const res = await toAtomicLevelSingle(await panelValue(UNKNOWN_MONOMER_HELM));
    expect(!!res.errorText && !res.mol, true, `the unknown monomer was built: ${res.errorText}`);
    await delay(500);
    expect(newConversionBalloons(shown).length, 0, 'the panel showed the warning balloon');
  });

  // bio-edge 2026-11-02: the panels said 'No structure generated' for every failed row, without its reason
  test('panel-failed-sequence-reason', async () => {
    for (const [seq, reason] of [
      [UNKNOWN_MONOMER_HELM, 'unknown monomer PEPTIDE:Unknown_Monomer_X'],
      ['PEPTIDE1{}|PEPTIDE2{A.C}$$$$', 'polymer PEPTIDE1 has no monomers'],
      ['PEPTIDE1{A.(G,C).A}$$$$', 'ambiguous monomer'],
      ['PEPTIDE1{A.C.D', 'HELM unterminated-brace'],
    ]) {
      const col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [seq]) as DG.Column<string>;
      col.semType = DG.SEMTYPE.MACROMOLECULE;
      col.meta.units = NOTATION.HELM;
      DG.DataFrame.fromColumns([col]);
      const res = await toAtomicLevelSingle(DG.SemanticValue.fromTableCell(col.dataFrame.cell(0, 'seq')));
      expect(res.mol, '', seq);
      expect(res.errorText.startsWith(`No structure generated: ${reason}`), true, `${seq}: ${res.errorText}`);
    }
  });

  test('column-failed-rows-balloon', async () => {
    const shown = balloonsNow();
    const seqCol = await seqColOf([UNKNOWN_MONOMER_HELM, 'PEPTIDE1{A.C}$$$$']);
    const res = await sequenceToMolfile(seqCol.dataFrame, seqCol, false, false, monomerLib, seqHelper);
    expect(res.warnings.length, 1);
    await awaitCheck(() => newConversionBalloons(shown).length > 0, 'no warning balloon for a column', 3000);
  });

  // 120 monomers; an all-RNA HELM row counts sugar, base and phosphate, so 360 positions
  test('panel-monomer-limit', async () => {
    const protein = (n: number) => 'MDYKETLLMPK'.repeat(12).slice(0, n);
    const fasta121 = await toAtomicLevelSingle(await panelValue(protein(121)));
    expect(fasta121.errorText, 'Maximum number of monomers for molecular conversion is 120');
    expect(fasta121.scrollableLength, 121);
    const fasta120 = await toAtomicLevelSingle(await panelValue(protein(120)));
    expect(!!fasta120.mol, true, `120 monomers: ${fasta120.errorText}`);
    const peptide121 = await toAtomicLevelSingle(await panelValue(`PEPTIDE1{${Array(121).fill('A').join('.')}}$$$$`));
    expect(peptide121.errorText, 'Maximum number of monomers for molecular conversion is 120');
    const rna = (n: number) => `RNA1{${Array(n).fill('r(A)p').join('.')}}$$$$`;
    const rna121 = await toAtomicLevelSingle(await panelValue(rna(121)));
    expect(rna121.errorText, 'Maximum number of monomers for molecular conversion is 360');
    const rna120 = await toAtomicLevelSingle(await panelValue(rna(120)));
    expect(!!rna120.mol, true, `120 nucleotides: ${rna120.errorText}`);
  });

  // bio-edge 2026-11-02: a fixed repeat is one position but many monomers: `A'3000'` passed the limit and the panel
  // built and drew a 15,000-atom molecule
  test('panel-monomer-limit-repeats', async () => {
    const limitMsg = 'Maximum number of monomers for molecular conversion is 120';
    for (const [helm, overLimit] of [
      [`PEPTIDE1{A'121'}$$$$`, true], [`PEPTIDE1{A'3000'}$$$$`, true], [`PEPTIDE1{(A.G)'61'}$$$$`, true],
      [`PEPTIDE1{A'120'}$$$$`, false], [`PEPTIDE1{G.(A.G)'59'.G}$$$$`, false],
    ] as [string, boolean][]) {
      const res = await toAtomicLevelSingle(await panelValue(helm));
      if (overLimit)
        expect(res.errorText, limitMsg, helm);
      else
        expect(!!res.mol && !res.errorText, true, `${helm}: ${res.errorText}`);
    }
    const rnaRes = await toAtomicLevelSingle(await panelValue(`RNA1{(r(A)p)'121'}$$$$`));
    expect(rnaRes.errorText, 'Maximum number of monomers for molecular conversion is 360');
  });

  // bio-edge 2026-11-02: an embedding that failed or timed out left the 3D panel an empty viewer, one that never
  // finished (RDKit mode's server script has no deadline) a spinner
  test('panel-3d-embedding-failed-message', async () => {
    const molBlock2D = 'CC(N)C(=O)O';
    for (const [what, embed] of [
      ['throws', async () => {
        throw new Error('chemEmbed3d: timed out');
      }],
      ['empty', async () => ''],
      ['blank', async () => '\n  \n'],
    ] as [string, (m: string) => Promise<string>][]) {
      const res = await embedMolecule3d(molBlock2D, embed);
      expect(res.error, EMBED_3D_FAILED, what);
      expect(res.molBlock3D, null, what);
    }
    const never = await embedMolecule3d(molBlock2D, () => new Promise<string>(() => {}), 200);
    expect(never.error, EMBED_3D_TOO_LONG, 'an embedding that never finishes');
    const ok = await embedMolecule3d(molBlock2D, async (m) => `3D of ${m}`);
    expect(ok.error, null);
    expect(ok.molBlock3D, `3D of ${molBlock2D}`);
  });

  // A click on the panel's molecule: drawn at twice the shown resolution, a 3x lens under the mouse
  test('dialog-resolution-and-lens', async () => {
    const dlg = await showMoleculeDialog('CC(C)CC(NC(=O)C(N)Cc1ccccc1)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)O');
    try {
      const [source, lens] = Array.from(dlg.root.querySelectorAll('canvas'));
      const shown = source.getBoundingClientRect();
      const ratio = source.width / shown.width / (window.devicePixelRatio || 1);
      expect(ratio >= 1.5, true, `the drawing has ${ratio.toFixed(2)}x the shown resolution`);
      const view = lens.parentElement!;
      const move = new MouseEvent('mousemove',
        {clientX: shown.left + shown.width / 2, clientY: shown.top + shown.height / 2, bubbles: true});
      const lensInked = () => {
        const d = lens.getContext('2d')!.getImageData(0, 0, lens.width, lens.height).data;
        for (let i = 0; i < d.length; i += 4) if (d[i] < 128) return true;
        return false;
      };
      // the drawing lands asynchronously: move again until the lens shows a part of it
      await awaitCheck(() => {
        view.dispatchEvent(move);
        return lens.style.display !== 'none' && lensInked();
      }, 'the lens shows no part of the molecule', 10000, 200);
      view.dispatchEvent(new MouseEvent('mouseleave'));
      expect(lens.style.display, 'none', 'the lens stays after the mouse left');
    } finally {
      dlg.close();
    }
  });

  async function _testToAtomicLevelFunc(
    seqCol: DG.Column<string>, nonlinear: boolean, tgt: TestDataTargetType,
  ): Promise<void> {
    const res = (await sequenceToMolfile(seqCol.dataFrame, seqCol, nonlinear, false, monomerLib, seqHelper))!;
    expect(res.molCol!.semType, DG.SEMTYPE.MOLECULE);
    const resMolStr = res.molCol!.get(0)!;
    expect(!!resMolStr, true, `No molecule generated: ${res.warnings.join('; ')}`);
    // heavy atoms and bonds, read from the molblock itself (no chemistry engine)
    const counts = heavyAtomBondCounts(resMolStr);
    expect(counts.atoms, tgt.atomCount);
    expect(counts.bonds, tgt.bondCount);
  }
});
