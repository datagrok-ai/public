import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

/** Heavy atoms and bonds of a V3000 molblock, read from its atom and bond blocks (no chemistry engine). */
function heavyAtomBondCounts(molfile: string): { atoms: number, bonds: number } {
  const atomEl = new Map<number, string>();
  let atoms = 0;
  let bonds = 0;
  let block: 'atom' | 'bond' | null = null;
  for (const line of molfile.split('\n')) {
    if (line.startsWith('M  V30 BEGIN ATOM')) block = 'atom';
    else if (line.startsWith('M  V30 BEGIN BOND')) block = 'bond';
    else if (line.startsWith('M  V30 END')) block = null;
    else if (block && line.startsWith('M  V30 ')) {
      const t = line.substring(7).trim().split(/\s+/);
      if (block === 'atom') {
        atomEl.set(parseInt(t[0]), t[1]);
        if (t[1] !== 'H') ++atoms;
      } else if (atomEl.get(parseInt(t[2])) !== 'H' && atomEl.get(parseInt(t[3])) !== 'H')
        ++bonds;
    }
  }
  return {atoms, bonds};
}

type TestDataTargetType = { atomCount: number, bondCount: number };
type TestDataType = {
  src: { seq: string, units: NOTATION },
  tgt: TestDataTargetType,
};

category('toAtomicLevel-ui', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  const tests: { [testName: string]: TestDataType } = {
    'fasta': {
      src: {seq: 'MDYKETLLMPK', units: NOTATION.FASTA},
      tgt: {atomCount: 94, bondCount: 95},
    },
    'fasta-with-gap': {
      src: {seq: 'MD-YKETLLMPK', units: NOTATION.FASTA},
      tgt: {atomCount: 94, bondCount: 95},
    },
    'helm': {
      src: {seq: 'PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2}$$$$', units: NOTATION.HELM},
      tgt: {atomCount: 68, bondCount: 68},
    },
    'helm-with-gap': {
      src: {seq: 'PEPTIDE1{meI.hHis.*.Aca.N.T.dK.Thr_PO3H2}$$$$', units: NOTATION.HELM},
      tgt: {atomCount: 68, bondCount: 68},
    },
  };

  const getDfAndSeqCol = async (testData: TestDataType): Promise<{df: DG.DataFrame, seqCol: DG.Column<string>}> => {
    const seq = testData.src.seq;
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [seq])]);
    await grok.data.detectSemanticTypes(df);
    return {df: df, seqCol: df.getCol('seq')};
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}-nonlinear`, async () => {
      const res = await getDfAndSeqCol(testData);
      await _testToAtomicLevelFunc(res.df, res.seqCol, true, testData.tgt);
    });
  }

  async function _testToAtomicLevelFunc(
    df: DG.DataFrame, seqCol: DG.Column<string>, nonlinear: boolean, tgt: TestDataTargetType,
  ): Promise<void> {
    await grok.functions.call('Bio:toAtomicLevel', {
      table: df, seqCol: seqCol, nonlinear: true, highlight: false
    });
    const molCol = df.col('molfile(seq)');
    expect(molCol?.semType, DG.SEMTYPE.MOLECULE);
    const resMolStr = molCol?.get(0)!;
    expect(!!resMolStr, true, 'No molecule generated');
    const counts = heavyAtomBondCounts(resMolStr);
    expect(counts.atoms, tgt.atomCount);
    expect(counts.bonds, tgt.bondCount);
  }
});
