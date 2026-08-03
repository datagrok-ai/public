import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, before, category, expect, expectObject, test, timeout} from '@datagrok-libraries/test/src/test';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {generateLongSequence} from '@datagrok-libraries/bio/src/utils/generator';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {initHelmMainPackage} from './utils';
import {getPropertiesDict, SeqPropertiesError} from '../widgets/properties-widget';

const enum MolProps {
  formula = 'formula',
  mw = 'molecular weight',
  extinctCoef = 'extinction coefficient',
}

const enum Tests {
  helm = 'Helm',
  separator = 'separator',
  separator2 = 'separator2', // with a few gaps in a row
  separatorGaps = 'separator-gaps',
  fastaPt = 'fasta-PT',
  fastaPtMsa = 'fasta-PT-MSA',
  fastaRna = 'fasta-RNA',
  fastaDna = 'fasta-DNA',
}

// hwe migration (Phase 7): molecular-weight goldens updated to the standalone
// `@datagrok-libraries/hwe` calculator, which uses IUPAC-2013 standard atomic
// weights (C 12.011, O 15.999, S 32.06, P 30.974) instead of the legacy
// JSDraw2 PT.Lite table (C 12.0107, O 15.9994, S 32.065, P 30.9738). Formulas
// (legacy C,H,N,O,S,P order) and extinction coefficients are unchanged; only
// MW shifts by ≤0.04 Da. Source: hwe src/calc/formula.ts.
const TestsData: {
  [testName: string]: { seq: string, units: NOTATION, separator?: string, alphabet?: ALPHABET, res: {} }
} = {
  [Tests.helm]: {
    seq: 'PEPTIDE1{I.H.A.C.T.[Thr_PO3H2].[Aca].[D-Tyr_Et]}|PEPTIDE2{E.A.[meG]}$PEPTIDE1,PEPTIDE2,4:R3-1:R1$$$V2.0',
    units: NOTATION.HELM,
    res: {[MolProps.formula]: 'C58H94N13O21SP', [MolProps.mw]: 1372.49, [MolProps.extinctCoef]: 1.55}
  },
  [Tests.separator]: {
    seq: 'I/H/A/C/T/Thr_PO3H2/Aca/D-Tyr_Et',
    units: NOTATION.SEPARATOR, separator: '/', alphabet: ALPHABET.UN,
    res: {[MolProps.formula]: 'C47H77N10O15SP', [MolProps.mw]: 1085.22, [MolProps.extinctCoef]: 1.55}
  },
  [Tests.separator2]: {
    seq: 'meI/hHis/Aca/N/T/dK/Thr_PO3H2/Aca/D-Tyr_Et/Aze/dV/E/N////Phe_4Me',
    units: NOTATION.SEPARATOR, separator: '/', alphabet: ALPHABET.UN,
    res: {[MolProps.formula]: 'C91H146N19O25P', [MolProps.mw]: 1937.25, [MolProps.extinctCoef]: 1.49},
  },
  [Tests.separatorGaps]: {
    seq: 'I/H/A//T/Thr_PO3H2/Aca/D-Tyr_Et',
    units: NOTATION.SEPARATOR, separator: '/', alphabet: ALPHABET.UN,
    res: {[MolProps.formula]: 'C44H72N9O14P', [MolProps.mw]: 982.08, [MolProps.extinctCoef]: 1.49},
  },
  [Tests.fastaPt]: {
    seq: 'IHACTTAY',
    units: NOTATION.FASTA, alphabet: ALPHABET.PT,
    res: {[MolProps.formula]: 'C38H58N10O12S', [MolProps.mw]: 879, [MolProps.extinctCoef]: 1.55}
  },
  [Tests.fastaRna]: {
    seq: 'ACGUAUUCC',
    units: NOTATION.FASTA, alphabet: ALPHABET.RNA,
    // hwe migration: EC now uses hwe's nearest-neighbor stacking model
    // (ε = 2·ΣNN − Σinternal-singletons) instead of the legacy Pistoia
    // sum-of-single-base ε (which gave 96.27). T is scored as U.
    res: {[MolProps.formula]: 'C84H107N30O65P9', [MolProps.mw]: 2855.69, [MolProps.extinctCoef]: 86.96}
  },
  [Tests.fastaDna]: {
    seq: 'ACGTATTCC',
    units: NOTATION.FASTA, alphabet: ALPHABET.DNA,
    // hwe NN model with T≈U (was 91.74 under the legacy sum-of-singles model).
    res: {[MolProps.formula]: 'C87H113N30O56P9', [MolProps.mw]: 2753.78, [MolProps.extinctCoef]: 86.96}
  }
};

category('properties-widget', () => {
  let seqHelper: ISeqHelper;
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: any = null;

  before(async () => {
    await initHelmMainPackage();
    seqHelper = await getSeqHelper();

    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  for (const [testName, testData] of Object.entries(TestsData)) {
    test(testName, async () => {
      await testPropertiesDict(testData.seq, testData.units, testData.separator, testData.alphabet, testData.res);
    });
  }

  test('too-long', async () => {
    const colList = generateLongSequence();
    const df = DG.DataFrame.fromColumns(colList);
    const col = df.getCol('MSA');

    const sh = seqHelper.getSeqHandler(col);
    let errCatched = true;
    try {
      const _actPropDict = await getPropertiesDict(col.get(0), sh);
      const k = 11;
    } catch (err: any) {
      if (err instanceof SeqPropertiesError) errCatched = true;
    }
    expect(errCatched, true);
  });

  async function testPropertiesDict(
    seq: string, units: NOTATION, separator: string | undefined, alphabet: ALPHABET | undefined, expPropDict: {}
  ) {
    const col = DG.Column.fromStrings('seq', [seq]) as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = units;
    if (separator) col.setTag(bioTAGS.separator, separator);
    if (alphabet) col.setTag(bioTAGS.alphabet, alphabet);
    const sh = seqHelper.getSeqHandler(col);
    const actPropDict = await getPropertiesDict(seq, sh);
    expectObject(actPropDict, expPropDict);
  }
});
