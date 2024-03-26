import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectObject, test} from '@datagrok-libraries/utils/src/test';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

import {getPropertiesDict} from '../widgets/properties-widget';

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
  fastaDna = 'fasta-DNA',
}

const TestsData: {
  [testName: string]: { seq: string, units: NOTATION, separator?: string, alphabet?: ALPHABET, res: {} }
} = {
  [Tests.helm]: {
    seq: 'PEPTIDE1{I.H.A.C.T.[Thr_PO3H2].[Aca].[D-Tyr_Et]}|PEPTIDE2{E.A.[meG]}$PEPTIDE1,PEPTIDE2,4:R3-1:R1$$$V2.0',
    units: NOTATION.HELM,
    res: {[MolProps.formula]: 'C58H94N13O21SP', [MolProps.mw]: 1372.48, [MolProps.extinctCoef]: 1.55}
  },
  [Tests.separator]: {
    seq: 'I/H/A/C/T/Thr_PO3H2/Aca/D-Tyr_Et',
    units: NOTATION.SEPARATOR, separator: '/', alphabet: ALPHABET.UN,
    res: {[MolProps.formula]: 'C47H77N10O15SP', [MolProps.mw]: 1085.21, [MolProps.extinctCoef]: 1.55}
  },
  [Tests.separator2]: {
    seq: 'meI/hHis/Aca/N/T/dK/Thr_PO3H2/Aca/D-Tyr_Et/Aze/dV/E/N////Phe_4Me',
    units: NOTATION.SEPARATOR, separator: '/', alphabet: ALPHABET.UN,
    res: {[MolProps.formula]: 'C91H146N19O25P', [MolProps.mw]: 1937.21, [MolProps.extinctCoef]: 1.49},
  },
  [Tests.separatorGaps]: {
    seq: 'I/H/A//T/Thr_PO3H2/Aca/D-Tyr_Et',
    units: NOTATION.SEPARATOR, separator: '/', alphabet: ALPHABET.UN,
    res: {[MolProps.formula]: 'C44H72N9O14P', [MolProps.mw]: 982.07, [MolProps.extinctCoef]: 1.49},
  },
  [Tests.fastaPt]: {
    seq: 'IHACTTAY',
    units: NOTATION.FASTA, alphabet: ALPHABET.PT,
    res: {[MolProps.formula]: 'C38H58N10O12S', [MolProps.mw]: 878.99, [MolProps.extinctCoef]: 1.55}
  },
  [Tests.fastaDna]: {
    seq: 'ACGTATTCC',
    units: NOTATION.FASTA, alphabet: ALPHABET.DNA,
    res: {[MolProps.formula]: 'C87H113N30O56P9', [MolProps.mw]: 2753.76, [MolProps.extinctCoef]: 91.74}
  }
};

category('properties-widget', () => {
  for (const [testName, testData] of Object.entries(TestsData)) {
    test(testName, async () => {
      testPropertiesDict(testData.seq, testData.units, testData.separator, testData.alphabet, testData.res);
    });
  }
});

function testPropertiesDict(
  seq: string, units: NOTATION, separator: string | undefined, alphabet: ALPHABET | undefined, expPropDict: {}
) {
  const col = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [seq]);
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.setTag(DG.TAGS.UNITS, units);
  if (separator) col.setTag(bioTAGS.separator, separator);
  if (alphabet) col.setTag(bioTAGS.alphabet, alphabet);
  const uh = UnitsHandler.getOrCreate(col);
  const actPropDict = getPropertiesDict(seq, uh);
  expectObject(actPropDict, expPropDict);
}
