import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {sequenceValidation} from '../validation';
import {SYNTHESIZERS, TECHNOLOGIES} from '../constants';

const ADDITIONAL_CODES: string[] = [];
const inputOutputObjs: {
  [inputSequence: string]: {
    indexOfFirstNotValidChar: number,
    synthesizer: string[] | null,
    technology: string[] | null
  }
} = {
  'fNAmCmGmAmCpsmU': {
    indexOfFirstNotValidChar: 0,
    synthesizer: null,
    technology: null,
  },
  // 'moeTmoeAmoe5mCfUfA': {
  //   indexOfFirstNotValidChar: 14,
  //   synthesizer: null,
  //   technology: null,
  // },
  // 'rArCrGRUrArCrG': {
  //   indexOfFirstNotValidChar: 6,
  //   synthesizer: null,
  //   technology: null,
  // },
  // 'mGrArCrU': {
  //   indexOfFirstNotValidChar: -1,
  //   synthesizer: [SYNTHESIZERS.GCRS],
  //   technology: null,
  // },
  // 'fUfA5mC': {
  //   indexOfFirstNotValidChar: -1,
  //   synthesizer: [SYNTHESIZERS.GCRS],
  //   technology: null,
  // },
  'mUmUmAmUmCmUmUmGmAmUmUmG': {
    indexOfFirstNotValidChar: -1,
    synthesizer: [SYNTHESIZERS.GCRS],
    technology: [TECHNOLOGIES.SI_RNA],
  },
  'mApsmApsfGmAmUmCfAfCfGmG': {
    indexOfFirstNotValidChar: -1,
    synthesizer: [SYNTHESIZERS.GCRS],
    technology: [TECHNOLOGIES.SI_RNA],
  },
  // 'ArA': {
  //   indexOfFirstNotValidChar: 1,
  //   synthesizer: null,
  //   technology: null,
  // },
  // 'AdA': {
  //   indexOfFirstNotValidChar: 1,
  //   synthesizer: null,
  //   technology: null,
  // },
  // 'rAA': {
  //   indexOfFirstNotValidChar: 2,
  //   synthesizer: null,
  //   technology: null,
  // },
  // 'dAA': {
  //   indexOfFirstNotValidChar: 2,
  //   synthesizer: null,
  //   technology: null,
  // },
  'fAGACGT': {
    indexOfFirstNotValidChar: -1,
    synthesizer: [SYNTHESIZERS.GCRS],
    technology: [TECHNOLOGIES.SI_RNA],
  },
  // 'AGfA': {
  //   indexOfFirstNotValidChar: -1,
  //   synthesizer: [SYNTHESIZERS.GCRS],
  //   technology: null,
  // },
  // 'TrA': {
  //   indexOfFirstNotValidChar: -1,
  //   synthesizer: [SYNTHESIZERS.GCRS],
  //   technology: [TECHNOLOGIES.ASO_GAPMERS],
  // },
  // 'rAT': {
  //   indexOfFirstNotValidChar: -1,
  //   synthesizer: [SYNTHESIZERS.GCRS],
  //   technology: [TECHNOLOGIES.ASO_GAPMERS],
  // },
  'mUrA': {
    indexOfFirstNotValidChar: -1,
    synthesizer: [SYNTHESIZERS.GCRS],
    technology: [TECHNOLOGIES.SI_RNA],
  },
  'fAGgfA': {
    indexOfFirstNotValidChar: 3,
    synthesizer: [SYNTHESIZERS.GCRS],
    technology: null,
  },
};

category('oligo-batch-calculator-sequence-validation', () => {
  for (const [inputSequence, outputObj] of Object.entries(inputOutputObjs)) {
    test(inputSequence, async () => expect(
      JSON.stringify(sequenceValidation(inputSequence, ADDITIONAL_CODES)),
      JSON.stringify(outputObj)),
    );
  }
});
