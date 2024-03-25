import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export {
  SplitterFunc,

} from './types';

export {
  TAGS,
  NOTATION,
  ALPHABET,
  ALIGNMENT,
  candidateAlphabets,
  positionSeparator
} from './consts';

export {
  splitterAsFasta,
  getSplitterWithSeparator,
  splitterAsHelm,
  getAlphabet,
  getAlphabetSimilarity,
  monomerToShort,
  pickUpPalette,
  pickUpSeqCol,
  getPaletteByType,
  MonomerToShortFunc
} from './utils';
