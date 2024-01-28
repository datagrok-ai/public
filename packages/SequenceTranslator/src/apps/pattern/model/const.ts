import {StrandType} from './types';

export const enum STRAND {
  SENSE = 'SS',
  ANTISENSE = 'AS',
};
export const STRANDS = [STRAND.SENSE, STRAND.ANTISENSE] as const;
export const STRAND_LABEL: Record<StrandType, string> = {
  [STRAND.SENSE]: 'Sense strand',
  [STRAND.ANTISENSE]: 'Anti sense',
};

export const enum STRAND_END {
  LEFT,
  RIGHT,
};

export const STRAND_ENDS = [STRAND_END.LEFT, STRAND_END.RIGHT] as const;

export const enum TERMINUS {
  FIVE_PRIME = '5\'',
  THREE_PRIME = '3\'',
};
export const TERMINI = [TERMINUS.THREE_PRIME, TERMINUS.FIVE_PRIME] as const;

export const STRAND_TO_END_TERMINUS_MAP = {
  [STRAND.SENSE]: {
    [STRAND_END.LEFT]: TERMINUS.THREE_PRIME,
    [STRAND_END.RIGHT]: TERMINUS.FIVE_PRIME
  },
  [STRAND.ANTISENSE]: {
    [STRAND_END.LEFT]: TERMINUS.FIVE_PRIME,
    [STRAND_END.RIGHT]: TERMINUS.THREE_PRIME
  }
} as const;

export const MAX_SEQUENCE_LENGTH = 35;

export const enum PATTERN_KEY {
  SS_BASES = 'ssBases',
  AS_BASES = 'asBases',
  SS_PTO = 'ssPtoLinkages',
  AS_PTO = 'asPtoLinkages',
  SS_3 = 'ssThreeModification',
  SS_5 = 'ssFiveModification',
  AS_3 = 'asThreeModification',
  AS_5 = 'asFiveModification',
  COMMENT = 'comment',
};

export const USER_STORAGE_KEY: string = 'SequenceTranslator';
export const EXAMPLE_MIN_WIDTH: string = '400px';

// todo: remove after refactoring
export const DEFAULT_PHOSPHOROTHIOATE: boolean = true;

export const OTHER_USERS = 'Other users';

export namespace DEFAULT_PATTERN_CONFIG {
  export const SEQUENCE_LENGTH = 23;
  export const PHOSPHOROTHIOATE = true;
  export const TERMINUS_MODIFICATION = '';
  export const COMMENT = '';
  export const IS_ANTISENSE_STRAND_VISIBLE = true;
  export const PATTERN_NAME = '';
  export const MODIFICATIONS_WITH_NUMERIC_LABELS = [] as string[];
}
