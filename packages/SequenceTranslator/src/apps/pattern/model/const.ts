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

// todo: remove after refactoring
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

export const STORAGE_NAME: string = 'TestOligoTools1';
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
  export const PATTERN_NAME = 'Pattern';
}

export namespace GRAPH_SETTINGS_KEYS {
  export const IS_ANTISENSE_STRAND_INCLUDED = 'isAntisenseStrandIncluded';
  export const NUCLEOTIDE_SEQUENCES = 'nucleotideSequences';
  export const PHOSPHOROTHIOATE_LINKAGE_FLAGS = 'phosphorothioateLinkageFlags';
  export const STRAND_TERMINUS_MODIFICATIONS = 'strandTerminusModifications';
}

export namespace LEGEND_SETTINGS_KEYS {
  export const PATTERN_NAME = 'patternName';
  export const PATTERN_COMMENT = 'patternComment';
  export const NUCLEOTIDES_WITH_NUMERIC_LABELS = 'nucleotidesWithNumericLabels';
}

export namespace PATTERN_RECORD_KEYS {
  export const PATTERN_CONFIG = 'patternConfig';
  export const AUTHOR_ID = 'authorID';
}

export const GRAPH_SETTINGS_KEY_LIST = [
  GRAPH_SETTINGS_KEYS.IS_ANTISENSE_STRAND_INCLUDED,
  GRAPH_SETTINGS_KEYS.NUCLEOTIDE_SEQUENCES,
  GRAPH_SETTINGS_KEYS.PHOSPHOROTHIOATE_LINKAGE_FLAGS,
  GRAPH_SETTINGS_KEYS.STRAND_TERMINUS_MODIFICATIONS
];

export const LEGEND_SETTINGS_KEY_LIST = [
  LEGEND_SETTINGS_KEYS.PATTERN_NAME,
  LEGEND_SETTINGS_KEYS.PATTERN_COMMENT,
  LEGEND_SETTINGS_KEYS.NUCLEOTIDES_WITH_NUMERIC_LABELS
];

export const PATTERN_RECORD_KEY_LIST = [
  PATTERN_RECORD_KEYS.PATTERN_CONFIG,
  PATTERN_RECORD_KEYS.AUTHOR_ID
];
