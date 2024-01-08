import {StrandType, TerminalType} from './types';

export const DEFAULT_PHOSPHOROTHIOATE: boolean = true;
export const DEFAULT_SEQUENCE_LENGTH: number = 23;
export const MAX_SEQUENCE_LENGTH: number = 35;
export const USER_STORAGE_KEY: string = 'SequenceTranslator';
export const EXAMPLE_MIN_WIDTH: string = '400px';

// todo: remove as legacy
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

export const SENSE_STRAND = 'SS' as const;
export const ANTISENSE_STRAND = 'AS' as const;
export const STRANDS = [SENSE_STRAND, ANTISENSE_STRAND];

export const STRAND_LABEL: Record<StrandType, string> = {
  [SENSE_STRAND]: 'Sense strand',
  [ANTISENSE_STRAND]: 'Anti sense',
};

export const FIVE_PRIME = '5\'';
export const THREE_PRIME = '3\'';
export const TERMINI = [FIVE_PRIME, THREE_PRIME] as const;

export const TERMINAL_KEYS: TerminalType[] = [THREE_PRIME, FIVE_PRIME];
export const TERMINAL: Record<TerminalType, number> = {
  [THREE_PRIME]: 3,
  [FIVE_PRIME]: 5,
};

export const OTHER_USERS = 'Other users';
